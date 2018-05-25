#include <set>
#include <vector>

#include <backends.hpp>
#include <cell_group.hpp>
#include <cell_group_factory.hpp>
#include <domain_decomposition.hpp>
#include <merge_events.hpp>
#include <simulation.hpp>
#include <recipe.hpp>
#include <util/filter.hpp>
#include <util/span.hpp>
#include <util/unique_any.hpp>
#include <profiling/profiler.hpp>

#include <hpx/include/apply.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/hpx_start.hpp>
#include <hpx/hpx_suspend.hpp>
#include <hpx/include/parallel_for_loop.hpp>

namespace arb {

template <typename T, typename U>
inline bool atomic_dec_test_update(std::atomic<T>& c, U def) {
    const bool zero = !(--c);
    if (zero) {
        c = def;
    }
    return zero;
}

simulation::simulation(const recipe& rec, const domain_decomposition& decomp):
    communicator_(rec, decomp),
    local_spikes_(decomp.groups.size()),
    event_lanes_(communicator_.num_local_cells()),
    cell_task_counter_(decomp.groups.size())
{
    hpx::resume();
    hpx::apply([&]() {
        const auto num_local_cells = communicator_.num_local_cells();

        // Cache the minimum delay of the network
        min_delay_ = communicator_.min_delay();

        // Initialize empty buffers for pending events for each local cell
        pending_events_.resize(num_local_cells);

        event_generators_.resize(num_local_cells);
        cell_local_size_type lidx = 0;
        const auto& grps = decomp.groups;
        for (auto i: util::make_span(0, grps.size())) {
            for (auto gid: grps[i].gids) {
                // Store mapping of gid to local cell index.
                gid_to_local_[gid] = lidx;

                // Set up the event generators for cell gid.
                auto rec_gens = rec.event_generators(gid);
                auto& gens = event_generators_[lidx];
                if (rec_gens.size()) {
                    // Allocate two empty event generators that will be used to
                    // merge events from the communicator and those already queued
                    // for delivery in future epochs.
                    gens.reserve(2+rec_gens.size());
                    gens.resize(2);
                    for (auto& g: rec_gens) {
                        gens.push_back(std::move(g));
                    }
                }
                ++lidx;
            }
        }

        // Generate the cell groups in parallel, with one task per cell group.
        cell_groups_.resize(decomp.groups.size());
        hpx::parallel::for_loop(hpx::parallel::execution::par, 0, cell_groups_.size(),
            [&](cell_gid_type i) {
                cell_groups_[i] = cell_group_factory(rec, decomp.groups[i]);
            });
    });
    hpx::suspend();
}

void simulation::reset() {
    t_ = 0.;

    // Reset cell group state.
    for (auto& group: cell_groups_) {
        group->reset();
    }

    // Clear all pending events in the event lanes.
    for (auto& lanes: event_lanes_) {
        for (auto& lane: lanes) {
            lane.clear();
        }
    }

    // Reset all event generators, and advance to t_.
    for (auto& lane: event_generators_) {
        for (auto& gen: lane) {
            gen.reset();
            gen.advance(t_);
        }
    }

    for (auto& lane: pending_events_) {
        lane.clear();
    }

    communicator_.reset();

    for (auto& lanes: local_spikes_) {
        for (auto& lane: lanes) {
            lane.clear();
        }
    }
}

epoch simulation::merge_lane(unsigned i, epoch ep) {
    //std::cout << " E  : " << i << " @ " << ep << " (event_lanes[" << ep.id-1 << "]) -> event_lanes[" << ep.id << "]\n";
    merge_events(
        ep.t0(), ep.t1(),
        event_lanes_[ep.id-1][i],   // in:  the current event lane
        pending_events_[i],         // in:  events from the communicator
        event_generators_[i],       // in:  event generators for this lane
        event_lanes_[ep.id][i]);    // out: the event lane for the next epoch
    pending_events_[i].clear();

    hpx::async(&simulation::launch_cell_update, this, i, ep);

    return ep;
}

epoch simulation::merge_lanes(epoch ep) {
    const unsigned n = communicator_.num_local_cells();
    for (unsigned i=0; i<n; ++i) {
        hpx::async(&simulation::merge_lane, this, i, ep);
    }

    // launch next epochs event merge tasks if the exchange has finished
    ++ep;
    /*
    if (atomic_dec_test_update(merge_task_counter_[ep.id], 2)) {
        std::cout << "launching " << ep << " from merge_lanes\n";
        merge_lanes(ep);
    }
    */

    return ep;
}

// task that performs spike exchange with the spikes generated in
// the previous integration period, generating the postsynaptic
// events that must be delivered at the start of the next
// integration period at the latest.
epoch simulation::exchange(epoch ep) {
    //std::cout << "exchange at " << ep << " (local_spikes[" << ep.id << "])\n";
    PE(communication_exchange_gatherlocal);

    /*******************************/
    std::size_t nlocal_spikes = 0;
    auto& input_spikes = local_spikes_[ep.id];
    for (auto& l: input_spikes) {
        nlocal_spikes += l.size();
    }
    std::vector<spike> local_spikes;
    local_spikes.reserve(nlocal_spikes);
    for (const auto& l: input_spikes) {
        local_spikes.insert(local_spikes.end(), l.begin(), l.end());
    }
    /*******************************/
    PL();

    auto global_spikes = communicator_.exchange(local_spikes);

    PE(communication_spikeio);
    local_export_callback_(local_spikes);
    global_export_callback_(global_spikes.values());
    PL();

    PE(communication_walkspikes);
    communicator_.make_event_queues(global_spikes, pending_events_);
    PL();

    // generate events to be delivered in 2 epoch's time.
    auto epm = ep + 2;
    //if (!epm.finished() && atomic_dec_test_update(merge_task_counter_[epm.id], 2)) {
    if (!epm.finished()) {
        //std::cout << "launching " << epm << "\n";
        merge_lanes(epm);
    }
    return ep;
};

epoch simulation::update_cell(unsigned i, epoch ep) {
    //std::cout << " C  : " << i << " @ " << ep << " (event_lanes[" << ep.id << "]) -> local_spikes[" << ep.id << "]\n";
    auto &group = cell_groups_[i];

    auto queues = util::subrange_view(
            event_lanes_[ep],
            communicator_.group_queue_range(i));
    group->advance(ep, dt_, queues);
    PE(advance_spikes);

    auto& spikes = group->spikes();
    auto& buffer = local_spikes_[ep][i];
    buffer.clear();
    buffer.insert(buffer.end(), spikes.begin(), spikes.end());

    group->clear_spikes();
    PL();

    //
    // WARNING: the exchange must be performed before launching
    // the following cell to avoid the situation where two exchanges
    // are performed simultaneously
    //

    // HERE: decrement exchange_task_counter_
    auto& c = exchange_task_counter_[ep.id];
    if (atomic_dec_test_update(c, num_groups())) {
        exchange(ep);
    }

    // HERE: update cell_task_counter_ and launch if needed
    auto f = hpx::async(&simulation::launch_cell_update, this, i, ep+1);

    return ep+1;
}

// launch update for cell i for epoch ep
epoch simulation::launch_cell_update(unsigned i, epoch ep) {
        auto& c = cell_task_counter_[ep.id-1][i];
        if (atomic_dec_test_update(c, 2) && !ep.finished()) {
            hpx::async(&simulation::update_cell, this, i, ep);
        }
        return ep;
}

time_type simulation::run(time_type t_final, time_type dt) {
    dt_ = dt;

    hpx::resume();
    hpx::apply([&]() {
        // Calculate the size of the largest possible time integration interval
        // before communication of spikes is required.
        // If spike exchange and cell update are serialized, this is the
        // minimum delay of the network, however we use half this period
        // to overlap communication and computation.
        const time_type t_interval = min_delay_/2;

        // set up the counters
        util::fill(exchange_task_counter_, num_groups());
        util::fill(merge_task_counter_, 2);
        for (auto& c: cell_task_counter_) {
            util::fill(c, 2);
        }
        util::fill(cell_task_counter_[-1], 1);

        auto e = epoch(0, t_, t_interval, t_final);
        merge_lanes(e);
        merge_lanes(e+1);
    });
    hpx::suspend();

    return t_;
}

// Populate the event lanes for epoch+1 (i.e event_lanes_[epoch+1)]
// Update each lane in parallel, if supported by the threading backend.
// On completion event_lanes[epoch+1] will contain sorted lists of events with
// delivery times due in or after epoch+1. The events will be taken from the
// following sources:
//      event_lanes[epoch]: take all events ≥ t_from
//      event_generators  : take all events < t_to
//      pending_events    : take all events
//void simulation::setup_events(time_type t_from, time_type t_to, std::size_t epoch) {

sampler_association_handle simulation::add_sampler(
        cell_member_predicate probe_ids,
        schedule sched,
        sampler_function f,
        sampling_policy policy)
{
    sampler_association_handle h = sassoc_handles_.acquire();

    threading::parallel_for::apply(0, cell_groups_.size(),
        [&](std::size_t i) {
            cell_groups_[i]->add_sampler(h, probe_ids, sched, f, policy);
        });

    return h;
}

void simulation::remove_sampler(sampler_association_handle h) {
    threading::parallel_for::apply(0, cell_groups_.size(),
        [&](std::size_t i) {
            cell_groups_[i]->remove_sampler(h);
        });

    sassoc_handles_.release(h);
}

void simulation::remove_all_samplers() {
    threading::parallel_for::apply(0, cell_groups_.size(),
        [&](std::size_t i) {
            cell_groups_[i]->remove_all_samplers();
        });

    sassoc_handles_.clear();
}

std::size_t simulation::num_spikes() const {
    return communicator_.num_spikes();
}

std::size_t simulation::num_groups() const {
    return cell_groups_.size();
}

void simulation::set_binning_policy(binning_kind policy, time_type bin_interval) {
    for (auto& group: cell_groups_) {
        group->set_binning_policy(policy, bin_interval);
    }
}

void simulation::set_global_spike_callback(spike_export_function export_callback) {
    global_export_callback_ = std::move(export_callback);
}

void simulation::set_local_spike_callback(spike_export_function export_callback) {
    local_export_callback_ = std::move(export_callback);
}

util::optional<cell_size_type> simulation::local_cell_index(cell_gid_type gid) {
    auto it = gid_to_local_.find(gid);
    return it==gid_to_local_.end()?
        util::nullopt:
        util::optional<cell_size_type>(it->second);
}

void simulation::inject_events(const pse_vector& events) {
    // Push all events that are to be delivered to local cells into the
    // pending event list for the event's target cell.
    for (auto& e: events) {
        if (e.time<t_) {
            throw std::runtime_error(
                "simulation::inject_events(): attempt to inject an event at time: "
                + std::to_string(e.time)
                + " ms, which is earlier than the current simulation time: "
                + std::to_string(t_)
                + " ms. Events must be injected on or after the current simulation time.");
        }
        // local_cell_index returns an optional type that evaluates
        // to true iff the gid is a local cell.
        if (auto lidx = local_cell_index(e.target.gid)) {
            pending_events_[*lidx].push_back(e);
        }
    }
}

} // namespace arb
