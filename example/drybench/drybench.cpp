/*
 * A miniapp that demonstrates how to use dry_run mode
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include <arbor/assert_macro.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/common_types.hpp>
#include <arbor/context.hpp>
#include <arbor/benchmark_cell.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/symmetric_recipe.hpp>
#include <arbor/recipe.hpp>
#include <arbor/version.hpp>

#include <arborenv/concurrency.hpp>

#include <sup/ioutil.hpp>
#include <sup/json_meter.hpp>
#include <sup/json_params.hpp>

struct run_params {
    std::string name = "default";
    unsigned num_cells_per_rank = 200;
    unsigned num_ranks = 1000;
    double min_delay = 10;        // ms
    double duration = 1000;       // ms
    double spike_frequency = 100.; // Hz
    double realtime_ratio = 0.01;
    unsigned num_synapses = 1000;
};

void print_params(run_params p) {
    const auto nr = p.num_ranks;
    const auto nc = p.num_cells_per_rank;
    const auto ts = p.duration/1000;  // time in seconds
    const auto ns_per_cell = p.spike_frequency*ts;
    std::cout << "================ params ================\n";
    std::cout << "expected thread time : " << nc*p.realtime_ratio*ts << "\n";
    std::cout << "expected spikes/cell : " << ns_per_cell << "\n";
    std::cout << "expected spikes total: " << nc*ns_per_cell*nr << "\n";
    std::cout << "========================================\n\n";
};

run_params read_options(int argc, char** argv);

using arb::cell_gid_type;
using arb::cell_lid_type;
using arb::cell_size_type;
using arb::cell_member_type;
using arb::cell_kind;
using arb::time_type;

class tile_desc: public arb::tile {
public:
    tile_desc(run_params params):
            params_(params),
            num_cells_(params.num_cells_per_rank),
            num_tiles_(params.num_ranks)
    {}
        //auto tile = std::make_unique<tile_desc>(params.num_cells_per_rank,
                //params.num_ranks, params.cell, params.min_delay);

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    cell_size_type num_tiles() const override {
        return num_tiles_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        using RNG = std::mt19937_64;
        auto gen = arb::poisson_schedule(params_.spike_frequency/1000, RNG(gid));
        return arb::benchmark_cell{std::move(gen), params_.realtime_ratio};
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::benchmark;
    }

    arb::util::any get_global_properties(arb::cell_kind) const override {
        arb::cable_cell_global_properties gprop;
        gprop.default_parameters = arb::neuron_parameter_defaults;
        return gprop;
    }

    // Each cell has one spike detector (at the soma).
    cell_size_type num_sources(cell_gid_type gid) const override {
        return 1;
    }

    // The cell has one target synapse
    cell_size_type num_targets(cell_gid_type gid) const override {
        return params_.num_synapses;
    }

    // Each cell has num_synapses incoming connections, from any cell in the
    // network spanning all ranks, src gid in {0, ..., num_cells_*num_tiles_ - 1}.
    std::vector<arb::cell_connection> connections_on(cell_gid_type gid) const override {
        std::uniform_int_distribution<cell_gid_type>
            source_distribution(0, num_cells_*num_tiles_ - 2);

        std::vector<arb::cell_connection> conns;
        auto src_gen = std::mt19937(gid);
        for (unsigned i=0; i<params_.num_synapses; ++i) {
            auto src = source_distribution(src_gen);
            if (src>=gid) ++src;
            conns.push_back(arb::cell_connection({src, 0}, {gid, 0}, 0, params_.min_delay));
        }

        return conns;
        //return {arb::cell_connection({src, 0}, {gid, 0}, 0, params_.min_delay)};
    }

    // Return an event generator on every 20th gid. This function needs to generate events
    // for ALL cells on ALL ranks. This is because the symmetric recipe can not easily
    // translate the src gid of an event generator
    std::vector<arb::event_generator> event_generators(cell_gid_type gid) const override {
        std::vector<arb::event_generator> gens;
        if (gid%20 == 0) {
            gens.push_back(arb::explicit_generator(arb::pse_vector{{{gid, 0}, 0.1, 1.0}}));
        }
        return gens;
    }

    // There is one probe (for measuring voltage at the soma) on the cell.
    cell_size_type num_probes(cell_gid_type gid)  const override {
        return 0;
    }

private:
    run_params params_;
    cell_size_type num_cells_;
    cell_size_type num_tiles_;
};

int main(int argc, char** argv) {
    try {
        bool root = true;
        auto params = read_options(argc, argv);

        print_params(params);

        auto resources = arb::proc_allocation();
        if (auto nt = arbenv::get_env_num_threads()) {
            resources.num_threads = nt;
        }
        else {
            resources.num_threads = arbenv::thread_concurrency();
        }
        auto ctx = arb::make_context(resources);

        ctx = arb::make_context(resources, arb::dry_run_info(params.num_ranks, params.num_cells_per_rank));
        arb_assert(arb::num_ranks(ctx)==params.num_ranks);

        arb::profile::profiler_initialize(ctx);
        std::cout << sup::mask_stream(root);

        // Print a banner with information about hardware configuration
        std::cout << "threads:  " << num_threads(ctx) << "\n";
        std::cout << "ranks:    " << num_ranks(ctx) << "\n" << std::endl;

        arb::profile::meter_manager meters;
        meters.start(ctx);

        // Create an instance of our tile and use it to make a symmetric_recipe.
        auto tile = std::make_unique<tile_desc>(params);
        arb::symmetric_recipe recipe(std::move(tile));

        auto decomp = arb::partition_load_balance(recipe, ctx);

        // Construct the model.
        arb::simulation sim(recipe, decomp, ctx);

        meters.checkpoint("model-init", ctx);

        // Run the simulation for 100 ms, with time steps of 0.025 ms.
        sim.run(params.duration, 0.025);

        meters.checkpoint("model-run", ctx);

        auto ns = sim.num_spikes();
        auto total_cells = params.num_ranks*params.num_cells_per_rank;
        std::cout << "\n" << ns << " spikes generated at rate of "
                  << ns/total_cells << " spikes per cell\n\n";

        auto profile = arb::profile::profiler_summary();
        std::cout << profile << "\n";

        auto report = arb::profile::make_meter_report(meters, ctx);
        std::cout << report;
    }
    catch (std::exception& e) {
        std::cerr << "exception caught in benchmark: \n" << e.what() << "\n";
        return 1;
    }

    return 0;
}

run_params read_options(int argc, char** argv) {
    using sup::param_from_json;

    run_params params;
    if (argc<2) {
        std::cout << "Using default parameters.\n";
        return params;
    }
    if (argc>2) {
        throw std::runtime_error("More than command line one option not permitted.");
    }

    std::string fname = argv[1];
    std::cout << "Loading parameters from file: " << fname << "\n";
    std::ifstream f(fname);

    if (!f.good()) {
        throw std::runtime_error("Unable to open input parameter file: "+fname);
    }

    nlohmann::json json;
    json << f;

    param_from_json(params.name, "name", json);
    param_from_json(params.num_cells_per_rank, "num-cells-per-rank", json);
    param_from_json(params.num_ranks, "num-ranks", json);
    param_from_json(params.num_synapses, "num-synapses", json);
    param_from_json(params.duration, "duration", json);
    param_from_json(params.min_delay, "min-delay", json);
    param_from_json(params.spike_frequency, "spike-frequency", json);
    param_from_json(params.realtime_ratio, "realtime-ratio", json);

    if (!json.empty()) {
        for (auto it=json.begin(); it!=json.end(); ++it) {
            std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
        }
        std::cout << "\n";
    }

    return params;
}
