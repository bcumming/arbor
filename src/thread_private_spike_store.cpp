#include <unordered_map>

#include <thread_private_spike_store.hpp>

namespace arb {

std::vector<spike> thread_private_spike_store::gather() const {
    std::vector<spike> spikes;
    unsigned num_spikes = 0u;
    for (auto& b : buffers_) {
        num_spikes += b.size();
    }
    spikes.reserve(num_spikes);

    for (auto& b : buffers_) {
        spikes.insert(spikes.begin(), b.begin(), b.end());
    }

    return spikes;
}

void thread_private_spike_store::clear() {
    for (auto& b: buffers_) {
        b.clear();
    }
}

void thread_private_spike_store::insert(const std::vector<spike>& spikes) {
    auto& buff = buffers_.get();
    buff.insert(buff.end(), spikes.begin(), spikes.end());
}

thread_private_spike_store::iterator
thread_private_spike_store::begin() { return buffers_.begin(); }

thread_private_spike_store::iterator
thread_private_spike_store::end() { return buffers_.begin(); }

thread_private_spike_store::const_iterator
thread_private_spike_store::begin() const { return buffers_.begin(); }

thread_private_spike_store::const_iterator
thread_private_spike_store::end() const { return buffers_.begin(); }

} // namespace arb


