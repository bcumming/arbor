#pragma once

#include <iostream>
#include <vector>
#include <mutex>

#include <common_types.hpp>
#include <spike.hpp>
#include <threading/threading.hpp>

namespace arb {

namespace impl {
    template <typename T>
    struct iterable_store {
        using size_type = std::size_t;

        using store_type = std::unordered_map<size_type, T>;
        using store_iterator = typename store_type::iterator;
        using const_store_iterator = typename store_type::const_iterator;

        struct iterator: public store_iterator {
            using base = store_iterator;
            iterator(): base() {}
            iterator(store_iterator it): base(it) {}

            T* operator->() {
                return static_cast<T*>(&(store_iterator::operator->()->second));
            }
            T& operator*(){
                return store_iterator::operator*().second;
            }
        };

        struct const_iterator: public const_store_iterator {
            using base = const_store_iterator;
            const_iterator(): base() {}
            const_iterator(base it): base(it) {}
            const_iterator(store_iterator it): base(it) {}

            const T* operator->() {
                return static_cast<const T*>(&(const_store_iterator::operator->()->second));
            }
            const T& operator*(){
                return const_store_iterator::operator*().second;
            }
        };

        std::unordered_map<size_type, T> store_;

        T& get() {
            std::lock_guard<std::mutex> g(mtx);
            const auto id = arb::threading::thread_id();
            return store_[id];
        }

        size_type size() const {
            return store_.size();
        }

        iterator begin() {
            return store_.begin();
        }
        iterator end() {
            return store_.end();
        }

        const_iterator begin() const {
            return store_.begin();
        }
        const_iterator end() const {
            return store_.end();
        }

        std::mutex mtx;
    };
} // namespace impl


/// Handles the complexity of managing thread private buffers of spikes.
/// Internally stores one thread private buffer of spikes for each hardware thread.
/// This can be accessed directly using the get() method, which returns a reference to
/// The thread private buffer of the calling thread.
/// The insert() and gather() methods add a vector of spikes to the buffer,
/// and collate all of the buffers into a single vector respectively.
class thread_private_spike_store {
private :
    /// thread private storage for accumulating spikes
    using local_spike_store_type = impl::iterable_store<std::vector<spike>>;

    local_spike_store_type buffers_;

public :
    using iterator = typename local_spike_store_type::iterator;
    using const_iterator = typename local_spike_store_type::const_iterator;

    /// Collate all of the individual buffers into a single vector of spikes.
    /// Does not modify the buffer contents.
    std::vector<spike> gather() const;

    /// Clear all of the thread private buffers
    void clear();

    /// Append the passed spikes to the end of the thread private buffer of the
    /// calling thread
    void insert(const std::vector<spike>& spikes);

    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
};

} // namespace arb
