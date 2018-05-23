#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include <epoch.hpp>

namespace arb {

template <typename T>
class epoch_buffer {
public:
    using index_type = epoch::index_type;
    using value_type = T;
    using iterator = typename std::array<T,2>::iterator;
    using const_iterator = typename std::array<T,2>::const_iterator;

    epoch_buffer() = default;

    epoch_buffer(epoch_buffer&&) = default;
    epoch_buffer(const epoch_buffer&) = default;
    epoch_buffer& operator=(const epoch_buffer&) = default;
    epoch_buffer& operator=(epoch_buffer&&) = default;

    value_type& operator[](index_type i) {
        return buffers_[index(i)];
    }

    const value_type& operator[](index_type i) const {
        return buffers_[index(i)];
    }

    iterator begin() { return buffers_.begin(); }
    iterator end()   { return buffers_.end(); }
    const_iterator begin() const { return buffers_.begin(); }
    const_iterator end()   const { return buffers_.end(); }

private:
    std::array<T, 2> buffers_;
    int index(index_type i) const {
        return i&1;
    }
};

template <typename T>
struct epoch_vector {
    using store_type = epoch_buffer<std::vector<T>>;
    using index_type = typename store_type::index_type;
    store_type store_;
    using iterator = typename store_type::iterator;
    using const_iterator = typename store_type::iterator;

    using value_type = T;

    epoch_vector(std::size_t n) {
        store_[0] = std::vector<T>(n);
        store_[1] = std::vector<T>(n);
    }

    std::vector<T>& operator[](index_type i) {
        return store_[i];
    }
    const std::vector<T>& operator[](index_type i) const {
        return store_[i];
    }
    std::vector<T>& operator[](epoch e) {
        return store_[e.id];
    }
    const std::vector<T>& operator[](epoch e) const {
        return store_[e.id];
    }
    iterator begin() { return store_.begin(); }
    iterator end()   { return store_.end(); }
    const_iterator begin() const { return store_.begin(); }
    const_iterator end()   const { return store_.end(); }
};

} // namespace arb

