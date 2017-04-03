#pragma once

/*
 * Given an underlying sequence S and an index map I,
 * present an indirect view V of S such that
 * V[i] = S[I[i]] for all valid indices i.
 *
 * Implemented as a transform view of the index map.
 */

#include <utility>

#include <util/deduce_return.hpp>
#include <util/transform.hpp>
#include <util/meta.hpp>

namespace nest {
namespace mc {
namespace util {

// Seq: random access sequence

namespace impl {
    template <typename Data>
    struct indirect_accessor {
        using reference = typename util::sequence_traits<Data>::reference;

        Data data;
        template <typename X>
        indirect_accessor(X&& data): data(std::forward<X>(data)) {}

        template <typename I>
        reference operator()(const I& i) const { return data[i]; }
    };
}

template <typename RASeq, typename Seq>
auto indirect_view(RASeq& data, const Seq& index_map)
DEDUCED_RETURN_TYPE(transform_view(index_map, impl::indirect_accessor<RASeq&>(data)));

template <typename RASeq, typename Seq>
auto indirect_view(RASeq&& data, const Seq& index_map)
DEDUCED_RETURN_TYPE(transform_view(index_map, impl::indirect_accessor<RASeq>(std::move(data))));

} // namespace util
} // namespace mc
} // namespace nest
