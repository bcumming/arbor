#pragma once

// 'Place' morphology in 3-d by applying an isometry to
// sample points and interpolating linearly.

#include <cmath>
#include <limits>
#include <utility>

#include <arbor/morph/morphology.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/math.hpp>

namespace arb {

struct place_pwlin_data;

struct place_pwlin {
    explicit place_pwlin(const morphology& m, const isometry& iso = isometry{});

    // Any point corresponding to the location loc.
    mpoint at(mlocation loc) const;

    // All points corresponding to the location loc.
    std::vector<mpoint> all_at(mlocation loc) const;

    // A minimal set of segments or part segments whose union is coterminous with extent.
    std::vector<msegment> segments(const mextent& extent) const;

    // Maximal set of segments or part segments whose union is coterminous with extent.
    std::vector<msegment> all_segments(const mextent& extent) const;

    // The closest location to p. Returns the location and its distance from the input coordinates.
    std::pair<mlocation, double> closest(double x, double y, double z) const;

private:
    std::shared_ptr<place_pwlin_data> data_;
};

} // namespace arb


