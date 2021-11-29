#pragma once

#include <algorithm>
#include <cstdlib>
#include <ostream>
#include <vector>

#include <arbor/math.hpp>
#include <arbor/util/hash_def.hpp>
#include <arbor/util/lexcmp_def.hpp>

//
//  Types used to identify concrete locations.
//

namespace arb {

using msize_t = std::uint32_t;
constexpr msize_t mnpos = msize_t(-1);

// a morphology sample point: a 3D location and radius.
struct mpoint {
    double x, y, z;  // [µm]
    double radius;   // [μm]

    friend bool operator==(const mpoint& l, const mpoint& r);
    friend std::ostream& operator<<(std::ostream&, const mpoint&);
    friend bool operator==(const mpoint& a, const mpoint& b) {
        return a.x==b.x && a.y==b.y && a.z==b.z && a.radius==b.radius;
    }
    friend bool operator!=(const mpoint& a, const mpoint& b) {
        return !(a==b);
    }
};

mpoint lerp(const mpoint& a, const mpoint& b, double u);
bool is_collocated(const mpoint& a, const mpoint& b);
double distance(const mpoint& a, const mpoint& b);

// Indicate allowed comparison operations for classifying regions
enum class comp_op {
    lt,
    le,
    gt,
    ge
};

// Describe a cable segment between two adjacent samples.
struct msegment {
    msize_t id;
    mpoint prox;
    mpoint dist;
    int tag;

    friend std::ostream& operator<<(std::ostream&, const msegment&);
};


// Describe a specific location on a morpholology.
struct mlocation {
    // The id of the branch.
    msize_t branch;
    // The relative position on the branch ∈ [0,1].
    double pos;

    friend std::ostream& operator<<(std::ostream&, const mlocation&);
};

// branch ≠ npos and 0 ≤ pos ≤ 1
bool test_invariants(const mlocation&);
ARB_DEFINE_LEXICOGRAPHIC_ORDERING(mlocation, (a.branch,a.pos), (b.branch,b.pos));

using mlocation_list = std::vector<mlocation>;
std::ostream& operator<<(std::ostream& o, const mlocation_list& l);

// Tests whether each location in the list satisfies the invariants for a location,
// and that the locations in the vector are ordered.
bool test_invariants(const mlocation_list&);

// Multiset operations on location lists.
mlocation_list sum(const mlocation_list&, const mlocation_list&);
mlocation_list join(const mlocation_list&, const mlocation_list&);
mlocation_list intersection(const mlocation_list&, const mlocation_list&);
mlocation_list support(mlocation_list);

// Describe an unbranched cable in the morphology.
//
// Cables are a representation of a closed interval of a branch in a morphology.
// They may be zero-length, and fork points in the morphology may have multiple,
// equivalent zero-length cable representations.

struct mcable {
    // The id of the branch on which the cable lies.
    msize_t branch;

    // Relative location of the end points on the branch.
    // 0 ≤ prox_pos ≤ dist_pos ≤ 1
    double prox_pos; // ∈ [0,1]
    double dist_pos; // ∈ [0,1]

    friend mlocation prox_loc(const mcable& c) {
        return {c.branch, c.prox_pos};
    }
    friend mlocation dist_loc(const mcable& c) {
        return {c.branch, c.dist_pos};
    }

    // branch ≠ npos, and 0 ≤ prox_pos ≤ dist_pos ≤ 1
    friend bool test_invariants(const mcable&);
    friend std::ostream& operator<<(std::ostream&, const mcable&);
};

ARB_DEFINE_LEXICOGRAPHIC_ORDERING(mcable, (a.branch,a.prox_pos,a.dist_pos), (b.branch,b.prox_pos,b.dist_pos));

using mcable_list = std::vector<mcable>;
std::ostream& operator<<(std::ostream& o, const mcable_list& c);
// Tests whether each cable in the list satisfies the invariants for a cable,
// and that the cables in the vector are ordered.
bool test_invariants(const mcable_list&);

struct isometry {
    // Represent a 3-d isometry as a rotation (via quaterion q_)
    // and a subsequent translation (tx_, ty_, tz_).

    isometry() = default;

    // Composition: rotations are interpreted as applied to intrinsic
    // coordinates (composed by right multiplication), and translations
    // as applied to absolute coordinates (composed by addition).

    friend isometry operator*(const isometry& a, const isometry& b) {
        return isometry(b.q_*a.q_, a.tx_+b.tx_, a.ty_+b.ty_, a.tz_+b.tz_);
    }

    template <typename PointLike>
    PointLike apply(PointLike p) const {
        auto w = quaternion(p.x, p.y, p.z)^q_;
        p.x = w.x+tx_;
        p.y = w.y+ty_;
        p.z = w.z+tz_;
        return p;
    }

private:
    using quaternion = arb::math::quaternion;
    quaternion q_{1, 0, 0, 0};
    double tx_ = 0, ty_ = 0, tz_ = 0;

    isometry(quaternion q, double tx, double ty, double tz):
        q_(std::move(q)), tx_(tx), ty_(ty), tz_(tz) {}

public:
    static isometry translate(double x, double y, double z) {
        return isometry(quaternion{1, 0, 0, 0}, x, y, z);
    }

    template <typename PointLike>
    static isometry translate(const PointLike& p) {
        return translate(p.x, p.y, p.z);
    }

    // Rotate about axis in direction (ax, ay, az) by theta radians.
    static isometry rotate(double theta, double ax, double ay, double az) {
        double norm = std::sqrt(ax*ax+ay*ay+az*az);
        double vscale = std::sin(theta/2)/norm;

        return isometry(quaternion{std::cos(theta/2), ax*vscale, ay*vscale, az*vscale}, 0, 0, 0);
    }

    template <typename PointLike>
    static isometry rotate(double theta, const PointLike& p) {
        return rotate(theta, p.x, p.y, p.z);
    }
};

} // namespace arb

ARB_DEFINE_HASH(arb::mcable, a.branch, a.prox_pos, a.dist_pos);
ARB_DEFINE_HASH(arb::mlocation, a.branch, a.pos);
