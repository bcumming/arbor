#pragma once

#include <cstdint>

#include <common_types.hpp>
#include <util/debug.hpp>

namespace arb {

// Information about a current time integration epoch.
// Each epoch has an integral id, that is incremented for successive epochs.
// An epoch describes the half open time interval:
//
//      time ∈ [t, min(tupper, t+dt))
//
// At the end of an epoch the solution is at min(tupper, t+dt), however events that are
// due for delivery at tupper_i are not delivered until epoch_i+1.

struct epoch {
    using index_type = std::int64_t;
    index_type id;          // integer id of the current epoch
    time_type t;            // start time of the current epoch
    time_type tint;         // length of each time interval
    time_type tupper;       // highest possible value of t
    epoch() = default;

    epoch(index_type id, time_type t, time_type tint, time_type tupper):
        id(id), t(t), tint(tint), tupper(tupper)
    {}

    epoch(const epoch& other):
        id(other.id), t(other.t), tint(other.tint), tupper(other.tupper)
    {}

    epoch& operator=(const epoch& other) {
        id = other.id;
        t = other.t;
        tint = other.tint;
        tupper = other.tupper;
        return *this;
    }

    time_type t0() const {
        return t;
    }

    time_type t1() const {
        return std::min(t+tint, tupper);
    }

    void advance() {
        t = t1();
        ++id;
    }

    epoch& operator++() {
        advance();
        return *this;
    }

    epoch operator++(int) {
        epoch e = *this;
        advance();
        return e;
    }

    bool finished() const {
        return t>=tupper;
    }

    friend epoch operator+(const epoch& ep, index_type n) {
        EXPECTS(n>=0);
        epoch e = ep;
        while (n--) ++e;
        return e;
    }

    friend std::ostream& operator<<(std::ostream& o, const epoch& e) {
        return o << "epoch(" << e.id << ", " << e.t0() << " to " << e.t1() << ")";
    }
};


} // namespace arb
