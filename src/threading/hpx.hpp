#pragma once

namespace arb {

struct hpx_guard {
    hpx_guard(int argc, char** argv);
    ~hpx_guard();
};

} // namespace arb
