#pragma once

#include <arbor/profile/meter.hpp>

namespace arb {
namespace profile {

meter_ptr make_power_meter(std::string fname);

} // namespace profile
} // namespace arb
