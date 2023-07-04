#include <string>
#include <vector>

#include <arbor/profile/meter.hpp>

#include "hardware/power.hpp"

namespace arb {
namespace profile {

class power_meter: public meter {
    std::vector<hw::energy_size_type> readings_;
    //std::string name_;
    std::string fname_;

public:
    power_meter(std::string name) {
        fname_=name;
        /*
        fname_ = name;
        if (name == std::string("energy")) {
            name_ = name;
        }
        else if (auto it = name.find("_energy"); it != std::string::npos) {
            name_ = name;
            name_.erase(name_.begin()+it, name_.end());
        }
        */
    }

    std::string name() override {
        return fname_;
    }

    std::string units() override {
        return "J";
    }

    std::vector<double> measurements() override {
        std::vector<double> diffs;

        for (auto i=1ul; i<readings_.size(); ++i) {
            diffs.push_back(readings_[i]-readings_[i-1]);
        }

        return diffs;
    }

    void take_reading() override {
        readings_.push_back(hw::energy(fname_.c_str()));
    }
};

meter_ptr make_power_meter(std::string fname) {
    if (!arb::hw::has_energy_measurement()) {
        return nullptr;
    }
    return meter_ptr(new power_meter(fname));
}

} // namespace profile
} // namespace arb
