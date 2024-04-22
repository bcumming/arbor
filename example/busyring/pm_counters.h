#pragma once

#include <memory>
#include <string>

struct instrument {
    // return true if successful
    virtual bool is_enabled() {return false;};

    // start instrumentation
    virtual void start() {};

    // stop instrumentation
    virtual void stop() {};

    // generate one line string that describes the result
    virtual void print_result(const std::string& prefix) {};

    virtual ~instrument() {};
};


std::unique_ptr<instrument> create_pm_instrument();
