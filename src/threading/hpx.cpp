#include "hpx.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <hpx/hpx_start.hpp>
#include <hpx/hpx_suspend.hpp>
#include <hpx/include/apply.hpp>
#include <hpx/util/yield_while.hpp>
#include <hpx/runtime/get_os_thread_count.hpp>

namespace arb {

hpx_guard::hpx_guard(int argc, char** argv) {
    std::vector<std::string> cfg = {
        "hpx.os_threads=4",
    };

    std::cout << ":::: hpx: starting" << std::endl;
    // -- hack -- a dummy init function myust be passed to hpx::start, to avoid memory corruption
    auto dummy = [](int, char**) { return 42; };
    hpx::start(dummy, 0, argv, cfg);

    // -- hack -- wait for hpx to start
    auto rt = hpx::get_runtime_ptr();
    hpx::util::yield_while([rt]() { return rt->get_state()<hpx::state_running; });

    auto num_threads = hpx::get_os_thread_count();
    std::cout << ":::: hpx: started with " << num_threads << " os threads" << std::endl;
}

hpx_guard::~hpx_guard() {
    // The hpx runtime has to be running before calling finalize.
    // The good news: calling hpx::resume is ok if the runtime is not suspended.
    hpx::resume();
    hpx::apply([](){hpx::finalize();});
    hpx::stop();
    std::cout << ":::: hpx: stopped\n";
}

} // namespace arb
