#include <iostream>

#include <hpx/hpx_start.hpp>
#include <hpx/hpx_suspend.hpp>
#include <hpx/include/apply.hpp>
#include <hpx/util/yield_while.hpp>

int main(int argc, char** argv) {
    hpx::start(nullptr, 0, argv);

    ///// workaround to wait for hpx runtime to start
    hpx::runtime* rt = hpx::get_runtime_ptr();
    hpx::util::yield_while([rt]() { return rt->get_state() < hpx::state_running; });
    /////////////////////////////////////////////////

    hpx::suspend();

    hpx::apply([](){return hpx::finalize();});
    std::cout << "hoi zäme" << std::endl; 
    hpx::stop();
}
