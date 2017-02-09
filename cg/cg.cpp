#include <random>

#include "cg_multicore.hpp"

#include "../tests/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


using namespace nest::mc;
namespace cpu = nest::mc::cg_cpu;

using T = cpu::value_type;
using I = cpu::size_type;
using array = cpu::array;
using matrix = cpu::matrix;
using CSR_matrix = cpu::CSR_matrix;
using cpu::norm2;

TEST(cpu, load) {
    matrix m;
    array x;
    array rhs;
    auto success = m.load_from_file("matrix_1000.json", x, rhs);
    EXPECT_TRUE(success) << "--- unable to open file for matrix input ---";
    EXPECT_EQ(m.size(), 301u);
    EXPECT_EQ(m.size(), x.size());
    EXPECT_EQ(m.size(), rhs.size());
}

TEST(cpu, tridiag_gemv) {
    // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
    // choose the rhs such that the solution to the linear system is "ones"
    std::vector<I> p = {0,0,1,2,3};
    const auto n = p.size();
    std::vector<T> d(n, 4);
    std::vector<T> u(n, -1);
    std::vector<T> rhs(n, 2);
    rhs.back() = rhs.front() = 3;
    array x(n, 1);
    matrix A(d, u, p);

    // perform matrix multiply
    array y(n);
    array result(n);
    gemv(result, A, x, 1, 0);

    for (auto i: util::make_span(0u, n))
        EXPECT_EQ(result[i], rhs[i]);
}

TEST(cpu, branch_gemv) {
    // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
    // choose the rhs such that the solution to the linear system is "ones"
    const auto n = 10u;
    std::vector<I> p(n, 0);
    std::vector<T> d(n, 4);
    std::vector<T> u(n, 0.1);
    array x(n, 1);
    std::vector<T> rhs(n, 4.1);
    rhs.front() = 4. + 0.1*(n-1);
    matrix A(d, u, p);

    // perform matrix multiply
    array y(n);
    array result(n);
    gemv(result, A, x, 1, 0);

    for (auto i: util::make_span(0u, n))
        EXPECT_NEAR(result[i], rhs[i], 1e-14);
}

const auto tol = 1e-8;
const auto max_iters = 100u;
const auto tridiag_n = 1000;

TEST(cpu, cg_tridiag) {
    // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
    // choose the rhs such that the solution to the linear system is "ones"
    const auto n = tridiag_n;
    std::vector<I> p(n);
    p[0] = 0;
    std::iota(p.begin()+1, p.end(), 0);
    std::vector<T> d(n, 4);
    std::vector<T> u(n, -1);

    // set the rhs vector
    array rhs(n, 2);
    rhs[0] = 3;
    rhs[n-1] = 3;

    array x(n, 0);
    matrix A(d, u, p);

    // perform matrix multiply
    auto success = cg(x, A, rhs, tol, max_iters);
    ASSERT_TRUE(success);

    for (auto i: util::make_span(0u, n))
        EXPECT_NEAR(x[i], 1.0, 10*tol);

    //print_vec(x, "cg solution");
    //A.hines_solve(x, rhs);
    //print_vec(x, "hines solution");
}

TEST(cpu, pcg_tridiag) {
    // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
    // choose the rhs such that the solution to the linear system is "ones"
    const auto n = tridiag_n;
    std::vector<I> p(n);
    p[0] = 0;
    std::iota(p.begin()+1, p.end(), 0);
    std::vector<T> d(n, 4);
    std::vector<T> u(n, -1);

    // set the rhs vector
    array rhs(n, 2);
    rhs[0] = 3;
    rhs[n-1] = 3;

    array x(n, 0);
    matrix A(d, u, p);

    // perform matrix multiply
    auto success = pcg(x, A, rhs, tol, max_iters);
    ASSERT_TRUE(success);

    for (auto i: util::make_span(0u, n))
        EXPECT_NEAR(x[i], 1.0, 10*tol);

    //print_vec(x, "cg solution");
    //A.hines_solve(x, rhs);
    //print_vec(x, "hines solution");
}

TEST(cpu, cg_applied) {
    // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
    // choose the rhs such that the solution to the linear system is "ones"
    matrix A;
    array x_cg;
    array rhs;
    auto load_status = A.load_from_file("matrix_1000.json", x_cg, rhs);
    ASSERT_TRUE(load_status) << "  -- unable to open file matrix_1000.json --";
    const auto n = A.size();

    //CSR_matrix C(A);
    //C.print();

    // this accelerates convergence very nicely...
    for (auto& d: A.d) d += 0.01;

    // perform matrix multiply
    auto solve_status = cg(x_cg, A, rhs, tol, max_iters);
    ASSERT_TRUE(solve_status);

    array x_hines(n, 0);
    A.hines_solve(x_hines, rhs);

    // compute the relative error
    for (auto i: util::make_span(0u,n))
        x_cg[i] -= x_hines[i];
    auto error = norm2(x_cg)/norm2(x_hines);
    std::cout << "difference between direct and cg : " << error << "\n";
}

TEST(cpu, power_method) {
    auto verbose = false;
    {
        // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
        std::vector<I> p = {0,0,1};
        const auto n = p.size();
        std::vector<T> d(n, 2);
        std::vector<T> u(n, -1);
        matrix A(d, u, p);

        auto radius = power_method(A, 10, 1e-4, verbose);

        EXPECT_NEAR(radius, 3.41421, 1e-4);
    }
    {
        // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
        std::vector<I> p = {0,0,1,2};
        const auto n = p.size();
        std::vector<T> d(n, 4);
        std::vector<T> u(n, -1);
        matrix A(d, u, p);

        auto radius = power_method(A, 100, 1e-5, verbose);

        EXPECT_NEAR(radius, 5.61803, 1e-4);
    }
}

TEST(cpu, sor_omega) {
    auto verbose = true;
    {
        // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
        std::vector<I> p = {0,0,1};
        const auto n = p.size();
        std::vector<T> d(n, 2);
        std::vector<T> u(n, -1);
        matrix A(d, u, p);

        auto omega = optimal_omega(A, verbose);
        std::cout << " omega : " << omega << "\n";
    }
}

TEST(cpu, pcg_applied) {

    // tridiagonal matrix with 4 on diagonal and -1 on off-diagonals
    // choose the rhs such that the solution to the linear system is "ones"
    matrix A;
    array x_cg;
    array rhs;
    for (auto i=0; i<4700; i+=100) {
        auto fname = "matrix_" + std::to_string(i) + ".json";
        auto load_status = A.load_from_file(fname, x_cg, rhs);
        ASSERT_TRUE(load_status) << "  -- unable to open file matrix_1000.json --";

        if (!i) {
            auto omega = optimal_omega(A, true);
            std::cout << " -- optimal omega " << omega << " --\n";
        }

        // this accelerates convergence very nicely...
        //for (auto& d: A.d) d += 0.1;

        // perform matrix multiply
        auto iters = pcg(x_cg, A, rhs, tol, max_iters);
        ASSERT_TRUE(iters!=-1);
        std::cout << "matrix " << i << " required " << iters << " iterations\n";
    }

    /*
    const auto n = A.size();
    array x_hines(n, 0);
    A.hines_solve(x_hines, rhs);

    // compute the relative error
    for (auto i: util::make_span(0u,n))
        x_cg[i] -= x_hines[i];
    auto error = norm2(x_cg)/norm2(x_hines);
    std::cout << "difference between direct and cg : " << error << "\n";
    */
}

TEST(cpu, CSR_init) {
    {
        std::cout << "*** tridiagonal matrix ***\n";
        const auto n = 7;
        std::vector<I> p(n);
        p[0] = 0;
        std::iota(p.begin()+1, p.end(), 0);
        std::vector<T> d(n, 4);
        std::vector<T> u(n, -1);

        matrix H(d, u, p);

        CSR_matrix A(H);
        A.print();
        std::cout << "\n";
    }
    {
        std::cout << "*** all roads lead to soma matrix ***\n";
        const auto n = 7;
        std::vector<I> p(n, 0);
        std::vector<T> d(n, 4);
        std::vector<T> u(n, 0.1);

        matrix H(d, u, p);

        CSR_matrix A(H);
        A.print();
        std::cout << "\n";
    }

}
