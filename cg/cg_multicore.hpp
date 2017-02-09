#include <fstream>
#include <map>
#include <random>
#include <string>

#include <algorithms.hpp>
#include <memory/memory.hpp>
#include <util/span.hpp>

#include <json/json.hpp>

template <typename C>
void print_vec(const C& container, const char* name) {
    std::cout << name << "  ";
    for (auto v: container) {
        std::cout << v << " ";
    }
    std::cout << "\n";
}

namespace nest {
namespace mc {
namespace cg_cpu {

/// define the real and index types
using value_type = double;
using size_type  = unsigned;

/// define storage types
using array  = memory::host_vector<value_type>;
using iarray = memory::host_vector<size_type>;

using view       = typename array::view_type;
using const_view = typename array::const_view_type;

using iview       = typename iarray::view_type;
using const_iview = typename iarray::const_view_type;

// return scalar inner product of x and y
inline value_type dot(const_view x, const_view y) {
    value_type sum = 0;
    for (auto i: util::make_span(0u, x.size())) {
        sum += x[i]*y[i];
    }
    return sum;
}

// return scalar value of the 2 norm of x
inline value_type norm2(const_view x) {
    value_type sum = 0;
    for (auto v: x) sum += v*v;
    return std::sqrt(sum);
}

// y <- y + alpha*x
inline void axpy(view y, const_view x, value_type alpha) {
    for (auto i: util::make_span(0u, y.size())) {
        y[i] += alpha*x[i];
    }
}

struct matrix {
    array d;
    array u;
    iarray p;

    matrix() = default;

    matrix(
        const std::vector<value_type>& d,
        const std::vector<value_type>& u,
        const std::vector<size_type>&  p
    ):
        d(d.begin(), d.end()),
        u(u.begin(), u.end()),
        p(p.begin(), p.end())
    {
        check_dimensions();
    }

    void check_dimensions() {
        const auto n = size();
        EXPECTS(n==u.size());
        EXPECTS(n==p.size());
    }

    void hines_solve(view x, const_view rhs) {
        const unsigned n = size();

        memory::copy(rhs, x);

        // backward sweep
        for(unsigned i=n-1; i>0u; --i) {
            auto factor = u[i] / d[i];
            x[p[i]] -= factor * x[i];
            d[p[i]] -= factor * u[i];
        }
        x[0] = x[0]/d[0];

        // forward sweep
        for(unsigned i=1u; i<n; ++i) {
            x[i] -= u[i] * x[p[i]];
            x[i] /= d[i];
        }
    }

    bool load_from_file(std::string fname, array& x, array& rhs) {
        auto fid = std::ifstream(fname);
        if (!fid.is_open()) {
            return false;
        }
        nlohmann::json as_json;
        fid >> as_json;

        unsigned n = as_json["size"];

        auto& rng_d = as_json["d"];
        d = array(std::begin(rng_d), std::end(rng_d));

        auto& rng_u = as_json["u"];
        u = array(std::begin(rng_u), std::end(rng_u));

        auto& rng_rhs = as_json["rhs"];
        rhs = array(std::begin(rng_rhs), std::end(rng_rhs));

        auto& rng_p = as_json["p"];
        p = iarray(std::begin(rng_p), std::end(rng_p));

        auto& rng_x = as_json["voltage"];
        x = array(std::begin(rng_x), std::end(rng_x));

        if (n!=d.size() || n!=u.size() || n!=rhs.size() || n!=x.size() || n!=p.size()) {
            return false;
        }

        return true;
    }

    unsigned size() const {
        return d.size();
    }
};

struct CSR_matrix {
    iarray row_ptr;
    iarray cols;
    array values;

    CSR_matrix(const matrix& other) {
        const auto n = other.size();
        const auto nnz = 3u*n-2u;

        //
        // find the row pointers
        //

        // the number of nonzeros in each row of L+D is trivial to find (no need to look at the indexes in p)
        // note that we count the nonzeros in row i in row_ptr[i+1]
        row_ptr = iarray(n+1, 2); // all row bar the first have 2 nonzeros in each row of L+D
        row_ptr[1] = 1; // the first row has only the one nonzero on D (none in L)
        row_ptr[0] = 0;
 
        // count the number of nonzeros in each row of U and add them to row_ptr
        for (auto i=1u; i<n; ++i) {
            row_ptr[other.p[i]+1]++;
        }

        // row_ptr can now be found using a scan
        std::partial_sum(row_ptr.begin(), row_ptr.end(), row_ptr.begin());

        // find the columns and nonzero values
        cols = iarray(nnz, -1);
        values = array(nnz);


        // step 1: add contributions from L and D
        std::vector<unsigned> counts(n, 2);
        cols[0] = 0;
        counts[0] = 1;
        values[0] = other.d[0];
        for (auto i=1u; i<n; ++i) {
            const auto start = row_ptr[i];
            values[start] = other.u[i];     // contribution from L
            cols[start]   = other.p[i];     // ... in column p[i]
            values[start+1] = other.d[i];   // contribution from D
            cols[start+1]   = i;            // ... in column i
        }

        // step 2: add contributions from U
        for (auto i=1u; i<n; ++i) {
            const auto row = other.p[i];
            const auto start = row_ptr[row];
            auto& row_count = counts[row];
            values[start + row_count] = other.u[i];
            cols[start + row_count]   = i;
            ++row_count;
        }
    }

    unsigned size() const {
        return row_ptr.size()-1;
    }

    // pretty printer... don't use with large matrices
    void print() {
        const auto n = size();
        for (auto i=0u; i<n; ++i) {
            auto pos = row_ptr[i];
            const auto end  = row_ptr[i+1];
            for (auto col=0u; col<n; ++col) {
                if (pos<end && col==cols[pos]) {
                    printf("%6g  ", values[pos++]);
                }
                else {
                    printf("     .  ");
                }
            }
            printf("\n");
        }
    }
};

struct SGS_precon {
    CSR_matrix A;
    value_type omega;

    SGS_precon(const matrix& other, value_type om=1):
        A(other), omega(om)
    {}

    // perform symmetric gauss seidel preconditioner sweep.
    //      z = inv(M)*r
    // This is performed in two sweeps
    //      forward sweep: only uses values in L+D because we assume that
    //      initially z=0. we make an optimization to this effect.
    void operator() (view z, const_view r) {
        const auto n = A.size();

        /// forward sweep ///

        // The code below takes advantage there being one and one only
        // nonzero on each row of the lower triangle L. This is because
        // each node has only one parent (except for the root node which has
        // no parent, and is handled seperately for the row=0 case.)

        // the first row is a special case (no contribution from L)
        z[0] = r[0]/A.values[0];
        for (auto row=1u; row<n; ++row) {
            const auto pos = A.row_ptr[row];
            // subtract lower triangle contribution
            z[row] = r[row] - omega*A.values[pos]*z[A.cols[pos]];
            // then divide by the diagonal
            z[row] /= A.values[pos+1];
        }

        /// backward sweep ///

        for (int row=n-1; row>=0; --row) {
            value_type diag = 0;
            value_type sum = r[row];

            for (auto pos=A.row_ptr[row]; pos<A.row_ptr[row+1]; pos++) {
                const auto col = A.cols[pos];
                if (col==unsigned(row)) {
                    diag=A.values[pos]; // record the value on the diagonal
                }
                else {
                    sum -= omega*A.values[pos]*z[col];
                }
            }
            z[row] = sum / diag;
        }
    }
};


// y <- beta*y + alpha*A*x
void gemv(view y, const matrix& A, const_view x, value_type alpha, value_type beta) {
    unsigned n = A.d.size();

    const auto& u = A.u;
    const auto& d = A.d;
    const auto& p = A.p;

    EXPECTS(y.size()==n);
    EXPECTS(x.size()==n);

    y[0] = beta*y[0] + alpha*d[0]*x[0];
    for (unsigned i=1u; i<n; ++i) {
        unsigned j = p[i];
        value_type part = alpha*u[i];
        // We take advantage of the invariant j<i for all i>0.
        // Hence, y[j] has already been updated with the beta*y[j]
        // contribution.
        y[i] = beta*y[i] + alpha*d[i]*x[i];
        y[i] += part*x[j];
        y[j] += part*x[i];
    }
}

// y <- beta*y + alpha*A*x
void gemv(view y, const CSR_matrix& A, const_view x, value_type alpha, value_type beta) {
    unsigned n = A.size();

    EXPECTS(y.size()==n);
    EXPECTS(x.size()==n);

    for (auto row=0u; row<n; ++row) {
        y[row] = beta*y[row];
        for (auto pos=A.row_ptr[row]; pos<A.row_ptr[row+1]; ++pos) {
            y[row] += alpha*A.values[pos]*x[A.cols[pos]];
        }
    }
}

void scale(view x, value_type s) {
    for (auto& v: x) {
        v *= s;
    }
}

value_type normalize(view x) {
    const auto nrm = norm2(x);
    scale(x, 1/nrm);
    return nrm;
}

value_type power_method(const CSR_matrix& A, unsigned max_iters, value_type tol, bool verbose=false) {
    const auto n = A.size();
    array x_old(n);
    std::default_random_engine gen;
    std::uniform_real_distribution<value_type> dist(0, 1);
    std::generate(
        std::begin(x_old), std::end(x_old),
        [&gen, &dist](){return dist(gen);});
    normalize(x_old);

    array x_new(n);

    value_type lambda = 0;

    for (auto i=0u; i<max_iters; ++i) {
        gemv(x_new, A, x_old, 1, 0);
        const auto lambda_new = normalize(x_new);
        const auto error = std::abs(lambda-lambda_new);
        if (verbose) {
            printf ("iter %4d error %10e : lambda %f\n", i, error, lambda_new);
        }
        if (error < tol) {
            return lambda_new;
        }
        lambda = lambda_new;
        std::swap(x_old, x_new);
    }

    // return NaN on failure to converge
    return std::numeric_limits<value_type>::quiet_NaN();
}

// calculates the optimal value of omega for the SOR method
// this formulation is for tridiagonal matrices, which we use as "near enough"
value_type optimal_omega(const matrix& A, bool verbose=false) {
    // explicitly form the jacobi iteration matrix
    CSR_matrix Gj(A);
    for (auto row=0u; row<A.size(); ++row) {
        value_type s = -1/A.d[row];
        for (auto pos=Gj.row_ptr[row]; pos<Gj.row_ptr[pos+1]; ++pos) {
            Gj.values[pos] *= s;
            if (Gj.cols[pos]==row) {
                Gj.values[pos] += 1;
            }
        }
    }

    auto rho = power_method(Gj, 100, 1e-3, verbose);

    return 2. / (1 + std::sqrt(1 - rho*rho));
}

bool cg(view x, const matrix& A, const_view rhs, value_type tol, unsigned max_iters) {
    unsigned n = A.size();

    auto test_convergence = [tol] (value_type rho) {return std::sqrt(rho)<tol;};

    // calculate the residual r = rhs - A*x
    array r = rhs;
    gemv(r, A, x, -1, 1);

    value_type rho_new = dot(r, r);
    value_type rho_old = rho_new;
    if (test_convergence(rho_new)) {
        return true;
    }

    // initially p = r
    array p = r;

    // scratch working space for matrix vector product A*p
    array Ap = array(n);

    bool converged = false;

    unsigned k=0;
    for (; k<max_iters; ++k) {
        printf("CG iterations  %-4d error: %.8e\n",
            int(k), float(std::sqrt(rho_new)));

        // Ap = A*p
        gemv(Ap, A, p, 1, 0);

        // alpha = rho_old / p'*Ap;
        value_type alpha = rho_old / dot(p, Ap);

        // x = x + alpha*p
        axpy(x, p, alpha);

        // r = r - alpha*Ap
        axpy(r, Ap, -alpha);

        rho_new = dot(r, r);

        if (test_convergence(rho_new)) {
            converged = true;
            break;
        }

        // p = r + rho_new/rho_old*p
        value_type ratio = rho_new/rho_old;
        for (auto i: util::make_span(0u, n)) {
            p[i] = r[i] + ratio*p[i];
        }

        rho_old = rho_new;
    }
    printf("CG iterations  %-4d error: %.8e\n",
        int(k+1), float(std::sqrt(rho_new)));

    return converged;
}

int pcg(view x, const matrix& A, const_view rhs, value_type tol, unsigned max_iters, bool verbose=false) {
    const auto& d = A.d;

    // the preconditioner
    SGS_precon apply_precon(A, 1.9);

    unsigned n = d.size();

    auto test_convergence = [tol] (value_type rho) {return std::sqrt(rho)<tol;};

    // scratch working space
    array p(n);     // search direction
    array w(n);     // A*p
    array z(n);     // inv(M)*r

    // calculate the residual: r_0 = rhs - A*x_0
    array r = rhs;
    gemv(r, A, x, -1, 1);

    // return early if the initial estimate satisfies the convergence criteria
    value_type rho = dot(r, r);
    if (test_convergence(rho)) {
        return 0;
    }

    // z_0 = inv(M) * r_0
    apply_precon(z, r);
    value_type tau = dot(r,z);

    // set intial serach direction p_0 = z_0
    memory::copy(z, p);

    bool converged = false;

    unsigned k = 0;
    for (; k<max_iters; ++k) {
        if (verbose) {
            printf("CG iterations  %-4d error: %.8e\n",
                int(k), float(std::sqrt(rho)));
        }

        // w = A*p
        gemv(w, A, p, 1, 0);

        // alpha = tau / p'*w;
        auto alpha = tau / dot(p,w);

        // x = x + alpha*p
        axpy(x, p, alpha);

        // r = r - alpha*w
        axpy(r, w, -alpha);

        rho = dot(r, r);

        if (test_convergence(rho)) {
            converged = true;
            break;
        }

        apply_precon(z, r);

        // p = z + beta*p, where beta=tau_{j+1}/tau_{j}
        auto tau_new = dot(r,z);
        auto beta = tau_new/tau;
        for (auto i: util::make_span(0u, n)) {
            p[i] = z[i] + beta*p[i];
        }

        tau = tau_new;
    }
    if (verbose) {
        printf("CG iterations  %-4d error: %.8e\n",
            int(k+1), float(std::sqrt(rho)));
    }

    return converged ? k+1 : -1;
}

} // namespace cg_cpu
} // namespace mc
} // namespace nest
