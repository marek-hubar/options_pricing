#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <span>
#include <fstream>
#include <string>

enum OptionType {call, put};

struct OptionSpec {
    double strike;
    double expiration;
    OptionType type;
};

struct MarketEnvironment {
    double underlying_price;
    double sigma;
    double risk_free_rate;
};

class EuropeanOptionPricer {

};

// Standard Normal PDF
double normal_pdf(double z) {
    static const double inv_sqrt_2pi = 0.3989422804014326779;
    return inv_sqrt_2pi * std::exp(-0.5 * z * z);
}
// Standard Normal CDF
double normal_cdf(double z) {
    return 0.5 * std::erfc(-z * std::sqrt(0.5));
}

double price_european_put_option(double strike_price, double expiration, double current_price, double risk_free_rate, double sigma) {
    double exponential = std::exp(-risk_free_rate * expiration);
    double normalizer = sigma * std::sqrt(expiration);
    double log_S_over_K = std::log(current_price / strike_price);
    double d1 = (log_S_over_K + (risk_free_rate - 0.5 * sigma * sigma) * expiration) / normalizer;
    double d2 = d1 + normalizer;
    return strike_price * normal_cdf(-d1) * exponential - current_price * normal_cdf(-d2);
}

double price_european_call_option(double strike_price, double expiration, double current_price, double risk_free_rate, double sigma) {
    double put_price = price_european_put_option(strike_price, expiration, current_price, risk_free_rate, sigma);
    return put_price + current_price - strike_price * std::exp(-risk_free_rate * expiration);
}

class ThomasSolver {
private:
    std::vector<double> c_prime; // Cached upper-diagonal multipliers
    std::vector<double> d_prime; // Reusable workspace for the RHS
    std::vector<double> denom;   // Cached denominators (b_i - a_i * c'_{i-1})

public:
    // Constructor: Pre-computes everything that doesn't change over time
    ThomasSolver(std::span<const double> a, 
                 std::span<const double> b, 
                 std::span<const double> c) 
    {
        const size_t N = b.size();
        c_prime.resize(N, 0.0);
        d_prime.resize(N, 0.0);
        denom.resize(N, 0.0);

        // Forward sweep for the static matrix bands
        denom[0] = b[0];
        c_prime[0] = c[0] / denom[0];


        for (size_t i = 1; i < N - 1; ++i) {
            denom[i] = b[i] - a[i] * c_prime[i - 1];
            c_prime[i] = c[i] / denom[i];
        }
        
        // The last element doesn't have a c_prime, but needs a denominator
        denom[N - 1] = b[N - 1] - a[N - 1] * c_prime[N - 2];
    }

    // Solve function: Called at every time step
    // d is the known vector from time i, x is the output vector for time i+1
    void solve(std::span<const double> a, 
               std::span<const double> d, 
               std::span<double> x) 
    {
        const size_t N = d.size();

        // 1. Forward sweep (only for the changing right-hand side 'd')
        d_prime[0] = d[0] / denom[0];

        for (size_t i = 1; i < N; ++i) {
            d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom[i];
        }

        // 2. Back substitution (solves directly into the output span 'x')
        x[N - 1] = d_prime[N - 1];

        for (int i = N - 2; i >= 0; --i) {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
        }
    }
};

class AmericanOptionPricer {
private:
    std::vector<double> x_grid;
    std::vector<double> base_intrinsic;
    int N;


    struct CrankNicolsonStepper {
        std::optional<ThomasSolver> solver;
        std::vector<double> a, b, c, d;
        int N;
        CrankNicolsonStepper(std::span<double> x_grid, double d_tau) {
            N = x_grid.size();
            double c_L = (x_grid[1] - x_grid[0]) / (x_grid[2] - x_grid[1]);
            double c_R = (x_grid[N-1] - x_grid[N-2]) / (x_grid[N-2] - x_grid[N-3]);

            a.resize(N-2, 0.0);
            b.resize(N-2, 0.0);
            c.resize(N-2, 0.0);
            d.resize(N-2, 0.0);

            for (int i=1; i<N-1; i++) {
                double h_L = x_grid[i] - x_grid[i-1];
                double h_R = x_grid[i+1] - x_grid[i];
                double h_L_plus_h_R = x_grid[i+1] - x_grid[i-1];
                a[i-1] = -d_tau / (h_L * h_L_plus_h_R);
                b[i-1] = 1.0 + d_tau / (h_L * h_R);
                c[i-1] = -d_tau / (h_R * h_L_plus_h_R);
            }
            b[0] += a[0] * (1 + c_L);
            c[0] -= a[0] * c_L;
            a[0] = 0.0;
            a[N-3] -= c[N-3] * c_R;
            b[N-3] += c[N-3] * (1 + c_R);
            c[N-3] = 0.0;

            solver.emplace(a, b, c);
        }

        // The fill price vector (of length N) should be passed into the price_curr arugment
        // and a buffer (of length N-2) for the solved interior of the next price vector should be passed to price_next
        void advance_one_step(std::span<double> price_curr, std::span<double> price_next) {
            for (int i=0; i<N-2; i++) {
                d[i] = -a[i] * price_curr[i] +  (2 - b[i]) * price_curr[i+1] - c[i] * price_curr[i+2];
            }

            solver->solve(a, d, price_next);
        }

    };

    void generate_x_grid_uniform(const OptionSpec &option, const MarketEnvironment &market, std::span<double> x_grid) {
        double x_center = std::log(market.underlying_price / option.strike);
        double std_devs = std::sqrt(option.expiration) * market.sigma * 3.0;

        double x_min = x_center - std_devs;
        double x_max = x_center + std_devs;

        double delta_x = (x_max - x_min) / (N - 1);

        for (int i=0; i<N; i++) {
            x_grid[i] = x_min + i * delta_x;
        }
    }

    void generate_payoff(const OptionSpec &option, double alpha, std::span<double> payoff) {
        if (option.type == call) {
            for (int i=0; i<N; i++) {
                payoff[i] = std::exp(-alpha * x_grid[i]) * std::max(std::exp(x_grid[i]) - 1.0, 0.0);
            }
        } else {
            for (int i=0; i<N; i++) {
                payoff[i] = std::exp(-alpha * x_grid[i]) * std::max(1.0 - std::exp(x_grid[i]), 0.0);
            }
        }
    }

public:
    AmericanOptionPricer(int spatial_resolution) {
        N = spatial_resolution;
        x_grid.resize(N);
        base_intrinsic.resize(N);
    }

    int spatial_res() const { return static_cast<int>(N); }

    // time_steps is the number of temporal points, not the intervals between them: |----|----| <- this is time_step=3
    void calculate_price_surface(OptionSpec &option_specification, MarketEnvironment &market_env, int time_steps, std::span<double> price_surface_buffer, std::span<double> s_grid_buffer, std::span<double> t_grid_buffer) {
        double tau_max = 0.5 * market_env.sigma * market_env.sigma * option_specification.expiration;
        double d_tau = tau_max / double(time_steps);

        generate_x_grid_uniform(option_specification, market_env, x_grid);

        double c_L = (x_grid[1] - x_grid[0]) / (x_grid[2] - x_grid[1]);
        double c_R = (x_grid[N-1] - x_grid[N-2]) / (x_grid[N-2] - x_grid[N-3]);

        double k = 2 * market_env.risk_free_rate / (market_env.sigma * market_env.sigma);
        double alpha = -0.5 * (k - 1);
        double beta = -0.25 * (k + 1) * (k + 1);

        generate_payoff(option_specification, alpha, base_intrinsic);

        auto first_row = price_surface_buffer.subspan(0, N);
        std::copy(base_intrinsic.begin(), base_intrinsic.end(), first_row.begin());

        CrankNicolsonStepper stepper(x_grid, d_tau);

        for (int i=0; i<time_steps-1; i++) {
            auto row_tau_i = price_surface_buffer.subspan(N * i, N);
            auto row_tau_i_plus_one = price_surface_buffer.subspan(N * (i+1), N);
            auto row_tau_i_plus_one_interior = price_surface_buffer.subspan(N * (i+1) + 1, N-2);

            stepper.advance_one_step(row_tau_i, row_tau_i_plus_one_interior);

            row_tau_i_plus_one[0] = (1 + c_L) * row_tau_i_plus_one[1] - c_L * row_tau_i_plus_one[2];
            row_tau_i_plus_one[N-1] = (1 + c_R) * row_tau_i_plus_one[N-2] - c_R * row_tau_i_plus_one[N-3];

            double current_tau = (i + 1) * d_tau;
            double time_scalar = std::exp(-beta * current_tau);
            
            for (int j=0; j<N; j++) {
                row_tau_i_plus_one[j] = std::max(row_tau_i_plus_one[j], time_scalar * base_intrinsic[j]);
            }
        }

        // 1. Populate the output Stock Price grid
        for (size_t j = 0; j < N; j++) {
            s_grid_buffer[j] = option_specification.strike * std::exp(x_grid[j]);
        }

        // 2. Populate the output Time to Maturity grid AND convert the surface
        for (int i = 0; i < time_steps; i++) {
            double current_tau = i * d_tau;
            
            // Map tau back to real physical time to maturity
            t_grid_buffer[i] = (2.0 * current_tau) / (market_env.sigma * market_env.sigma);

            // FIXED: Multiply by positive beta to invert the transformation
            double time_growth = std::exp(beta * current_tau);

            for (size_t j = 0; j < N; j++) {
                // FIXED: Multiply by positive alpha to invert the transformation
                double space_growth = std::exp(alpha * x_grid[j]);
                double abstract_u = price_surface_buffer[i * N + j];
                
                // V(S, t) = K * u(x, tau) * exp(+alpha*x + beta*tau)
                price_surface_buffer[i * N + j] = option_specification.strike * abstract_u * space_growth * time_growth;
            }
        }

    }

    double price(const OptionSpec &option_specification, const MarketEnvironment &market_env, int time_steps) {
        // 1. Setup Time and Mesh parameters
        double tau_max = 0.5 * market_env.sigma * market_env.sigma * option_specification.expiration;
        int time_intervals = time_steps > 1 ? time_steps - 1 : 1;
        // Generate grid and payoff
        generate_x_grid_uniform(option_specification, market_env, x_grid);

        double d_tau = tau_max / double(time_intervals);
        double min_dx = x_grid[1] - x_grid[0];
        for (int i=1; i<N-1; i++)
            if (x_grid[i + 1] - x_grid[i] < min_dx)
                min_dx = x_grid[i + 1] - x_grid[i];

        double gamma = d_tau / (min_dx * min_dx);

        if (gamma > 1.0) {
            std::cerr << "\n==============================================================\n"
                      << " WARNING: NUMERICAL INSTABILITY DETECTED\n"
                      << " Grid Ratio (Gamma) = " << gamma << "\n"
                      << " Gamma exceeds 1.0. The Crank-Nicolson method may produce\n"
                      << " spurious oscillations, triggering the American ratchet effect.\n"
                      << " \n"
                      << " RECOMMENDATION: To maintain stability, increase time_steps\n"
                      << " exponentially relative to spatial_resolution. If you double\n"
                      << " the spatial nodes, you must quadruple the time steps.\n"
                      << "==============================================================\n\n";
        }


        double c_L = (x_grid[1] - x_grid[0]) / (x_grid[2] - x_grid[1]);
        double c_R = (x_grid[N-1] - x_grid[N-2]) / (x_grid[N-2] - x_grid[N-3]);

        double k = 2.0 * market_env.risk_free_rate / (market_env.sigma * market_env.sigma);
        double alpha = -0.5 * (k - 1.0);
        double beta = -0.25 * (k + 1.0) * (k + 1.0);

        generate_payoff(option_specification, alpha, base_intrinsic);

        // 2. Memory Optimization: We only need two rows
        std::vector<double> u_prev(N, 0.0);
        std::vector<double> u_curr = base_intrinsic; // Starts at tau = 0

        CrankNicolsonStepper stepper(x_grid, d_tau);

        // 3. Time Loop (Swapping rows instead of building a matrix)
        for (int i = 0; i < time_intervals; i++) {
            
            // The current row becomes the previous row for this step
            std::swap(u_prev, u_curr); 

            auto u_curr_inside = std::span<double>(u_curr).subspan(1, N - 2);

            // Advance the PDE
            stepper.advance_one_step(u_prev, u_curr_inside);

            // Apply boundaries
            u_curr[0] = (1.0 + c_L) * u_curr[1] - c_L * u_curr[2];
            u_curr[N-1] = (1.0 + c_R) * u_curr[N-2] - c_R * u_curr[N-3];

            // Apply American early-exercise constraint
            double current_tau = (i + 1) * d_tau;
            double time_scalar = std::exp(-beta * current_tau);
            
            for (size_t j = 0; j < N; j++) {
                u_curr[j] = std::max(u_curr[j], time_scalar * base_intrinsic[j]);
            }
        }

        // 4. Extract the exact Option Price at the Current Stock Price
        // In log-space, the current spot price S_0 maps exactly to x_target
        double target_x = std::log(market_env.underlying_price / option_specification.strike);
        
        // Find where target_x falls on our grid
        size_t idx = 0;
        for (size_t j = 0; j < N - 1; j++) {
            if (target_x >= x_grid[j] && target_x <= x_grid[j+1]) {
                idx = j;
                break;
            }
        }

        // Linearly interpolate the abstract heat value between the two closest nodes
        double x0 = x_grid[idx];
        double x1 = x_grid[idx+1];
        double u0 = u_curr[idx];
        double u1 = u_curr[idx+1];
        
        double weight = (target_x - x0) / (x1 - x0);
        double interpolated_u = u0 + weight * (u1 - u0);

        // 5. Apply the Inverse Transformation (Only ONCE!)
        double time_growth = std::exp(beta * tau_max);
        double space_growth = std::exp(alpha * target_x);
        
        double final_price = option_specification.strike * interpolated_u * space_growth * time_growth;

        return final_price;
    }
};


// Takes real stock prices, real times, and real option prices
#ifndef __EMSCRIPTEN__
void export_surface_to_json(std::span<const double> s_grid, 
                            std::span<const double> t_grid, 
                            std::span<const double> price_surface, 
                            const std::string& filename) 
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    size_t N_x = s_grid.size();
    size_t N_y = t_grid.size();

    file << "{\n";

    // 1. Export X-Axis: Stock Prices
    file << "  \"x_stock_price\": [";
    for (size_t j = 0; j < N_x; ++j) {
        file << s_grid[j];
        if (j < N_x - 1) file << ", ";
    }
    file << "],\n";

    // 2. Export Y-Axis: Time to Maturity
    file << "  \"y_time_to_maturity\": [";
    for (size_t i = 0; i < N_y; ++i) {
        file << t_grid[i];
        if (i < N_y - 1) file << ", ";
    }
    file << "],\n";

    // 3. Export Z-Axis: Option Prices
    file << "  \"z_option_price\": [\n";
    for (size_t i = 0; i < N_y; ++i) {
        file << "    [";
        for (size_t j = 0; j < N_x; ++j) {
            file << price_surface[i * N_x + j];
            if (j < N_x - 1) file << ", ";
        }
        file << "]";
        if (i < N_y - 1) file << ",";
        file << "\n";
    }
    file << "  ]\n";
    
    file << "}\n";
    std::cout << "\nSuccessfully exported fully-transformed data to " << filename << std::endl;
}
#endif

#ifndef __EMSCRIPTEN__
int main() {
    float S = 100, K = 100, r = 0.05, sigma = 0.55, T = 1.0;
    std::cout << "Call Option Price: " << price_european_call_option(K, T, S, r, sigma) << std::endl;
    std::cout << "Put Option Price:  " << price_european_put_option(K, T, S, r, sigma) << std::endl;

    OptionSpec option = {100.0, 1, put};
    MarketEnvironment market = {100.0, 0.2, 0.05};

    int spatial_resolution =  500;
    int temporal_resolution = 10000;

    
    std::vector<double> price_surface_buffer(spatial_resolution * temporal_resolution, 0.0);
    std::vector<double> s_grid_buffer(spatial_resolution, 0.0);
    std::vector<double> t_grid_buffer(temporal_resolution, 0.0);

    AmericanOptionPricer pricer(spatial_resolution);

    pricer.calculate_price_surface(option, market, temporal_resolution, price_surface_buffer, s_grid_buffer, t_grid_buffer);
    std::cout << "Option Price: " << pricer.price(option, market, temporal_resolution) << std::endl;

    export_surface_to_json(s_grid_buffer, t_grid_buffer, price_surface_buffer, "test/price_surface.json");

    return 0;
}
#endif
