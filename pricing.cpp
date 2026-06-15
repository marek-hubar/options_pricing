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
    double dividend_yield; // As fraction of the stock price
    double sigma;
    double risk_free_rate;
};


class EuropeanOptionPricer {
private:
    // Standard Normal CDF
    double normal_cdf(double z) {
        return 0.5 * std::erfc(-z * std::sqrt(0.5));
    }

    void generate_s_grid_exponential(const OptionSpec &option, const MarketEnvironment &market, std::span<double> s_grid) {
        int N = s_grid.size();
        double x_center = std::log(market.underlying_price / option.strike);
        double std_devs = std::sqrt(option.expiration) * market.sigma * 3.0;

        double x_min = x_center - std_devs;
        double x_max = x_center + std_devs;

        double delta_x = (x_max - x_min) / (N - 1);

        for (int i=0; i<N; i++) {
            s_grid[i] = option.strike * std::exp(x_min + i * delta_x);
        }
    }

    void generate_payoff(const OptionSpec &option, std::span<double> s_grid, std::span<double> payoff) {
        int N = s_grid.size();
        if (option.type == call) {
            for (int i=0; i<N; i++) {
                payoff[i] = std::max(s_grid[i] - option.strike, 0.0);
            }
        } else {
            for (int i=0; i<N; i++) {
                payoff[i] = std::max(option.strike - s_grid[i], 0.0);
            }
        }
    }

public:
    double price(const OptionSpec &option_specification, const MarketEnvironment &market_env) {
        // First we calculate the put option price with the specified parameters
        if (option_specification.expiration <= 0.0) {
            if (option_specification.type == put)
                return std::max(option_specification.strike - market_env.underlying_price, 0.0);
            return std::max(market_env.underlying_price - option_specification.strike, 0.0);
        }
        double exponential = std::exp(-market_env.risk_free_rate * option_specification.expiration);
        double normalizer = market_env.sigma * std::sqrt(option_specification.expiration);
        double log_S_over_K = std::log(market_env.underlying_price / option_specification.strike);
        double d1 = (log_S_over_K + (market_env.risk_free_rate - market_env.dividend_yield + 0.5 * market_env.sigma * market_env.sigma) * option_specification.expiration) / normalizer;
        double d2 = d1 - normalizer;
        double option_price = option_specification.strike * normal_cdf(-d2) * exponential - market_env.underlying_price * std::exp(-market_env.dividend_yield * option_specification.expiration) * normal_cdf(-d1);

        // Then we use put-call parity to calculate the call price (if the option is a call)
        if (option_specification.type == call) {
            option_price += market_env.underlying_price * std::exp(-market_env.dividend_yield * option_specification.expiration) - option_specification.strike * std::exp(-market_env.risk_free_rate * option_specification.expiration);
        }
        return option_price;
    }
    //
    // time_steps is the number of temporal points, not the intervals between them: |----|----| <- this is time_step=3
    void calculate_price_surface(const OptionSpec &option_specification, const MarketEnvironment &market_env, int time_steps, int spatial_resolution, std::span<double> price_surface_buffer, std::span<double> s_grid_buffer, std::span<double> t_grid_buffer) {
        generate_s_grid_exponential(option_specification, market_env, s_grid_buffer);
        // Populate the first row of the price_surface_buffer with the payoff
        generate_payoff(option_specification, s_grid_buffer, price_surface_buffer);

        // First calculate the price surface for a put option
        std::vector<double> log_s_over_K_vector(spatial_resolution, 0.0);
        for (int j=0; j<spatial_resolution; j++)
            log_s_over_K_vector[j] = std::log(s_grid_buffer[j] / option_specification.strike);

        t_grid_buffer[0] = 0.0;
        int time_intervals = time_steps > 1 ? time_steps - 1 : 1;
        double dt = option_specification.expiration / double(time_intervals);
        for (int i=1; i<time_steps; i++) {
            double t = i * dt;
            t_grid_buffer[i] = t;
            double exponential = std::exp(-market_env.risk_free_rate * t);
            double normalizer = market_env.sigma * std::sqrt(t);

            for (int j=0; j<spatial_resolution; j++) {
                double d1 = (log_s_over_K_vector[j] + (market_env.risk_free_rate - market_env.dividend_yield + 0.5 * market_env.sigma * market_env.sigma) * t) / normalizer;
                double d2 = d1 - normalizer;
                price_surface_buffer[i * spatial_resolution + j] = option_specification.strike * normal_cdf(-d2) * exponential - s_grid_buffer[j] * std::exp(-market_env.dividend_yield * t) * normal_cdf(-d1);
            }
        }
    
        // If a call option is specified, use the put-call parity to calculate the call option surface
        if (option_specification.type == call) {
            for (int i=1; i<time_steps; i++) {
                double time_component = option_specification.strike * std::exp(-market_env.risk_free_rate * t_grid_buffer[i]);
                for (int j=0; j<spatial_resolution; j++) {
                    price_surface_buffer[i * spatial_resolution + j] += s_grid_buffer[j] * std::exp(-market_env.dividend_yield * t_grid_buffer[i]) - time_component;
                }
            }
        }
    }
};

struct TridiagonalMatrix {
    std::span<const double> a; // Lower Diagonal (1 to N-1)
    std::span<const double> b; // Diagonal       (0 to N-1)
    std::span<const double> c; // Upper Diagonal (0 to N-2)
};

class ThomasSolver {
private:
    std::vector<double> c_prime; // Cached upper-diagonal multipliers
    std::vector<double> d_prime; // Reusable workspace for the RHS
    std::vector<double> denom;   // Cached denominators (b_i - a_i * c'_{i-1})
    size_t N;

public:
    // Constructor: Pre-computes everything that doesn't change over time
    ThomasSolver(std::span<const double> a, 
                 std::span<const double> b, 
                 std::span<const double> c) 
    {
        N = b.size();
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
    struct CrankNicolsonStepper {
        std::optional<ThomasSolver> solver;
        std::vector<double> a, b, c, d;
        int N;
        
        // Asymptotic boundary weights
        double c_L1, c_L2, c_R1, c_R2;

        CrankNicolsonStepper(std::span<const double> x_grid, double d_tau, double alpha) {
            N = static_cast<int>(x_grid.size());
            
            double dx = x_grid[1] - x_grid[0]; 
            double E = std::exp(dx);

            // Left Boundary Multipliers (V_SS = 0 as S -> 0)
            c_L1 = (1.0 + 1.0 / E) * std::exp(alpha * dx);
            c_L2 = std::exp((2.0 * alpha - 1.0) * dx);

            // Right Boundary Multipliers (V_SS = 0 as S -> infinity)
            c_R1 = (1.0 + E) * std::exp(-alpha * dx);
            c_R2 = std::exp((1.0 - 2.0 * alpha) * dx);

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

            b[0] += a[0] * c_L1;
            c[0] -= a[0] * c_L2;
            a[0] = 0.0;

            a[N-3] -= c[N-3] * c_R2;
            b[N-3] += c[N-3] * c_R1;
            c[N-3] = 0.0;

            solver.emplace(a, b, c);
        }

        void advance_one_step(std::span<const double> price_curr, std::span<double> price_next, int current_step_count, double d_tau, double alpha, double beta, std::span<const double> base_intrinsic) {
            for (int i=0; i<N-2; i++) {
                d[i] = -a[i] * price_curr[i] + (2.0 - b[i]) * price_curr[i+1] - c[i] * price_curr[i+2];
            }

            solver->solve(a, d, price_next.subspan(1, N-2));

            price_next[0] = c_L1 * price_next[1] - c_L2 * price_next[2];
            price_next[N-1] = c_R1 * price_next[N-2] - c_R2 * price_next[N-3];

            double current_tau = (current_step_count + 1) * d_tau;
            double time_scalar = std::exp(-beta * current_tau);

            for (int j=0; j<N; j++) {
                price_next[j] = std::max(price_next[j], time_scalar * base_intrinsic[j]);
            }
        }
    };

    static void generate_x_grid_uniform(const OptionSpec &option, const MarketEnvironment &market, std::span<double> x_grid) {
        double x_center = std::log(market.underlying_price / option.strike);
        double std_devs = std::sqrt(option.expiration) * market.sigma * 3.0;

        double x_min = x_center - std_devs;
        double x_max = x_center + std_devs;

        int N = static_cast<int>(x_grid.size());
        double delta_x = (x_max - x_min) / (N - 1);

        for (int i=0; i<N; i++) {
            x_grid[i] = x_min + i * delta_x;
        }
    }

    static void generate_payoff(const OptionSpec &option, double alpha, std::span<const double> x_grid, std::span<double> payoff) {
        int N = static_cast<int>(x_grid.size());
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
    void calculate_price_surface(const OptionSpec &option_specification, const MarketEnvironment &market_env, int time_steps, std::span<double> price_surface_buffer, std::span<double> s_grid_buffer, std::span<double> t_grid_buffer) {
        int spatial_resolution = static_cast<int>(s_grid_buffer.size());

        std::vector<double> x_grid(spatial_resolution);
        std::vector<double> base_intrinsic(spatial_resolution);

        double tau_max = 0.5 * market_env.sigma * market_env.sigma * option_specification.expiration;
        int time_intervals = time_steps > 1 ? time_steps - 1 : 1;
        double d_tau = tau_max / double(time_intervals);

        generate_x_grid_uniform(option_specification, market_env, x_grid);

        double c_L = (x_grid[1] - x_grid[0]) / (x_grid[2] - x_grid[1]);
        double c_R = (x_grid[spatial_resolution-1] - x_grid[spatial_resolution-2]) / (x_grid[spatial_resolution-2] - x_grid[spatial_resolution-3]);

        double k1 = 2 * (market_env.risk_free_rate - market_env.dividend_yield) / (market_env.sigma * market_env.sigma);
        double k2 = 2.0 * market_env.risk_free_rate / (market_env.sigma * market_env.sigma);

        double alpha = -0.5 * (k1 - 1);
        double beta = -alpha * alpha - k2;

        generate_payoff(option_specification, alpha, x_grid, base_intrinsic);

        auto first_row = price_surface_buffer.subspan(0, spatial_resolution);
        std::copy(base_intrinsic.begin(), base_intrinsic.end(), first_row.begin());

        CrankNicolsonStepper stepper(x_grid, d_tau, alpha);

        for (int i=0; i<time_steps-1; i++) {
            auto row_tau_i = price_surface_buffer.subspan(spatial_resolution * i, spatial_resolution);
            auto row_tau_i_plus_one = price_surface_buffer.subspan(spatial_resolution * (i+1), spatial_resolution);

            stepper.advance_one_step(row_tau_i, row_tau_i_plus_one, i, d_tau, alpha, beta, base_intrinsic);
        }

        for (int j = 0; j < spatial_resolution; j++) {
            s_grid_buffer[j] = option_specification.strike * std::exp(x_grid[j]);
        }

        for (int i = 0; i < time_steps; i++) {
            double current_tau = i * d_tau;

            t_grid_buffer[i] = (2.0 * current_tau) / (market_env.sigma * market_env.sigma);

            double time_growth = std::exp(beta * current_tau);

            for (int j = 0; j < spatial_resolution; j++) {
                double space_growth = std::exp(alpha * x_grid[j]);
                double abstract_u = price_surface_buffer[i * spatial_resolution + j];

                price_surface_buffer[i * spatial_resolution + j] = option_specification.strike * abstract_u * space_growth * time_growth;
            }
        }
    }

    double price(const OptionSpec &option_specification, const MarketEnvironment &market_env, int time_steps, int spatial_resolution) {

        std::vector<double> x_grid(spatial_resolution);
        std::vector<double> base_intrinsic(spatial_resolution);

        double tau_max = 0.5 * market_env.sigma * market_env.sigma * option_specification.expiration;
        int time_intervals = time_steps > 1 ? time_steps - 1 : 1;

        generate_x_grid_uniform(option_specification, market_env, x_grid);

        double d_tau = tau_max / double(time_intervals);
        double min_dx = x_grid[1] - x_grid[0];
        for (int i=1; i<spatial_resolution-1; i++)
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
        double c_R = (x_grid[spatial_resolution-1] - x_grid[spatial_resolution-2]) / (x_grid[spatial_resolution-2] - x_grid[spatial_resolution-3]);

        double k1 = 2.0 * (market_env.risk_free_rate - market_env.dividend_yield) / (market_env.sigma * market_env.sigma);
        double k2 = 2.0 * market_env.risk_free_rate / (market_env.sigma * market_env.sigma);

        double alpha = -0.5 * (k1 - 1.0);
        double beta = -alpha * alpha - k2;

        generate_payoff(option_specification, alpha, x_grid, base_intrinsic);

        std::vector<double> u_prev(spatial_resolution, 0.0);
        std::vector<double> u_curr = base_intrinsic;

        CrankNicolsonStepper stepper(x_grid, d_tau, alpha);

        for (int i = 0; i < time_intervals; i++) {
            std::swap(u_prev, u_curr);
            auto u_curr_inside = std::span<double>(u_curr).subspan(1, spatial_resolution - 2);
            stepper.advance_one_step(u_prev, u_curr, i, d_tau, alpha, beta, base_intrinsic);
        }

        double target_x = std::log(market_env.underlying_price / option_specification.strike);

        size_t idx = 0;
        for (size_t j = 0; j < static_cast<size_t>(spatial_resolution - 1); j++) {
            if (target_x >= x_grid[j] && target_x <= x_grid[j+1]) {
                idx = j;
                break;
            }
        }

        double x0 = x_grid[idx];
        double x1 = x_grid[idx+1];
        double u0 = u_curr[idx];
        double u1 = u_curr[idx+1];

        double weight = (target_x - x0) / (x1 - x0);
        double interpolated_u = u0 + weight * (u1 - u0);

        double time_growth = std::exp(beta * tau_max);
        double space_growth = std::exp(alpha * target_x);

        return option_specification.strike * interpolated_u * space_growth * time_growth;
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
    OptionSpec call_option = {100.0, 1.0, call};
    OptionSpec put_option = {100.0, 1.0, put};
    MarketEnvironment market = {100.0, 0.0, 0.2, 0.05};

    int spatial_resolution =  437;
    int temporal_resolution = 4750;

    
    std::vector<double> price_surface_buffer(spatial_resolution * temporal_resolution, 0.0);
    std::vector<double> price_surface_buffer1(spatial_resolution * temporal_resolution, 0.0);
    std::vector<double> price_surface_buffer2(spatial_resolution * temporal_resolution, 0.0);
    std::vector<double> s_grid_buffer(spatial_resolution, 0.0);
    std::vector<double> t_grid_buffer(temporal_resolution, 0.0);

    AmericanOptionPricer american_pricer;
    EuropeanOptionPricer european_pricer;

    std::cout << "American Call Option Price: " << american_pricer.price(call_option, market, temporal_resolution, spatial_resolution) << std::endl;
    std::cout << "American Put Option Price:  " << american_pricer.price(put_option, market, temporal_resolution, spatial_resolution) << std::endl;
    std::cout << "European Call Option Price: " << european_pricer.price(call_option, market) << std::endl;
    std::cout << "European Put Option Price:  " << european_pricer.price(put_option, market) << std::endl;

    european_pricer.calculate_price_surface(call_option, market, temporal_resolution, spatial_resolution, price_surface_buffer1, s_grid_buffer, t_grid_buffer);
    american_pricer.calculate_price_surface(call_option, market, temporal_resolution, price_surface_buffer2, s_grid_buffer, t_grid_buffer);

    for (int i=0; i<spatial_resolution * temporal_resolution; i++)
        price_surface_buffer[i] = price_surface_buffer2[i] - price_surface_buffer1[i];

    export_surface_to_json(s_grid_buffer, t_grid_buffer, price_surface_buffer2, "test/price_surface.json");

    return 0;
}
#endif
