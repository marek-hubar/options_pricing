#include "pricing.cpp"

extern "C" {

void* pricer_create(int spatial_resolution) {
    return new (std::nothrow) AmericanOptionPricer(spatial_resolution);
}

void pricer_destroy(void* handle) {
    delete static_cast<AmericanOptionPricer*>(handle);
}

int pricer_spatial_res(void* handle) {
    return static_cast<AmericanOptionPricer*>(handle)->spatial_res();
}

double pricer_price(void* handle,
                    int option_type,
                    double strike,
                    double expiration,
                    double underlying_price,
                    double sigma,
                    double risk_free_rate,
                    int time_steps) {
    auto* p = static_cast<AmericanOptionPricer*>(handle);
    OptionSpec spec = {strike, expiration, option_type == 0 ? call : put};
    MarketEnvironment market = {underlying_price, sigma, risk_free_rate};
    return p->price(spec, market, time_steps);
}

void pricer_surface(void* handle,
                    int option_type,
                    double strike,
                    double expiration,
                    double underlying_price,
                    double sigma,
                    double risk_free_rate,
                    int time_steps,
                    double* surface_buf,
                    double* s_grid_buf,
                    double* t_grid_buf) {
    auto* p = static_cast<AmericanOptionPricer*>(handle);
    OptionSpec spec = {strike, expiration, option_type == 0 ? call : put};
    MarketEnvironment market = {underlying_price, sigma, risk_free_rate};

    int N = p->spatial_res();
    std::span<double> surface_span(surface_buf, static_cast<size_t>(time_steps) * N);
    std::span<double> s_grid_span(s_grid_buf, N);
    std::span<double> t_grid_span(t_grid_buf, time_steps);

    p->calculate_price_surface(spec, market, time_steps,
                               surface_span, s_grid_span, t_grid_span);
}

}
