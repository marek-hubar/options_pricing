#include "pricing.cpp"

extern "C" {

void* pricer_create() {
    static int dummy;
    return &dummy;
}

void pricer_destroy(void* handle) {
    (void)handle;
}

double pricer_price(void* handle,
                    int pricing_model,
                    int option_type,
                    double strike,
                    double expiration,
                    double underlying_price,
                    double dividend_yield,
                    double sigma,
                    double risk_free_rate,
                    int time_steps,
                    int spatial_resolution) {
    (void)handle;
    OptionSpec spec = {strike, expiration, option_type == 0 ? call : put};
    MarketEnvironment market = {underlying_price, dividend_yield, sigma, risk_free_rate};
    if (pricing_model == 0) {
        EuropeanOptionPricer pricer;
        return pricer.price(spec, market);
    } else {
        AmericanOptionPricer pricer;
        return pricer.price(spec, market, time_steps, spatial_resolution);
    }
}

void pricer_surface(void* handle,
                    int pricing_model,
                    int option_type,
                    double strike,
                    double expiration,
                    double underlying_price,
                    double dividend_yield,
                    double sigma,
                    double risk_free_rate,
                    int time_steps,
                    int spatial_resolution,
                    double* surface_buf,
                    double* s_grid_buf,
                    double* t_grid_buf) {
    (void)handle;
    OptionSpec spec = {strike, expiration, option_type == 0 ? call : put};
    MarketEnvironment market = {underlying_price, dividend_yield, sigma, risk_free_rate};

    int N = spatial_resolution;
    std::span<double> surface_span(surface_buf, static_cast<size_t>(time_steps) * N);
    std::span<double> s_grid_span(s_grid_buf, N);
    std::span<double> t_grid_span(t_grid_buf, time_steps);

    if (pricing_model == 0) {
        EuropeanOptionPricer pricer;
        pricer.calculate_price_surface(spec, market, time_steps, N,
                                       surface_span, s_grid_span, t_grid_span);
    } else {
        AmericanOptionPricer pricer;
        pricer.calculate_price_surface(spec, market, time_steps,
                                       surface_span, s_grid_span, t_grid_span);
    }
}

}
