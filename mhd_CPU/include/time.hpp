#pragma once
#include <cassert>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

template <typename Real>
struct Time{
    Real time, tend;
    Real dt_output;
    Real dt;
    int n_step, n_output, n_output_digits;


    Time(Real tend_, Real dt_output_, int n_output_digits_)
        :tend(tend_), dt_output(dt_output_), n_output_digits(n_output_digits_) {
            assert(tend > 0);
            assert(dt_output > 0);

            dt = 0.0;
            time = 0;
            n_step = 0;
            n_output = 0;
        }
            
    void update() {
        time += dt;
        n_step++;
    };
    

    static Time from_config(const json& config_json) {
        return Time(
            config_json["Time"]["tend"],
            config_json["Time"]["dt_output"],
            config_json["Time"]["n_output_digits"]);
    }
};