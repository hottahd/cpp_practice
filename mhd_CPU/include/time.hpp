#pragma once
#include <cassert>
#include <nlohmann/json.hpp>

#include "utility.hpp"

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
    
    void save(const Config& config) const {
        std::ostringstream fname;
        fname << config.save_dir << "/time." << util::zfill(n_step, n_output_digits) << ".txt";
        std::ofstream ofs(fname.str());
        assert(ofs.is_open());
        ofs << time << "\n";
        ofs << n_output << "\n";
        ofs << n_step << "\n";

        std::ofstream ofs_step(config.save_dir + "/time_output.txt");
        assert(ofs_step.is_open());
        ofs_step << n_output << "\n";
    }

    static Time from_config(const json& json_obj) {
        return Time(
            json_obj["time"]["tend"],
            json_obj["time"]["dt_output"],
            json_obj["time"]["n_output_digits"]);
    }
};