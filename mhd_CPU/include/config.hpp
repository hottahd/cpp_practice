#pragma once

#include <string>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <type_traits>

#include <nlohmann/json.hpp>

#include "types.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

struct Config {
    std::string load_filepath;
    json json_obj;
    std::string save_dir;

    Config(const std::string& load_filepath_)
        : load_filepath(load_filepath_) {
            assert(!load_filepath.empty());
            if (!std::filesystem::exists(load_filepath)) {
                throw std::runtime_error("Config file not found: " + load_filepath);
            }
            std::ifstream ifs(load_filepath);
            ifs >> json_obj;

            if constexpr (std::is_same_v<Real, float>) {
                json_obj["type"]["Real"] = "float";
            } else if constexpr (std::is_same_v<Real, double>) {
                json_obj["type"]["Real"] = "double";
            } else {
                throw std::runtime_error("Unsupported Real type");
            }
            save_dir = json_obj.at("base").at("save_dir").get<std::string>();
        }


    void create_save_directory() const;
    void save() const;
};
