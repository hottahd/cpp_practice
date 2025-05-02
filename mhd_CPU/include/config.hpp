#pragma once

#include <string>
#include <cassert>
#include <filesystem>
#include <fstream>

#include <nlohmann/json.hpp>

namespace fs = std::filesystem;
using json = nlohmann::json;

struct Config {
    std::string load_filepath;
    json config_json;
    std::string save_dir;

    Config(const std::string& load_filepath_)
        : load_filepath(load_filepath_) {
            assert(!load_filepath.empty());
            if (!std::filesystem::exists(load_filepath)) {
                throw std::runtime_error("Config file not found: " + load_filepath);
            }
            std::ifstream ifs(load_filepath);
            ifs >> config_json;
            save_dir = config_json["base"]["save_dir"].get<std::string>();
        }

    

    void create_save_directory() const;
    void save() const;
};
