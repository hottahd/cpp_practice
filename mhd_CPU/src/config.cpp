#include "config.hpp"

#include <cassert>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

void Config::create_save_directory() const {
    fs::path dir(save_dir);
    if (!fs::exists(dir)) {
        fs::create_directories(dir);
    }
}

void Config::save() const {
    json j;

    std::ofstream ofs(save_dir + "/config.json");
    assert(ofs.is_open());
    ofs << j.dump(4);
}

