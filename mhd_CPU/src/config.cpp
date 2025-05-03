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
    std::ofstream ofs(save_dir + "/config.json");
    assert(ofs.is_open());
    ofs << json_obj.dump(4);
}

