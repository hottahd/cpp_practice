#pragma once
#include <string>
#include <filesystem>
#include <fstream>
#include <cassert>


namespace fs = std::filesystem;

struct Config {
    std::string save_dir;

    Config(const std::string& save_dir_)
        : save_dir(save_dir_) {
            assert(!save_dir.empty());
        }

    void create_save_directory() const;
    void save() const;
};
