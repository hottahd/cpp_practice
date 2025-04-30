#include <fstream>
#include <nlohmann/json.hpp>
#include "grid.hpp"
#include "config.hpp"

using json = nlohmann::json;

template <typename Real>
Grid<Real>::Grid(int i_size_, int margin_, Real xmin_, Real xmax_)
    : i_size(i_size_), margin(margin_), 
      xmin(xmin_), xmax(xmax_) {
    i_total = i_size + 2*margin;
    dx = (xmax - xmin) / i_size;

    assert(i_size > 0);
    assert(margin >= 0);
    assert(xmax > xmin);

    x.resize(i_total);
    x[margin] = xmin + 0.5*dx;
    for (int i = margin + 1; i < i_total; ++i) {
        x[i] = x[i - 1] + dx;
    }
    for (int i = margin - 1; i >= 0; --i) {
        x[i] = x[i + 1] - dx;
    }
}

template <typename Real>
void Grid<Real>::save(const Config& config) const {
    json j;
    j["i_size"] = i_size;
    j["margin"] = margin;
    j["xmin"] = xmin;
    j["xmax"] = xmax;
    j["dx"] = dx;

    std::ofstream ofs_json(config.save_dir + "/grid.json");
    assert(ofs_json.is_open());
    ofs_json << j.dump(4);

    std::ofstream ofs_bin(config.save_dir + "/grid.bin", std::ios::binary);
    assert(ofs_bin.is_open());
    ofs_bin.write(reinterpret_cast<const char*>(x.data()), sizeof(Real)*x.size());
}

template struct Grid<double>;
template struct Grid<float>;