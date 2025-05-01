#pragma once
#include <vector>

template <typename T>
class Array3D {
    private:
        int i_total, j_total, k_total;
        std::vector<T> array;
    public:
        Array3D(int i_total_, int j_total_, int k_total_):
            i_total(i_total_),
            j_total(j_total_),
            k_total(k_total_),
            array(i_total_ * j_total_ * k_total_) {}

        T& operator()(int i, int j, int k) {
            return array[i * j_total * k_total + j * k_total + k];
        }
        const T& operator()(int i, int j, int k) const {
            return array[i * j_total * k_total + j * k_total + k];
        }

        T* data() { return array.data(); }
        const T* data() const { return array.data(); }

        int size_x () const { return i_total; }
        int size_y () const { return j_total; }
        int size_z () const { return k_total; }
        int size () const { return i_total * j_total * k_total; }
};