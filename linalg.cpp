#include "linalg.hpp"

void Array2::fill(double value) {
    for (int i = 0; i < nx*ny; i++) {
        data[i] = value;
    } 
}


void Array2::to_file(std::string filename) {
    std::ofstream file;
    file.open(filename);
    file << *this;
    file.close();
}
