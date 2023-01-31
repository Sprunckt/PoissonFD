#include <iostream>
#include "linalg.hpp"

int main() {
    int ny =500;
    int nx = 750;
    std::cout << "Nb cells: " << ny*nx << std::endl;

    // create a 2D array of size nyXnx
    Array2 a(ny, nx);   
    a.fill(0.0);
    
    // create the associated boundary condition array (0 for interior, 1 for border)
    Array2 bc(ny, nx);
    bc.fill(0.0);

    // fill every border with 1.0
    for (int i = 0; i < ny; i++) { // i is the y coordinate
        a(i, 0) = 1.0;
        a(i, nx-1) = 1.0;
        bc(i, 0) = 1.0;
        bc(i, nx-1) = 1.0;
    }

    for (int j = 0; j < nx; j++) {
        a(0, j) = 1.0;
        a(ny-1, j) = 1.0;
        bc(0, j) = 1.0;
        bc(ny-1, j) = 1.0;
    }

    a.solve(bc, bc, 1e-4, 100000, 1.9, true);
    a.to_file("output.txt");

    
}