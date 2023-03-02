#ifndef LINALG_H
#define LINALG_H

#include <cmath>
#include <iostream>
#include <fstream>


class Array2 {
public:
    Array2() {data = nullptr;}

    Array2(int n1, int n2): nx(n2), ny(n1) {
        data = new double[nx*ny];
    }

    Array2(Array2 const& a) {
        int * s = a.shape();
        nx = s[0];
        ny = s[1];
        data = new double[nx*ny];
        for (int i = 0; i < nx*ny; i++) {
            data[i] = a.data[i];
        }
    }

    Array2& operator=(Array2 const& a) {
        *this = Array2(a);
        return *this;
    }

    ~Array2() {
        delete[] data;
    }

    inline double& operator()(int i, int j) {return data[i*nx + j];}

    inline double operator()(int i, int j) const {return data[i*nx + j];}
    
    inline int * shape() const {int * s = new int[2]; s[1] = nx; s[0] = ny; return s;}
    inline int size() const {return nx*ny;}

    void fill(double value);

    friend std::ostream& operator<<(std::ostream& os, const Array2& a) {
        for (int i = 0; i < a.ny; i++) {
            for (int j = 0; j < a.nx; j++) {
                os << a(i, j) << " ";
            }
            os << std::endl;
        }
                    os << std::endl;

        return os;
    }
    
    void to_file(std::string filename);


protected:
    int nx, ny;
    double* data;
};

#endif
