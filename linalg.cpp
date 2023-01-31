#include "linalg.hpp"

void Array2::fill(double value) {
    for (int i = 0; i < nx*ny; i++) {
        data[i] = value;
    }
}

Array2 Array2::solve(Array2 bc, Array2 permittivity, double tol, int maxiter, double omega, bool verbose){
    // apply gauss-seidel method to solve in place for the electric field (permittivity is ignored for now)
    double error = tol + 1.;
    int iter = 0;
    double currval;
    while (error > tol && iter < maxiter)
    {
        iter++;
        if (verbose && (iter % 500 == 0))
            std::cout << "Iteration " << iter << " error: " << error << std::endl;
            
        error = 0.;
        for (int i=1; i < ny-1; i++)
        {
            for (int j=1; j < nx-1; j++)
            {
                if (bc(i,j) == 0)
                {
                    currval = (*this)(i,j);
                    (*this)(i,j) = (1-omega)*currval + omega*0.25*((*this)(i+1,j) + (*this)(i-1,j) + (*this)(i,j+1) + (*this)(i,j-1));
                    
                    error += std::abs(currval - (*this)(i,j));
                }
            }
        }
    }
    return (*this);
}

void Array2::to_file(std::string filename) {
    std::ofstream file;
    file.open(filename);
    file << *this;
    file.close();
}
