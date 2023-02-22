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
    // precompute the scheme coefficients depending on permittivity
    Array2 up_coeff(ny, nx);
    Array2 down_coeff(ny, nx);
    Array2 left_coeff(ny, nx);
    Array2 right_coeff(ny, nx);

    for (int i=1; i < ny-1; i++)
    {
        for (int j=1; j < nx-1; j++)
        {
            float center_coeff = 1./(permittivity(i,j) + permittivity(i-1,j) + 
                                    permittivity(i,j-1) + permittivity(i-1,j-1));
            right_coeff(i,j) = 0.5*(permittivity(i,j) + permittivity(i-1,j))*center_coeff;
            up_coeff(i,j) = 0.5*(permittivity(i,j) + permittivity(i,j-1))*center_coeff;
            left_coeff(i,j) = 0.5*(permittivity(i,j-1) + permittivity(i-1,j-1))*center_coeff;
            down_coeff(i,j) = 0.5*(permittivity(i-1,j) + permittivity(i-1,j-1))*center_coeff;
        }
    }



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
                if (bc(i,j) == 0)  // interior node
                {
                    currval = (*this)(i,j);
                    (*this)(i,j) = ((1-omega)*currval + 
                                     omega*((*this)(i+1,j)*up_coeff(i,j) + (*this)(i-1,j)*down_coeff(i,j) + 
                                     (*this)(i,j+1)*right_coeff(i,j) + (*this)(i,j-1)*left_coeff(i,j)));
                    
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
