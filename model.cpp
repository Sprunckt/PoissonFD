#include "model.hpp"


void Model::solve(double tol, int maxiter, double omega, bool verbose) {
    electric_field = electric_field.solve(bc, permittivity, tol, maxiter, omega, verbose);
}