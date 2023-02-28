#include "model.hpp"


Model::Model(double xmin_, double xmax_, double ymin_, double ymax_, double dx_,
         std::vector<Shape*> dielectrics, std::vector<Shape*> conductors,
         std::vector<double> permittivities, std::vector<double> potentials): 
         mesh(xmin_, xmax_, ymin_, ymax_, dx_, dielectrics, conductors, permittivities, potentials) {
    this->permittivities = permittivities;
    this->potentials = potentials;
    this->dielectrics = dielectrics;
    this->conductors = conductors;
}



void Model::solve(double tol, int maxiter, double omega, bool verbose) {
    this->mesh.solve(tol, maxiter, omega, verbose);
}