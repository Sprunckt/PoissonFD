#include "mesh.hpp"


class Model {
public:
    Model(double xmin_, double xmax_, double ymin_, double ymax_, double dx_):
        mesh(xmin_, xmax_, ymin_, ymax_, dx_) {
        int * shape = mesh.get_shape();
        permittivity = Array2(shape[0], shape[1]);
        permittivity.fill(1.);
        bc = Array2(shape[0], shape[1]);  // boundary conditions
        bc.fill(0.);
        electric_field = Array2(shape[0], shape[1]);
        electric_field.fill(0.);
    }

    double * get_position(int i, int j) const;  // get position from cell index

    int * get_cell(double x, double y) const;  // get cell index from position

    void solve(double tol, int maxiter, double omega, bool verbose);  // solve for electric field

protected:    
    Mesh mesh;
    Array2 permittivity;  // permittivity array
    Array2 bc;  // boundary conditions : 0 for interior node, 1 for Dirichlet
    Array2 electric_field;  // electric field
};