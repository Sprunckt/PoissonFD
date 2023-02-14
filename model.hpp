#include "linalg.hpp"

class Model {
public:
    Model(double xmin_, double xmax_, double ymin_, double ymax_, double dx_):
        xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), dx(dx_) {
        nx =(int) (xmax - xmin)/dx + 1;
        ny = (int) (ymax - ymin)/dx + 1;
        permittivity = Array2(ny, nx);
        permittivity.fill(1.);
        bc = Array2(ny, nx);  // boundary conditions
        bc.fill(0.);
        electric_field = Array2(ny, nx);
        electric_field.fill(0.);
    }

    Model(double xmin_, double xmax_, double ymin_, double ymax_, Array2 permittivity_):
        xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_) {
        int * shape = permittivity_.shape();
        nx = shape[1];
        ny = shape[0];
        dx = (xmax - xmin)/(nx - 1);
        permittivity = permittivity_;
        bc = Array2(ny, nx);  // boundary conditions
        bc.fill(0.);
        electric_field = Array2(ny, nx);
        electric_field.fill(0.);
    }

    double * get_position(int i, int j) const;  // get position from cell index

    int * get_cell(double x, double y) const;  // get cell index from position

protected:
    double xmin, xmax, ymin, ymax, dx;
    int nx, ny;
    Array2 permittivity;
    Array2 bc;
    Array2 electric_field;
};