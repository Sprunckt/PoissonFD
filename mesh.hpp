#ifndef MESH_H
#define MESH_H

#include "linalg.hpp"
#include <vector>
#include <tuple>
#include <algorithm>
#include <queue>


class Shape;
std::pair<std::vector<int>, std::vector<int>> get_interior_indices(std::vector<int>, std::vector<int>);


// mesh class that handles nodes positions and indexing
class Mesh {
public:
    Mesh(double xmin_, double xmax_, double ymin_, double ymax_, double dx_);

    Mesh(double xmin_, double xmax_, double ymin_, double ymax_, double dx_,
         std::vector<Shape*> dielectrics, std::vector<Shape*> conductors,
         std::vector<double> permittivities, std::vector<double> potentials);

    double * get_position(int i, int j) const;  // get position from cell index
    void get_position(int i, int j, double * pos) const;  // get position from cell index (inplace)

    int * get_cell(double x, double y) const;  // get cell index from position
    void get_cell(double x, double y, int * cell) const;  // get cell index from position (inplace)

    int * get_shape() const;  // get shape of mesh
    inline double get_dx() const {return dx;}  // get mesh spacing

    // solve for electric field
    Array2 solve(double tol, int maxiter, double omega, bool verbose);

    Array2 & get_electric_field() {return electric_field;}  // get electric field

protected:
    double xmin, xmax, ymin, ymax, dx;
    int nx, ny;

    Array2 permittivity;  // permittivity values at nodes
    Array2 electric_field;  // potential values at boundary nodes
    Array2 node_type;  // 1 for interior nodes, 0 for boundary nodes, -1 for exterior nodes
};




class Shape {
public:
    // abstract method to compute border indices on a mesh
    virtual std::pair<std::vector<int>, std::vector<int>> get_border_indices(Mesh mesh) = 0;
};

class Circle: public Shape {
public:
    Circle(double x_, double y_, double r_): x(x_), y(y_), r(r_) {}
    std::pair<std::vector<int>, std::vector<int>> get_border_indices(Mesh mesh);

    protected:
    double x, y, r;
};

class Polygon: public Shape {  // todo: implement border indices for polygon
public:
    Polygon(std::vector<double> x_, std::vector<double> y_): x(x_), y(y_) {}

    std::pair<std::vector<int>, std::vector<int>> get_border_indices(Mesh mesh);

    protected:
    std::vector<double> x, y;
};

#endif