
#include "linalg.hpp"
#include <vector>
#include <tuple>


class Shape;

// mesh class that handles nodes positions and indexing
class Mesh {
public:
    Mesh(double xmin_, double xmax_, double ymin_, double ymax_, double dx_):   // default: interior is not computed
        xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), dx(dx_) {
        nx =(int) (xmax - xmin)/dx + 1;
        ny = (int) (ymax - ymin)/dx + 1;
        interior = Array2(ny, nx);
    }

    Mesh(double xmin_, double xmax_, double ymin_, double ymax_, double dx_,
         std::vector<Shape> dielectrics, std::vector<Shape> conductors,
         std::vector<double> permittivities, std::vector<double> potentials) { 

        Mesh(xmin_, xmax_, ymin_, ymax_, dx_);
        //todo: compute border indices, interior of domain (first conductor is enclosure), add dielectrics and conductors

    }

    // todo : add constructor with enclosure shape, interior computing method

    double * get_position(int i, int j) const;  // get position from cell index
    void get_position(int i, int j, double * pos) const;  // get position from cell index (inplace)

    int * get_cell(double x, double y) const;  // get cell index from position
    void get_cell(double x, double y, int * cell) const;  // get cell index from position (inplace)

    int * get_shape() const;  // get shape of mesh
    inline double get_dx() const {return dx;}  // get mesh spacing



protected:
    Array2 interior;  // 1 for interior nodes, 0 for boundary 
    Array2 bc_val;  // potential values at boundary nodes

    double xmin, xmax, ymin, ymax, dx;
    int nx, ny;
};




class Shape {
public:
    // abstract method to compute border indices on a mesh
    virtual std::tuple<std::vector<int>, std::vector<int>> get_border_indices(Mesh mesh) = 0;
};

class Circle: public Shape {
public:
    Circle(double x_, double y_, double r_): x(x_), y(y_), r(r_) {}
    std::tuple<std::vector<int>, std::vector<int>> get_border_indices(Mesh mesh);

    protected:
    double x, y, r;
};

class Polygon: public Shape {
public:
    Polygon(std::vector<double> x_, std::vector<double> y_): x(x_), y(y_) {}

    protected:
    std::vector<double> x, y;
};