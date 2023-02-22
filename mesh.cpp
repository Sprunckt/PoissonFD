#include "mesh.hpp"


void Mesh::get_cell(double x, double y, int * cell) const {
    double yval = (y - ymin)/dx;
    double xval = (x - xmin)/dx;
    cell[0] = (int) static_cast<int>(yval < 0 ? yval - 0.5 : yval + 0.5);
    cell[1] = (int) static_cast<int>(xval < 0 ? xval - 0.5 : xval + 0.5);
}


int * Mesh::get_cell(double x, double y) const {
    int * cell = new int[2];
    get_cell(x, y, cell);
    return cell;
}


void Mesh::get_position(int i, int j, double * position) const {
    position[0] = xmin + j*dx;
    position[1] = ymin + i*dx;
}


double * Mesh::get_position(int i, int j) const {
    double * position = new double[2];
    get_position(i, j, position);
    return position;
}


int * Mesh::get_shape() const {
    int * shape = new int[2];
    shape[0] = ny;
    shape[1] = nx;
    return shape;
}


std::tuple<std::vector<int>, std::vector<int>> Circle::get_border_indices(Mesh mesh) {
    // compute the border indices of a circle on the given mesh using the midpoint algorithm
    std::vector<int> bindx;
    std::vector<int> bindy;

    // project circle center on mesh
    int * cell = mesh.get_cell(x, y);
    // project circle radius on mesh
    int * cellr = mesh.get_cell(x + r, y + r);

    int currx = cellr[1];
    int curry = 0;  // start at the rightmost point of the circle

    // compute the first octant
    int d = 1 - cellr[1];

    while (currx > curry) {  // compute the 8 octants in clockwise order
        curry++;
        if (d > 0) { 
            currx--;
            d += 2*(curry - currx) + 1;
        }  
        else {
            d += 2*curry + 1;
        }

        // add the 8 octants to the border indices
        bindx.push_back(cell[1] + currx);
        bindy.push_back(cell[0] + curry);

        bindx.push_back(cell[1] - currx);
        bindy.push_back(cell[0] + curry);

        bindx.push_back(cell[1] + currx);
        bindy.push_back(cell[0] - curry);

        bindx.push_back(cell[1] - currx);
        bindy.push_back(cell[0] - curry);

        if (currx != curry) {
            bindx.push_back(cell[1] + curry);
            bindy.push_back(cell[0] + currx);

            bindx.push_back(cell[1] - curry);
            bindy.push_back(cell[0] + currx);

            bindx.push_back(cell[1] + curry);
            bindy.push_back(cell[0] - currx);

            bindx.push_back(cell[1] - curry);
            bindy.push_back(cell[0] - currx);
        }
    }

    return make_tuple(bindx, bindy);
}
