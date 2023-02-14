#include "model.hpp"

int * Model::get_cell(double x, double y) const {
    int * cell = new int[2];
    cell[0] = (int) (y - ymin)/dx;
    cell[1] = (int) (x - xmin)/dx;
    return cell;
}

double * Model::get_position(int i, int j) const {
    double * position = new double[2];
    position[0] = xmin + j*dx;
    position[1] = ymin + i*dx;
    return position;
}