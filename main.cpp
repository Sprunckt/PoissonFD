#include <iostream>
#include "linalg.hpp"
#include "model.hpp"


int main() {
    double dx = 0.01;
    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;

    double rinside = 0.5;
    double routside = 1.0;
    double rconduct = 0.25; 

    // create a circle
    Circle enclosure = Circle(0.0, 0.0, routside);
    Circle inner_dielectric = Circle(0.0, 0.0, rinside);
    Circle conductor = Circle(0.0, 0.0, rconduct);

    Mesh mesh(xmin, xmax, ymin, ymax, dx);

    auto borderind = enclosure.get_border_indices(mesh);
    auto interiorind = get_interior_indices(borderind.first, borderind.second);

    // create array to test the border indices
    int shape[2] = {mesh.get_shape()[0], mesh.get_shape()[1]};
    Array2 test(shape[0], shape[1]);
    test.fill(0.0);
    for (uint i=0; i < borderind.first.size(); i++) {
        test(borderind.first[i], borderind.second[i]) = 1.0;
    }
    for (uint i=0; i < interiorind.first.size(); i++) {
        test(interiorind.first[i], interiorind.second[i]) = 2.0;
    }

    // std::cout << test << std::endl;

    // create model
    Model model(xmin, xmax, ymin, ymax, dx, {&inner_dielectric}, {&enclosure, &conductor}, {5.0}, {0.0, 1.0});
    
    int* shapem = model.get_mesh().get_shape();
    std::cout << "Nb cells: " << shapem[0]*shapem[1] << std::endl;

    // solve for electric field
    model.solve(1e-6, 100000, 1.9, true);
 
    model.get_mesh().get_electric_field().to_file("output.txt");
}