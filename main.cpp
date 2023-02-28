#include <iostream>
#include "linalg.hpp"
#include "model.hpp"


int main() {
    double dx = 0.1;
    double xmin = -1.0;
    double xmax = 1.0;
    double ymin = -1.0;
    double ymax = 1.;

    double rinside = 0.5;
    double routside = 1.0;
    double rconduct = 0.25; 

    // create a circle
    Circle enclosure = Circle(0.0, 0.0, routside);
    Circle inner_dielectric = Circle(0.0, 0.0, rinside);
    Circle conductor = Circle(0.0, 0.0, rconduct);

    // create model
    Model model(xmin, xmax, ymin, ymax, dx, {&inner_dielectric}, {&enclosure, &conductor}, {5.0}, {0.0, 1.0});
    
    int* shape = model.get_mesh().get_shape();
    std::cout << "Nb cells: " << shape[0]*shape[1] << std::endl;


    // solve for electric field
    // model.solve(1e-6, 100000, 1.9, true);
 
    // model.get_mesh().get_electric_field().to_file("output.txt");
}