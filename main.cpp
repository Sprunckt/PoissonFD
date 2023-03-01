#include <iostream>
#include "linalg.hpp"
#include "model.hpp"


double exact_sol(double x, double y, double rin, double rout, double vin, double vout,
                 double rdiel, double permin, double permout) {
    /* Exact solution of the Poisson equation for two concentric conductors of radiuses rint and rout at potentials 
    vin and vout, with two dielectrics of radiuses rdiel and rout at permittivities permin and permout.*/

    double r = sqrt(x*x + y*y);
    double c1 = (vout - vin)/((permin/permout)*log(rout/rdiel) + log(rdiel/rin));
    double c2 = vin - log(rin)*c1;
    double c3 = (permin/permout)*c1;
    double c4 = vout - log(rout)*c3;
    if (r <= rin) {
        return vin;
    }
    else if (r < rdiel) {
        return c1 * log(r) + c2;
    }
    else if (r < rout) {
        return c3 * log(r) + c4;
    }
    else {
        return vout;
    }
}


int main(int argc, char *argv[]) {
    // read input parameters, first argument is test case id
    int test_case = 0;
    if (argc > 1) {
        test_case = std::stoi(argv[1]);
    }

    double dx = 0.003;
    

    if (test_case == 0)
    {
        double xmin = -1.;
        double xmax = 1.;
        double ymin = -1.;
        double ymax = 1.;
        double rinside = 0.5;
        double routside = 1.0;
        double rconduct = 0.25; 
        double perm_diel = 5.0;
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
        Model model(xmin, xmax, ymin, ymax, dx, {&inner_dielectric}, {&enclosure, &conductor}, {perm_diel}, {0.0, 1.0});
        
        int* shapem = model.get_mesh().get_shape();
        std::cout << "Nb cells: " << shapem[0]*shapem[1] << std::endl;

        // solve for electric field
        model.solve(1e-6, 100000, 1.9, true);
        model.get_mesh().get_electric_field().to_file("computed_sol.txt");

        // compute exact solution at mesh points
        Array2 exact_sol_array(shape[0], shape[1]);
        exact_sol_array.fill(0.0);
        double * pos = new double[2];
        for (int i=0; i < shape[0]; i++) {
            for (int j=0; j < shape[1]; j++) {
                model.get_mesh().get_position(i, j, pos);
                exact_sol_array(i, j) = exact_sol(pos[0], pos[1], rconduct, routside, 1., 0., rinside, perm_diel, 1.0);
            }
        }
        exact_sol_array.to_file("exact_sol.txt");
    }
   

    
    
}