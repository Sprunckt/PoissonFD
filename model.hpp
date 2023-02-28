#include "mesh.hpp"


class Model {
public:
    Model(double xmin_, double xmax_, double ymin_, double ymax_, double dx_,
         std::vector<Shape*> dielectrics, std::vector<Shape*> conductors,
         std::vector<double> permittivities, std::vector<double> potentials);

    double * get_position(int i, int j) const;  // get position from cell index

    int * get_cell(double x, double y) const;  // get cell index from position

    void solve(double tol, int maxiter, double omega, bool verbose);  // solve for electric field

    Mesh & get_mesh() {return mesh;}  // get mesh

protected:    
    Mesh mesh; // in future, we may have multiple meshes
    std::vector<double> permittivities;  // permittivity values for each dielectric
    std::vector<double> potentials;  // potential values for each conductor
    std::vector<Shape*> dielectrics;  // shapes of dielectrics
    std::vector<Shape*> conductors;  // shapes of conductors
};