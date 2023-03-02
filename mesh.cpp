#include "mesh.hpp"


Mesh::Mesh(double xmin_, double xmax_, double ymin_, double ymax_, double dx_):   // default: interior is not computed
    xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), dx(dx_), nx((int) (xmax - xmin)/dx + 2), ny((int) (ymax - ymin)/dx + 2),
    permittivity(ny,nx), electric_field(ny,nx), node_type(ny,nx) {
    node_type.fill(-1.);
    permittivity.fill(1.);
    electric_field.fill(0.);
}


Mesh::Mesh(double xmin_, double xmax_, double ymin_, double ymax_, double dx_,
           std::vector<Shape*> dielectrics, std::vector<Shape*> conductors,
           std::vector<double> permittivities, std::vector<double> potentials): Mesh(xmin_, xmax_, ymin_, ymax_, dx_){ 
    
    // check if potential and conductors have the same size
    if ((potentials.size() != conductors.size()) || conductors.size() == 0) {
        throw std::invalid_argument("potentials and conductors must have the same size and be non-empty");
    }
    if (permittivities.size() != dielectrics.size()) {
        throw std::invalid_argument("permittivities and dielectrics must have the same size");
    }

    // compute interior indices
    for (uint i=0; i < conductors.size(); i++) {
        std::pair<std::vector<int>, std::vector<int>> border_indices = conductors[i]->get_border_indices(*this);
        std::pair<std::vector<int>, std::vector<int>> interior_indices = get_interior_indices(border_indices.first, border_indices.second);
        for (uint j=0; j < interior_indices.first.size(); j++) {
            node_type(interior_indices.first[j], interior_indices.second[j]) = 1;
        }

        // set boundary values
        for (uint j=0; j < border_indices.first.size(); j++) {
            node_type(border_indices.first[j], border_indices.second[j]) = 0;
            electric_field(border_indices.first[j], border_indices.second[j]) = potentials[i];
        }
    }

    // set permittivities inside dielectrics
    for (uint i=0; i < dielectrics.size(); i++) {
        std::pair<std::vector<int>, std::vector<int>> border_indices = dielectrics[i]->get_border_indices(*this);
        std::pair<std::vector<int>, std::vector<int>> interior_indices = get_interior_indices(border_indices.first, border_indices.second);
        for (uint j=0; j < interior_indices.first.size(); j++) {
            permittivity(interior_indices.first[j], interior_indices.second[j]) = permittivities[i];
        }
    }
}


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


std::pair<std::vector<int>, std::vector<int>> Circle::get_border_indices(const Mesh mesh) {
    // compute the border indices of a circle on the given mesh using the midpoint algorithm
    std::vector<int> bindx;
    std::vector<int> bindy;
    double * xbounds = mesh.get_xbounds();
    double * ybounds = mesh.get_ybounds();
    // project circle center on mesh
    int * cell = mesh.get_cell(x, y);
    // project circle radius on mesh
    int * cellr = mesh.get_cell(xbounds[0] + r, ybounds[0] + r);

    int currx = cellr[1];
    int curry = 0;  // start at the rightmost point of the circle

    // compute the first octant
    int d = 1 - cellr[1];
    
    // add left, right top and bottom points to the border indices
    bindx.push_back(cell[1] + cellr[1]);
    bindy.push_back(cell[0]);

    bindx.push_back(cell[1]);
    bindy.push_back(cell[0] + cellr[0]);

    bindx.push_back(cell[1] - cellr[1]);
    bindy.push_back(cell[0]);

    bindx.push_back(cell[1]);
    bindy.push_back(cell[0] - cellr[0]);

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

    return make_pair(bindy, bindx);
}

std::pair<std::vector<int>, std::vector<int>> draw_line(double x1, double y1, double x2, double y2, Mesh mesh) {
    /*Draw a line using the Bresenham algorithm*/
    std::vector<int> xind;
    std::vector<int> yind;
    int * cell1 = mesh.get_cell(x1, y1);
    int * cell2 = mesh.get_cell(x2, y2);
    int xspan = cell2[1] - cell1[1];
    int yspan = cell2[0] - cell1[0];
    int xstep = xspan > 0 ? 1 : -1;
    int ystep = yspan > 0 ? 1 : -1;
    xspan = abs(xspan);
    yspan = abs(yspan);
    int x = cell1[1];
    int y = cell1[0];
    int d = 0;
    int err = xspan - yspan;
    while (true) {
        xind.push_back(x);
        yind.push_back(y);
        if (x == cell2[1] && y == cell2[0]) {
            break;
        }
        int e2 = 2*err;
        if (e2 > -yspan) {
            err -= yspan;
            x += xstep;
        }
        if (e2 < xspan) {
            err += xspan;
            y += ystep;
        }
    }
    return make_pair(yind, xind);
}



std::pair<std::vector<int>, std::vector<int>> Polygon::get_border_indices(const Mesh mesh) {
    std::vector<int> bindx;
    std::vector<int> bindy;
    for (uint i = 0; i < x.size() - 1; i++) {  // draw lines between the vertices
        auto line = draw_line(x[i], y[i], x[i+1], y[i+1], mesh);
        bindx.insert(bindx.end(), line.second.begin(), line.second.end());
        bindy.insert(bindy.end(), line.first.begin(), line.first.end());
    }
    auto line = draw_line(x[x.size() - 1], y[y.size() - 1], x[0], y[0], mesh);
    bindx.insert(bindx.end(), line.second.begin(), line.second.end());
    bindy.insert(bindy.end(), line.first.begin(), line.first.end());

    return make_pair(bindy, bindx);
}


std::pair<std::vector<int>, std::vector<int>> get_interior_indices(std::vector<int> yind, std::vector<int> xind){
    /* compute the interior indices of a polygon on the given mesh using the flood fill algorithm, 
    assume the shape has only one connected component*/

    const auto xbounds = std::minmax_element(xind.begin(), xind.end());
    const auto ybounds = std::minmax_element(yind.begin(), yind.end());
    int xmin = *xbounds.first;
    int xmax = *xbounds.second;
    int ymin = *ybounds.first;
    int ymax = *ybounds.second;


    // encasing rectangle around the border, enlarged by 2 in each direction
    Array2 labels = Array2(ymax - ymin + 3, xmax - xmin + 3);  
    int * shape = labels.shape();

    labels.fill(2);  // 0 for exterior, 1 for border, 2 for interior
    for (uint i=0; i < xind.size(); i++) {  // mark border
        labels(yind[i] - ymin + 1, xind[i] - xmin + 1) = 1;
    }

    // flood fill from exterior, starting point is the lower left corner
    std::queue<std::pair<int, int>> q;
    q.push(std::make_pair(0, 0));  
    while (!q.empty()) {
        std::pair<int, int> p = q.front();
        q.pop();
        int i = p.first;
        int j = p.second;
        if (labels(i, j) == 2) {  // if unvisited, mark as exterior
            labels(i, j) = 0;
            // add all neighbors to the queue
            if (i > 0) {
                q.push(std::make_pair(i - 1, j));
            }
            if (i < shape[0] - 1) {
                q.push(std::make_pair(i + 1, j));
            }
            if (j > 0) {
                q.push(std::make_pair(i, j - 1));
            }
            if (j < shape[1] - 1) {
                q.push(std::make_pair(i, j + 1));
            }
        }
    }
    
    // collect interior indices
    std::vector<int> iindx;
    std::vector<int> iindy;
    for (int i=1; i < shape[0] - 1; i++) {
        for (int j=1; j < shape[1] - 1; j++) {
            if ((labels(i, j) == 2)) {
                iindx.push_back(j + xmin - 1);
                iindy.push_back(i + ymin - 1);
            }
        }
    }
    
    return make_pair(iindy, iindx);
}



Array2 Mesh::solve(double tol, int maxiter, double omega, bool verbose){
    if (verbose) {
        std::cout << "Solving for the electric field at tol " << tol << std::endl;
    }

    // apply gauss-seidel method to solve in place for the electric field
    double error = tol + 1.;
    int iter = 0;
    double currval;
    // precompute the scheme coefficients depending on permittivity
    Array2 up_coeff(ny, nx);
    Array2 down_coeff(ny, nx);
    Array2 left_coeff(ny, nx);
    Array2 right_coeff(ny, nx);

    for (int i=1; i < ny-1; i++)
    {
        for (int j=1; j < nx-1; j++)
        {
            float center_coeff = 1./(permittivity(i,j) + permittivity(i-1,j) + 
                                     permittivity(i,j-1) + permittivity(i-1,j-1));
            right_coeff(i,j) = 0.5*(permittivity(i,j) + permittivity(i-1,j))*center_coeff;
            up_coeff(i,j) = 0.5*(permittivity(i,j) + permittivity(i,j-1))*center_coeff;
            left_coeff(i,j) = 0.5*(permittivity(i,j-1) + permittivity(i-1,j-1))*center_coeff;
            down_coeff(i,j) = 0.5*(permittivity(i-1,j) + permittivity(i-1,j-1))*center_coeff;
        }
    }



    while (error > tol && iter < maxiter)
    {
        iter++;
        if (verbose && (iter % 500 == 0))
            std::cout << "Iteration " << iter << " error: " << error << std::endl;
            
        error = 0.;
        for (int i=1; i < ny-1; i++)
        {
            for (int j=1; j < nx-1; j++)
            {
                if (node_type(i,j) == 1)  // interior node
                {
                    currval = electric_field(i,j);
                    electric_field(i,j) = ((1-omega)*currval + 
                                           omega*(electric_field(i+1,j)*up_coeff(i,j) + 
                                                  electric_field(i-1,j)*down_coeff(i,j) + 
                                                  electric_field(i,j+1)*right_coeff(i,j) + 
                                                  electric_field(i,j-1)*left_coeff(i,j)));
                    
                    error += std::abs(currval - electric_field(i,j));
                }
            }
        }
    }
    return electric_field;
}
