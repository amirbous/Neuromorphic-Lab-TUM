#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <numeric>
#include <algorithm>


#include "include/IO.hpp"
#include "include/model.hpp"
#include "include/ComputeModel.hpp"

int main(int argc, char* argv[])  {

    using index_type = int;
    using value_type = float;

    bool write_matrix = true;

    std::string model_name{};

    Model<index_type, value_type> poissfem_model;
    CSR_matrix<index_type, value_type> A;

    std::vector<value_type> b; // we still don't know the size, given the sparsity pattern but will be decided in initialize_CSR_indices


    model_name = (argc > 1 ? argv[1] : "Sphere_00");
    write_matrix = (argc > 2 ? (std::string(argv[2]) == "1" ? true : false) : false);

    ReadVTK(model_name, poissfem_model);

    initialize_CSR_indices<index_type, value_type>(poissfem_model, A);
    fill_FEM_CSR<index_type, value_type>(poissfem_model, A, b);

    if (write_matrix) {
        WriteCSRMatrix<index_type, value_type>(A, model_name);
        WriteRHSVector<index_type, value_type>(b, model_name);
    }

    // TODO implement solver ******************
    
    // *********************************************

    if (!write_matrix) {

        std::vector<value_type> x_solution = ReadVector<index_type, value_type>(model_name, "sol");
        std::vector<index_type> boundary_nodes = extract_boundary_nodes<index_type, value_type>(poissfem_model);
        
        // Create a fast lookup mask
        std::vector<bool> is_boundary(poissfem_model.n_vertices, false);
        for (index_type idx : boundary_nodes) {
            is_boundary[idx] = true;
        }

        // 2. Assign potentials
        // We maintain a separate counter for the solution vector, which tracks internal nodes only.
        index_type internal_idx = 0; 

        for (int i = 0; i < poissfem_model.n_vertices; ++i) {
            if (!is_boundary[i]) {
                // CASE A: Internal Node -> Value comes from the solver
                poissfem_model.vertices[i].potential = x_solution[internal_idx];
                internal_idx++;
            } else {
                // CASE B: Boundary Node -> Value comes from the known boundary condition
                // (The solver did not calculate this, so we must re-calculate it here)
                value_type bx = poissfem_model.vertices[i].x;
                value_type by = poissfem_model.vertices[i].y;
                value_type bz = poissfem_model.vertices[i].z;
                
                poissfem_model.vertices[i].potential = analytical_solution(bx, by, bz);
            }
        }

        //assign analytical solution to densities
        for (int i = 0; i < poissfem_model.n_vertices; ++i) {
            value_type x = poissfem_model.vertices[i].x;
            value_type y = poissfem_model.vertices[i].y;
            value_type z = poissfem_model.vertices[i].z;
            poissfem_model.vertices[i].density = analytical_solution<value_type>(x, y, z);
        }

        //analytical solution 
        std::vector ana_sol = std::vector<value_type>(poissfem_model.n_vertices);
        for (int i = 0; i < poissfem_model.n_vertices; ++i) {
            value_type x = poissfem_model.vertices[i].x;
            value_type y = poissfem_model.vertices[i].y;
            value_type z = poissfem_model.vertices[i].z;
            ana_sol[i] = analytical_solution<value_type>(x, y, z);
        }

        // l2 norm error
        value_type total_l2_error_sq{0.0};
        // max tetrahedra volume for convergence check
        value_type max_edge_length{0.0};

        std::vector<Edge<index_type>> edges = get_mesh_edges<index_type, value_type>(poissfem_model);
        max_edge_length = compute_max_edge_length<index_type, value_type>(poissfem_model, edges);

        for (int i = 0; i < poissfem_model.n_elements; ++i) {

            index_type v1 = poissfem_model.elements[i].v1;
            index_type v2 = poissfem_model.elements[i].v2;
            index_type v3 = poissfem_model.elements[i].v3;
            index_type v4 = poissfem_model.elements[i].v4;

            value_type x1 = poissfem_model.vertices[v1].x, y1 = poissfem_model.vertices[v1].y, z1 = poissfem_model.vertices[v1].z;
            value_type x2 = poissfem_model.vertices[v2].x, y2 = poissfem_model.vertices[v2].y, z2 = poissfem_model.vertices[v2].z;
            value_type x3 = poissfem_model.vertices[v3].x, y3 = poissfem_model.vertices[v3].y, z3 = poissfem_model.vertices[v3].z;
            value_type x4 = poissfem_model.vertices[v4].x, y4 = poissfem_model.vertices[v4].y, z4 = poissfem_model.vertices[v4].z;

            // have to move this to a seperate method
            value_type vol = std::abs(
                (x2 - x1) * ((y3 - y1) * (z4 - z1) - (y4 - y1) * (z3 - z1)) -
                (x3 - x1) * ((y2 - y1) * (z4 - z1) - (y4 - y1) * (z2 - z1)) +
                (x4 - x1) * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))
            ) / 6.0f;


            // centroid for error at centroid
            value_type cx = (x1 + x2 + x3 + x4) * 0.25f;
            value_type cy = (y1 + y2 + y3 + y4) * 0.25f;
            value_type cz = (z1 + z2 + z3 + z4) * 0.25f;

            // average solution to get FEM solution at centroid
            value_type u_fem_centroid = (poissfem_model.vertices[v1].potential + poissfem_model.vertices[v2].potential + 
                                            poissfem_model.vertices[v3].potential + poissfem_model.vertices[v4].potential) * 0.25f;
            
            // Exact solution at centroid
            value_type u_exact_centroid = analytical_solution(cx, cy, cz);

            // error at centroid
            value_type diff = u_fem_centroid - u_exact_centroid;
            // accumulate weighted squared error
            total_l2_error_sq += (diff * diff) * vol;
        }

        // final l2 error
        value_type l2_error = std::sqrt(total_l2_error_sq);



        print_log<index_type, value_type>(model_name, poissfem_model, A, max_edge_length, l2_error, "");

        write_vtu<index_type, value_type>(model_name + "_solution", poissfem_model);
    }


    return 0;
}
