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

    bool solution_exists = false;

    std::string model_name{};

    Model<index_type, value_type> poissfem_model;
    CSR_matrix<index_type, value_type> A;

    std::vector<value_type> b; // we still don't know the size, given the sparsity pattern but will be decided in initialize_CSR_indices

    DisplayWelcomeHeader();

    model_name = (argc > 1 ? argv[1] : "Sphere_00");

    ReadVTK(model_name, poissfem_model);

    std::cout << "Model: " << model_name << std::endl;
    std::cout << "Number of vertices: " << poissfem_model.n_vertices << std::endl;
    std::cout << "Number of elements: " << poissfem_model.n_elements << std::endl;

    initialize_CSR_indices<index_type, value_type>(poissfem_model, A);
    fill_FEM_CSR<index_type, value_type>(poissfem_model, A, b);

    if (!solution_exists) {
        WriteCSRMatrix<index_type, value_type>(A, model_name);
        WriteRHSVector<index_type, value_type>(b, model_name);
    }

    // TODO implement solver ******************
    
    // *********************************************



    if (solution_exists) {

        std::vector<value_type> x_solution = ReadVector<index_type, value_type>(model_name, "sol");

        // Update model with results and save
        for (int i = 0; i < poissfem_model.n_vertices; ++i) {
            poissfem_model.vertices[i].potential = x_solution[i];
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
        value_type max_vol{0.0};
        // total volume for average volume calculation
        value_type total_vol{0.0};
        // average volume for convergence check
        value_type avg_vol{0.0}; 


        avg_vol= total_vol / poissfem_model.n_elements;
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

            if (vol > max_vol) max_vol = vol;
            total_vol += vol;

            // centroid for error at centroid
            value_type cx = (x1 + x2 + x3 + x4) * 0.25f;
            value_type cy = (y1 + y2 + y3 + y4) * 0.25f;
            value_type cz = (z1 + z2 + z3 + z4) * 0.25f;

            // average solution to get FEM solution at centroid
            value_type u_fem_centroid = (x_solution[v1] + x_solution[v2] + x_solution[v3] + x_solution[v4]) * 0.25f;
            
            // Exact solution at centroid
            value_type u_exact_centroid = analytical_solution(cx, cy, cz);

            // error at centroid
            value_type diff = u_fem_centroid - u_exact_centroid;
            // accumulate weighted squared error
            total_l2_error_sq += (diff * diff) * vol;
        }

        // final l2 error
        value_type l2_error = std::sqrt(total_l2_error_sq);
        avg_vol= total_vol / poissfem_model.n_elements;



        std::cout << "Max volume: " << max_vol << std::endl;
        std::cout << "Avg volume: " << avg_vol << std::endl;
        std::cout << "L2 Norm (Integral): " << l2_error << std::endl;
        // write geometry solution

        write_vtu<index_type, value_type>(model_name + "_solution", poissfem_model);
    }
    else {
        std::cout << "No solution available." << std::endl;
    }

    return 0;
}
