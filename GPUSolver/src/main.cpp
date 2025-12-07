#include <iostream>

#include <string>

#include "include/IO.hpp"
#include "include/model.hpp"
#include "include/ComputeModel.hpp"


int main(int argc, char* argv[])  {

    using index_type = int;
    using value_type = float;

    std::string model_name{};

    Model poissfem_model;
    CSR_matrix<index_type, value_type> A;
    std::vector<value_type> b;

    DisplayWelcomeHeader();

    model_name = (argc > 1 ? argv[1] : "Sphere_00");

    ReadVTK(model_name, poissfem_model);

    std::cout << "Model: " << model_name << std::endl;
    std::cout << "Number of vertices: " << poissfem_model.n_vertices << std::endl;
    std::cout << "Number of elements: " << poissfem_model.n_elements << std::endl;
    assignBoundaryConditions(poissfem_model);


    initialize_CSR_indices<index_type, value_type>(poissfem_model, A);
    fill_FEM_CSR<index_type, value_type>(poissfem_model, A, b);

    write_vtu(model_name, poissfem_model);

    return 0;

}