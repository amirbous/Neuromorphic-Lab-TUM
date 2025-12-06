#include <iostream>

#include <string>

#include "include/IO.hpp"
#include "include/geometry.hpp"
#include "include/ComputeModel.hpp"


int main(int argc, char* argv[])  {

    std::string model_name{};

    ProblemProperties* poissfem_properties = new ProblemProperties();
    std::vector<Vertex> vertices;
    std::vector<Element> conn;


    DisplayWelcomeHeader();

    model_name = (argc > 1 ? argv[1] : "Sphere");

    ReadVTK(model_name, vertices, conn, *poissfem_properties);

    std::cout << "Model: " << model_name << std::endl;
    std::cout << "Number of vertices: " << poissfem_properties->n_vertices << std::endl;
    std::cout << "Number of elements: " << poissfem_properties->n_faces << std::endl;

    assignBoundaryConditions(vertices, conn, *poissfem_properties);

    write_vtu(model_name, vertices, conn, *poissfem_properties);

    return 0;

}