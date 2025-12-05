#include <iostream>

#include <string>

#include "include/IO.hpp"


int main(int argc, char* argv[])  {

    std::string model_name{};

    ProblemProperties* electrostatics_problem_properties = new ProblemProperties();
    std::vector<Vertex> vertices;
    std::vector<Face> faces;


    DisplayWelcomeHeader();

    model_name = (argc > 1 ? argv[1] : "standard_model");

    ReadVTK(model_name, vertices, faces, *electrostatics_problem_properties);

    return 0;

}