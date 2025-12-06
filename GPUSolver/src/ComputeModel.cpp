#include <vector>
#include <string>
#include <fstream>

#include "math.h"

#include "include/geometry.hpp"

void assignBoundaryConditions(std::vector<Vertex> &vertices, const std::vector<Element> &conn, ProblemProperties &problem_properties) {
    for (int i = 0; i < problem_properties.n_vertices; ++i) {
        vertices[i].potential = sqrt(vertices[i].x * vertices[i].x + vertices[i].y * vertices[i].y + vertices[i].z * vertices[i].z);
    }


}
