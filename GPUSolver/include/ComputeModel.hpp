#include <vector>
#include <string>
#include <fstream>

#include "include/geometry.hpp"

void assignBoundaryConditions(std::vector<Vertex> &vertices, const std::vector<Element> &conn, ProblemProperties &problem_properties);
