#ifndef IO_HPP
#define IO_HPP

#include "include/geometry.hpp"

#include <stdio.h>
#include <string>
#include <vector>

void DisplayWelcomeHeader();
void ReadVTK(std::string model_name, std::vector<Vertex>& vertices, std::vector<Element>& conn,
				ProblemProperties& problem_properties);
void ReadObj(std::string model_name, std::vector<Vertex>& vertices, std::vector<Element>& conn,
				ProblemProperties& problem_properties);

void write_vtu(const std::string problem_name, const std::vector<Vertex>& vertices, const std::vector<Element>& conn,
														ProblemProperties& problem_properties);

#endif



