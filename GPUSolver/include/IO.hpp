#ifndef IO_HPP
#define IO_HPP

#include "include/model.hpp"

#include <stdio.h>
#include <string>
#include <vector>

void DisplayWelcomeHeader();
void ReadVTK(std::string model_name, Model &model);


void write_vtu(const std::string problem_name,
               Model &model);

#endif



