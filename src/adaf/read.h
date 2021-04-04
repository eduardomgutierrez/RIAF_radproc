#pragma once

#include <H5Cpp.h>
#include <fparameters/ParamSpaceValues.h>

void readThermalProcesses(int []);
void readEandRParamSpace(const std::string& filename, ParamSpaceValues& data, int t, int vol);
int OpenAndCloseFile();