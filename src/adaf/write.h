#pragma once

// External libraries
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <H5Cpp.h>

// Local libraries
#include <fmath/physics.h>
#include <fparticle/Particle.h>

// Local headers
#include "State.h"

class File;

typedef std::map<std::string, File*> OFM;

double safeLog10(double x);

class File {
public:
	std::string name;
	std::ofstream file;
	File(std::vector<File*>& all, OFM& byName, std::string name) :name(name) {
		file.open(name + ".txt");
		byName[name] = this;
		all.push_back(this);
	}
	void check() {
		std::ifstream result;
		result.open(name + ".txt");

		std::ifstream reference;
		reference.open(name + ".ref");

		char a, b;
		int c = 0;
		while (!result.eof()) {
			c++;
			result >> a;
			reference >> b;
			if (a != b) {
				std::cout << name << +" doesn't match reference." << std::endl;
				throw;
			}
		}
		if (!reference.eof()) {
			std::cout << name << +" reference content larger than result." << std::endl;
			throw;
		}
		result.close();
		reference.close();
		std::cout << name << +" matches reference (" << c << "b)." << std::endl;
	}
};

/*It writes data for all param space*/
void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data);

/*It writes data for all param space, but fixing r*/
void writeEandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int r); 

/*Idem previous, but fixing E*/
void writeRandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int E);

/*Idem previous, but fixing t*/
void writeEandRParamSpace(const std::string& filename, const ParamSpaceValues& data, int t, int vol);

/* Idem previous, but integrating in energy */
void writeRParamSpace(const std::string& filename, const ParamSpaceValues& data, int t, int s);

/*It writes data as a function of energy for a given r and t*/
void writeEnergyFunction(const std::string& filename, const ParamSpaceValues& data, int r, int t);


void writeMatrix(const std::string& filename, Particle& p, Matrix& a);
void writeFields(State& st);
std::string dataName(std::string);
int createH5file(const H5std_string FILE_NAME, const State& st);