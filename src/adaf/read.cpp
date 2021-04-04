#include <iostream>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <string>

#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <fparameters/ParamSpaceValues.h>
#include <fparameters/SpaceIterator.h>

#include "globalVariables.h"
#include "read.h"
#include "adafFunctions.h"
#include "write.h"

void readThermalProcesses(int flags[])
{
	for (int i=0;i<numProcesses;i++) {
		flags[i] = GlobalConfig.get<int>("thermal.processNumber."+std::to_string(i));
	}
}

void readEandRParamSpace(const std::string& filename, ParamSpaceValues& data, int t, int vol)
{
	std::ifstream file;
	file.open(dataName(filename).c_str(), std::ios::in);
	
	data.ps.iterate([&](const SpaceIterator& i){
		double voll = (vol == 1) ? volume(i.val(DIM_R)) : 1.0;
		double aux1,aux2,aux3,dist;
		file >> aux1 >> aux2 >> aux3 >> dist;
		data.set(i,pow(10,dist)/voll);
	},{-1,-1,t});  
	file.close();
}

using namespace H5;

int OpenAndCloseFile() {

	const H5std_string FILE_NAME("h5tutr_dset.h5");
	const H5std_string DATASET_NAME("dset");
	const int          NX   = 4; // dataset dimensions
	const int          NY   = 6;
	const int          RANK = 2;

	// Create the data space for the dataset.
    hsize_t dims[2]; // dataset dimensions
    dims[0] = NX;
    dims[1] = NY;

	DataSpace dataspace(RANK, dims);
	//H5File file(FILE_NAME, H5F_ACC_TRUNC);

	//DataSet dataset = file.createDataSet(DATASET_NAME, PredType::STD_I32BE, dataspace);
	//file.close();
	// Data initialization.
   // int i, j;
    //int data[NX][NY]; // buffer for data to write

    //for (j = 0; j < NX; j++)
    //    for (i = 0; i < NY; i++)
    //        data[j][i] = i * 6 + j + 1;

	// Open an existing file and dataset.
    //H5File file2(FILE_NAME, H5F_ACC_RDWR);
    //DataSet dataset2 = file2.openDataSet(DATASET_NAME);

    // Write the data to the dataset using default memory space, file
    // space, and transfer properties.
    //dataset2.write(data, PredType::NATIVE_INT);
    /* Terminate access to the file. */
    //file.close();

	// Open an existing file and dataset.
	//file2.close();
    
	H5File  file3(FILE_NAME, H5F_ACC_RDWR);
    DataSet dataset = file3.openDataSet(DATASET_NAME);
	//DataSet dataset = file3.createDataSet("AVERR", PredType::STD_I32BE, dataspace);
	cout << "Opened dataset" << endl;
    // Create the data space for the attribute.
	hsize_t dims2[1]      = {2};
    DataSpace attr_dataspace = DataSpace(1, dims2);

    // Create a dataset attribute.
	const H5std_string ATTR_NAME("Units");
	StrType str_type(PredType::C_S1, H5T_VARIABLE);
	const int RANK2 = 1;
    hsize_t dims3[RANK2];
	vector<string>    attr_data;
	attr_data.push_back("erg");
	attr_data.push_back("erg s-1");
	//string * attr_data;
	//attr_data[0] = "erg";
	//attr_data[1] = "erg s-1";
    //dims3[0] = attr_data.size();  //The attribute will have 3 strings
    dims3[0] = 2;
	DataSpace att_datspc(RANK2, dims3);
    Attribute attribute = dataset.createAttribute("BOCA8", str_type, att_datspc);

	//Convert the vector into a C string array.
    //Because the input function ::write requires that.
    vector<const char *> cStrArray;
    for (int index=0; index<attr_data.size(); ++index) {
        cStrArray.push_back(attr_data[index].c_str());
    }
	vector<const char *> cStrArray2 = {"erg", "erg s-1"};

    //WRITE DATA
    //att_vector must not change during this operation
    attribute.write(str_type, (void*)&cStrArray2[0]);
    //attribute.write(str_type, attr_data[0]);
	//file3.close();
	return 0;
}
