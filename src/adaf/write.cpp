// External libraries
#include <boost/property_tree/ptree.hpp>

// Local libraries
#include <fmath/RungeKutta.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

// Local headers
#include "write.h"
#include "State.h"
#include "modelParameters.h"
#include "globalVariables.h"
#include "adafFunctions.h"

//namespace {
inline double safeLog10( double x )	{ return x > 0.0 ? log10(x) : -100.0; }
//}
std::string dataName(std::string id) {
	return id + ".dat";
}

void generateViewScript(std::string path) {
	std::string filename = path.substr(path.find("\\") + 1);
	std::string folder = path.substr(0,path.find("\\"));
	std::ofstream file;
	file.open((folder+"/plots/plot-"+filename+".bat").c_str(), std::ios::out);
	file << "@../../plot-svg-and-view.bat " + filename;
	file.close();
}

void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);
	//const ParamSpace* a = &(data.ps);
	
	// a.iterate;  no me deja hacer esta operacion
	
	data.ps.iterate([&file, &data](const SpaceIterator& i){
		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double logR = (i.val(DIM_R));
		double logT = (i.val(DIM_Rcd));
		double logQ = safeLog10(data.get(i)); //log10(salida.values(i));  // values(i));
//		salida.values(i);


		file << logE << '\t' << logR << '\t' << logT << '\t' << 
			logQ << std::endl;
			//logQ << std::endl;
	});

	file.close();
	generateViewScript(filename);
}

void writeEandRParamSpace(const std::string& filename, const ParamSpaceValues& data, int t,int vol)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	//double time = log10(data.ps[1][t]);
	
	data.ps.iterate([&](const SpaceIterator& i) {
		
		double logE = log10(i.val(DIM_E) / EV_TO_ERG);
		double r = i.val(DIM_R);
		double voll = (vol == 1) ? volume(r) : 1.0;
		double logQ = safeLog10(data.get(i)*voll);
		file << logE << '\t' << i.coord[DIM_R] << '\t' << safeLog10(r/schwRadius) << '\t' << logQ << std::endl;
			
	}, { -1, -1, t });  

	file.close();
	generateViewScript(filename);
}

void writeRParamSpace(const std::string& filename, const ParamSpaceValues& data, int t, int s)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);
	double Emin = data.ps[DIM_E].first();
	double Emax = data.ps[DIM_E].last();
	data.ps.iterate([&](const SpaceIterator& iR) {
		double tot = integSimpson(log(Emin),log(Emax),
						[&](double loge)
						{
							double e = exp(loge);
							double nPh = data.interpolate({{DIM_E,e}},&iR.coord);
							return e*nPh;
						},100);
		file << safeLog10(iR.val(DIM_R)/schwRadius) << '\t' << safeLog10(tot) << std::endl;
	},{t,-1,s});  
	file.close();
	generateViewScript(filename);
}

void writeEandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int r)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	
	file << "log(r)=" << logR << '\t' ;

	for (size_t t_ix = 0; t_ix < data.ps[2].size(); t_ix++) {
		double time = data.ps[2][t_ix];
		file << "t=" << log10(time) << '\t';
	}


	for (int E_ix = 0; E_ix < data.ps[0].size(); E_ix++) {

		file << std::endl;

		double logE = log10(data.ps[0][E_ix]/1.6e-12);

		file << logE << '\t';

		data.ps.iterate([&file, &data](const SpaceIterator& i){

			//double logR = log10(i.val(DIM_R));
			//double time = i.val(DIM_T);
			double logQ = safeLog10(data.get(i));

			file << logQ << '\t';
			;
		}, { E_ix, r, -1 });  //el -1 indica que las E se recorren, no quedan fijas
		//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )
	}
	file.close();
	generateViewScript(filename);
}


void writeRandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int E)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);


	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double r = i.val(DIM_R);
		//double theta = i.val(DIM_THETA);
		
		file << r << '\t' << '\t' << data.get(i) << std::endl;
		
	}, { E, -1, -1 });  //el -1 indica que las E se recorren, no quedan fijas
	//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )
	
	file.close();
	generateViewScript(filename);
}


void writeEnergyFunction(const std::string& filename, const ParamSpaceValues& data, int r, int t)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	// version larga
	double logT = log10(data.ps.dimensions[2]->values[t]);

	file << "log(r)=" << logR << '\t' << "log(t)=" << logT << std::endl;
	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double logQ = safeLog10(data.get(i));

		file << logE << '\t' << logQ << std::endl;
		;
	}, { -1, r, t });  //el -1 indica que las E se recorren, no quedan fijas
	//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )

	file.close();
	generateViewScript(filename);
}




void writeMatrix(const std::string& filename, Particle& p, Matrix& a)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	int nR = p.ps[DIM_R].size();  //ver el -1

	file << '\t';
	
	for (size_t z_j = 0; z_j < nR; z_j++) { 
		const double r_j = p.ps[DIM_R][z_j];
		file << r_j << '\t' ;
	}
	
	file << std::endl; 
			
	for (size_t z_i = 0; z_i < nR; z_i++) { 
		
		const double r_i = p.ps[DIM_R][z_i];
				
		file << r_i << '\t';
		
		for (size_t z_j = 0; z_j < nR; z_j++) { 
			const double r_j = p.ps[DIM_R][z_j];
	
			file << a[z_i][z_j] << '\t' ;
		}
		
		
		
		file << std::endl;  
	}
	
	file.close();
	generateViewScript(filename);
}

void writeFields(State& st) {
	std::ofstream fields;
	fields.open("fields.dat",std::ios::out);
	fields		<< "r_in [Rs] = " << exp(logr.front()) << endl;
	fields		<< "r_out [Rs] = " << exp(logr.back()) << endl << endl;
	fields  << "rB1 [Rs]" 	<< "\t"
			<< "r [Rs]" 	<< "\t"
			<< "rB2 [Rs]" 	<< "\t"
			<< "MdotRIAF" 	<< "\t"
			<< "MdotCD" 	<< "\t"
			<< "Te"			<< "\t"
			<< "Ti"			<< "\t"
			<< "H/R"		<< "\t"
			<< "ne"			<< "\t"
			<< "v/c"		<< "\t"
			<< "B"			<< endl;

	st.photon.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		fields 	<< r/schwRadius/sqrt(paso_r) << "\t"
				<< r/schwRadius << "\t"
				<< r/schwRadius * sqrt(paso_r) << "\t"
				<< gAcc(r) << "\t"
				<< accRateColdDisk(r)/accRateOut << "\t"
				<< st.tempElectrons.get(iR) << "\t"
				<< st.tempIons.get(iR) << "\t"
				<< height_fun(r)/r << "\t"
				<< electronDensity(r) << "\t"
				<< abs(radialVel(r))/cLight <<"\t"
				<< st.magf.get(iR) << endl;
	},{0,-1,0});
	fields.close();
	
	std::ofstream SSD;
	SSD.open("SSD.dat",std::ios::out);
	SSD		<< "r_tr [Rs] = " << rTr/schwRadius << endl;
	SSD		<< "r_outCD [Rs] = " << rOutCD/schwRadius << endl << endl;
	SSD  	<< "rB1 [Rs]"			<< "\t"
			<< "r [Rs]" 			<< "\t"
			<< "rB2 [Rs]"			<< "\t"
			<< "MdotCD/MdotOut" 	<< endl;
	
	st.photon.ps.iterate([&](const SpaceIterator& iRcd) {
		double rCD = iRcd.val(DIM_Rcd);
		SSD 	<< rCD/sqrt(paso_rCD)/schwRadius << "\t"
				<< rCD/schwRadius << "\t"
				<< rCD*sqrt(paso_rCD)/schwRadius << "\t"
				<< accRateColdDisk(rCD)/accRateOut << endl;
	},{0,0,-1});
	SSD.close();
}

using namespace H5;

float *toFloatArray(const double *arr, size_t n) {
  float *ret = new float[n];
  for (int i = 0; i < n; i++) {
    ret[i] = float(arr[i]);
  }
  return ret;
}

void MultiplyVectorByScalar(Vector v, double s, Vector& w) {
	// Multiply vector v by scalar s and stores the data in w
    transform(v.begin(), v.end(), w.begin(), [s](double &x){ return x*s; });
}

void DivideVectorByVector(Vector v1, Vector v2, Vector& w) {
	// Multiply vector v by scalar s and stores the data in w
    transform(v1.begin(), v1.end(), v2.begin(), w.begin(), divides<double>());
}

void EvaluateFunctionInVector(fun1 f, Vector v, Vector& w) {
	// Evaluate f in vector v and stores the data in w
    transform(v.begin(), v.end(), w.begin(), [f](double &x){ return f(x); });
}

void myH5_write_single_float(const double x, H5File *h5f, H5std_string DATASET_NAME)
{
	float y[1];	*y = float(x);
	hsize_t dimsZero[1] = { 1 };
	DataSpace *dataspace = new DataSpace(1, dimsZero);
	DataSet *dataset = new DataSet(h5f->createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, *dataspace));
	dataset->write(y, PredType::NATIVE_FLOAT);
	delete dataspace;	delete dataset;
}

void myH5_write_single_float(const double x, Group *h5g, H5std_string DATASET_NAME)
{
	float y[1];	*y = float(x);
	hsize_t dimsZero[1] = { 1 };
	DataSpace *dataspace = new DataSpace(1, dimsZero);
	DataSet *dataset = new DataSet(h5g->createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, *dataspace));
	dataset->write(y, PredType::NATIVE_FLOAT);
	delete dataspace;	delete dataset;
}

void myH5_write_1d_float(const Vector v, H5File *h5f, H5std_string DATASET_NAME)
{
	size_t n = v.size();
	hsize_t dims[1] = { n };
	DataSpace *dataspace = new DataSpace(1, dims);
	DataSet *dataset = new DataSet(h5f->createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, *dataspace));
	dataset->write(toFloatArray(v.data(), n), PredType::NATIVE_FLOAT);
	delete dataspace;	delete dataset;
}

void myH5_write_1d_float(const Vector v, Group *h5g, H5std_string DATASET_NAME)
{
	size_t n = v.size();
	hsize_t dims[1] = { n };
	DataSpace *dataspace = new DataSpace(1, dims);
	DataSet *dataset = new DataSet(h5g->createDataSet(DATASET_NAME, PredType::NATIVE_FLOAT, *dataspace));
	dataset->write(toFloatArray(v.data(), n), PredType::NATIVE_FLOAT);
	delete dataspace;	delete dataset;
}

void myH5_write_single_double(const double x, H5File *h5f, H5std_string DATASET_NAME)
{
	float y[1];	*y = float(x);
	hsize_t dimsZero[1] = { 1 };
	DataSpace *dataspace = new DataSpace(1, dimsZero);
	DataSet *dataset = new DataSet(h5f->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace));
	dataset->write(y, PredType::NATIVE_DOUBLE);
	delete dataspace;	delete dataset;
}

void myH5_write_single_double(const double x, Group *h5g, H5std_string DATASET_NAME)
{
	double y[1];	*y = double(x);
	hsize_t dimsZero[1] = { 1 };
	DataSpace *dataspace = new DataSpace(1, dimsZero);
	DataSet *dataset = new DataSet(h5g->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace));
	dataset->write(y, PredType::NATIVE_DOUBLE);
	delete dataspace;	delete dataset;
}

void myH5_write_1d_double(const Vector v, H5File *h5f, H5std_string DATASET_NAME)
{
	size_t n = v.size();
	hsize_t dims[1] = { n };
	DataSpace *dataspace = new DataSpace(1, dims);
	DataSet *dataset = new DataSet(h5f->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace));
	dataset->write(v.data(), PredType::NATIVE_DOUBLE);
	delete dataspace;	delete dataset;
}

void myH5_write_1d_double(const Vector v, Group *h5g, H5std_string DATASET_NAME)
{
	size_t n = v.size();
	hsize_t dims[1] = { n };
	DataSpace *dataspace = new DataSpace(1, dims);
	DataSet *dataset = new DataSet(h5g->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace));
	dataset->write(v.data(), PredType::NATIVE_DOUBLE);
	delete dataspace;	delete dataset;
}

int createH5file(const H5std_string FILE_NAME, const State& st)
{
	// Try block to detect exceptions raised by any of the calls inside it
    try {
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        Exception::dontPrint();

        // Create a new file using default property lists.
        H5File *h5file = new H5File(FILE_NAME, H5F_ACC_TRUNC);

		// HEADER ///////////////////////////////////////////////////////////////////////
		Group *header_id = new Group(h5file->createGroup("Header"));

		myH5_write_single_float(blackHoleMass/solarMass, header_id, "M");
		myH5_write_single_float(accRateOut / (1.39e18 * blackHoleMass/solarMass), header_id, "Mdot_out");
		myH5_write_single_float(s, header_id, "sWind");
		myH5_write_single_float(delta, header_id, "delta_e");
		myH5_write_single_float(magFieldPar, header_id, "beta");
		myH5_write_single_float(alpha, header_id, "alpha");
		myH5_write_single_float(jAngMom, header_id, "jAngMom");

		delete header_id;
		/////////////////////////////////////////////////////////////////////////////////

		// RIAF /////////////////////////////////////////////////////////////////////////

		Group *RIAF_id = new Group(h5file->createGroup("RIAF"));

		Vector rVec = st.photon.ps.dimensions[DIM_R]->values;
		Vector w(nR, 0.0);
	
		MultiplyVectorByScalar(rVec, 1.0/schwRadius, w);
		myH5_write_1d_float(w, RIAF_id, "r");
		w.assign(w.size(), 0.0);

		EvaluateFunctionInVector(volume, rVec, w);
		myH5_write_1d_double(w, RIAF_id, "volume");
		w.assign(w.size(), 0.0);

		EvaluateFunctionInVector([](double r) { return accRateADAF(r)/accRateOut; }, rVec, w);
		myH5_write_1d_float(w, RIAF_id, "Mdot");
		w.assign(w.size(), 0.0);
		
		EvaluateFunctionInVector([](double r)
					{ return boltzmann/electronRestEnergy * electronTemp(r); }, rVec, w);
		myH5_write_1d_float(w, RIAF_id, "Te");

		EvaluateFunctionInVector([](double r)
					{ return boltzmann/protonRestEnergy * ionTemp(r); }, rVec, w);
		myH5_write_1d_float(w, RIAF_id, "Ti");
		w.assign(w.size(), 0.0);

		EvaluateFunctionInVector([](double r) { return height_fun(r)/r; }, rVec, w);
		myH5_write_1d_float(w, RIAF_id, "H_R");
		w.assign(w.size(), 0.0);

		EvaluateFunctionInVector(electronDensity, rVec, w);
		myH5_write_1d_float(w, RIAF_id, "n_e");
		w.assign(w.size(), 0.0);

		EvaluateFunctionInVector([](double r) {return abs(radialVel(r))/cLight; }, rVec, w);
		myH5_write_1d_float(w, RIAF_id, "v_c");
		w.assign(w.size(), 0.0);

		EvaluateFunctionInVector(magneticField, rVec, w);
		myH5_write_1d_float(w, RIAF_id, "B");

		myH5_write_1d_float(redshift_to_inf, RIAF_id, "redshift_to_inf");

		RIAF_id->close();

		h5file->createGroup("SSD");

		for (const auto &p : st.particles) {
   			h5file->createGroup(p->id);
		}
		h5file->close();

    } // end of try block
    // catch failure caused by the H5File operations
    catch (FileIException error) {
        error.printErrorStack();
        return -1;
    }
    // catch failure caused by the Group operations
    catch (GroupIException error) {
        error.printErrorStack();
        return -1;
    }
	catch (DataSetIException error) {
        error.printErrorStack();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch (DataSpaceIException error) {
        error.printErrorStack();
        return -1;
    }
	return 0;
}