#include "adafParameters.h"
#include "adafFunctions.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <fmath/physics.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
#include <fparameters/SpaceIterator.h>
#include "State.h"

using namespace std;

///////////////////////

void adafParameters() 
{	
	ifstream adafFile, adafParams;

	adafFile.open("adafFile.txt"); adafParams.open("adafParameters.txt");
	adafFile >> nRaux;
	adafParams >> blackHoleMass >> accRateOut >> s >> magFieldPar >> alpha >> jAngMom >> delta >> rTr >> powerIndex;
	adafParams.close();
	
	logr.resize(nRaux,0.0);
	logTi.resize(nRaux,0.0);
	logTe.resize(nRaux,0.0);
	logv.resize(nRaux,0.0);
	
	for (size_t i=0;i<nRaux;i++) {
		adafFile >> logr[i] >> logTi[i] >> logTe[i] >> logv[i];
	}
	adafFile.close();
    
	blackHoleMass *= solarMass;
	schwRadius = 2.0*gravitationalConstant*blackHoleMass / cLight2;
	accRateOut = accRateOut * 1.39e18 * blackHoleMass/solarMass;
	cout << "Total outer accretion power = " << accRateOut*cLight2 << " erg s^{-1}" << endl;
	cout << "Total inner accretion power = " << accRateOut*pow(1.0/exp(logr.back()), s)*cLight2 
		 << " erg s^{-1}" << endl;
    
	rTr = rTr * schwRadius;
	rOutCD = GlobalConfig.get<double>("rOutCD") * schwRadius;
	eMeanMolecularWeight = GlobalConfig.get<double>("mu_e");
	iMeanMolecularWeight = GlobalConfig.get<double>("mu_i");

	nR = GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
	nE = GlobalConfig.get<int>("model.particle.photon.dim.energy.samples");
	nRcd = GlobalConfig.get<int>("model.particle.default.dim.radius_cd.samples");
    
    paso_r = pow(exp(logr.back())/exp(logr.front()),1.0/(nR));
	paso_rCD = pow(rOutCD/rTr,1.0/(nRcd));

	logMinEnergy = GlobalConfig.get<double>("model.particle.photon.dim.energy.min");
	logMaxEnergy = GlobalConfig.get<double>("model.particle.photon.dim.energy.max");

	inclination = GlobalConfig.get<double>("inclination");
    
    calculateComptonScatt = GlobalConfig.get<int>("calculateComptonScatt");
	height_method = GlobalConfig.get<int>("height_method");
    
    calculateThermal = GlobalConfig.get<int>("calculateThermal");
	numProcesses = GlobalConfig.get<int>("thermal.numProcesses");
    if (calculateThermal) {
        calculateComptonRedMatrix = GlobalConfig.get<int>("thermal.compton.calculateRedMatrix");
        comptonMethod = GlobalConfig.get<int>("thermal.compton.method");
        if (1) {

            nGammaCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nGammaCompton");
            nTempCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nTempCompton");
            nNuPrimCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nNuPrimCompton");
            nNuCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nNuCompton");
            gammaMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.gammaMinCompton");
            gammaMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.gammaMaxCompton");
            tempMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.tempMinCompton");
            tempMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.tempMaxCompton");
            nuPrimMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuPrimMinCompton");
            nuPrimMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuPrimMaxCompton");
            nuMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuMinCompton");
            nuMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuMaxCompton");
            
            ofstream fileSizes;
            fileSizes.open("comptonRedMatrix/sizesVecCompton.dat",ios::out);
            fileSizes << nGammaCompton << "\t" << gammaMinCompton << "\t" << gammaMaxCompton << endl
                      << nTempCompton << "\t" << tempMinCompton << "\t" << tempMaxCompton << endl
                      << nNuPrimCompton << "\t" << nuPrimMinCompton << "\t" << nuPrimMaxCompton << endl
                      << nNuCompton << "\t" << nuMinCompton << "\t" << nuMaxCompton << endl;
            fileSizes.close();
        } else {
            ifstream fileSizes;
            fileSizes.open("comptonRedMatrix/sizesVecCompton.dat",ios::in);
            fileSizes >> nGammaCompton >> gammaMinCompton >> gammaMaxCompton
                      >> nTempCompton >> tempMinCompton >> tempMaxCompton
                      >> nNuPrimCompton >> nuPrimMinCompton >> nuPrimMaxCompton
                      >> nNuCompton >> nuMinCompton >> nuMaxCompton;
            fileSizes.close();
        }
        calculatePhotonDensityGap = GlobalConfig.get<int>("thermal.calculatePhotonDensityGap");
    }
	calculateJetEmission = GlobalConfig.get<int>("nonThermal.calculateJetEmission");
	
    calculateNonThermal = GlobalConfig.get<int>("calculateNonThermal");
    if (calculateNonThermal) {
		accMethod = GlobalConfig.get<double>("nonThermal.acc_method");
		calculateNTprotons = GlobalConfig.get<int>("nonThermal.protons");
		calculateNTelectrons = GlobalConfig.get<int>("nonThermal.electrons");
        calculateLosses = GlobalConfig.get<int>("nonThermal.calculateLosses");
		calculateFlare = GlobalConfig.get<int>("nonThermal.calculateFlare");
        calculateNTdistributions = GlobalConfig.get<int>("nonThermal.calculateDistributions");
        calculateNonThermalLum = GlobalConfig.get<int>("nonThermal.calculateLuminosities");
		calculateNonThermalHE = GlobalConfig.get<int>("nonThermal.calculateHighEnergyProcesses");
        calculateNeutronInj = GlobalConfig.get<int>("nonThermal.neutrons.calculateInjection");
		calculateNeutronDis = GlobalConfig.get<int>("nonThermal.neutrons.calculatePropagation");
		calculateJetDecay = GlobalConfig.get<int>("nonThermal.neutrons.calculateJetDecay");
		calculateSecondaries = GlobalConfig.get<int>("nonThermal.calculateSecondaries");
		calculateNeutrinos = GlobalConfig.get<int>("nonThermal.calculateNeutrinos");
	}
}