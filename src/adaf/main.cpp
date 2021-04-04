// External libraries
#include <iostream>
#include <stdio.h>
#include <boost/property_tree/ptree.hpp>

// Local libraries
#include <inout/ioutil.h>
#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>

// Local headers
#include "read.h"
#include "messages.h"
#include "modelParameters.h"
#include "State.h"
#include "write.h"
#include "comptonScattMatrix.h"
#include "thermalProcesses.h"
#include "adafFunctions.h"
#include "globalVariables.h"
#include "flareProcesses.h"
#include "redshiftFunction.h"
#include "thermalDistribution.h"
#include "NTtimescales.h"
#include "NTinjection.h"
#include "injectionNeutrons.h"
#include "distributionNeutrons.h"
#include "NTdistribution.h"
#include "NTradiation.h"
#include "absorption.h"
#include "secondariesProcesses.h"
#include "jetEmission.h"

using namespace std;

int main()
{
	string folder{ prepareOutputfolder() };
	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		
		State model(GlobalConfig.get_child("model"));
		redshiftFactor(model);
		createH5file("h5prueba.h5", model);
	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
	}
	return 0;
}

/*	string folder{ prepareOutputfolder() };
	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		
		State model(GlobalConfig.get_child("model"));
		writeFields(model);
		redshiftFactor(model);

		if (calculateThermal) {
			if (calculateComptonScatt)
				comptonScattMatrix(model);
			else
				comptonScattMatrixRead(model);
			thermalRadiation(model, "lumThermal.dat");
		}
		else
			readEandRParamSpace("photonDensity", model.photon.distribution, 0, 0);
		
		//writeEandRParamSpace("photonDensity_z",model.photon.injection,0,0);
		
		if (calculateJetEmission)
			jetProcesses(model,"lumJet.txt");
		
		if (calculateFlare)
			flareProcesses(model);
		
//***********nonthermal particles**************		
		if (calculateNonThermal) {
            
			if (calculateLosses) {
				if (calculateNTelectrons)
					nonThermalTimescales(model.ntElectron, model, "electronCoolingTimes.dat");
				if (calculateNTprotons)
					nonThermalTimescales(model.ntProton, model, "protonCoolingTimes.dat");
			}
		
			if (calculateNTdistributions) {

				if (calculateNTelectrons) {
					injection(model.ntElectron, model);
					if (accMethod == 0)
						//distributionSpatialDiffusion(model.ntElectron, model);
						distributionMultiZone(model.ntElectron, model);
						//distributionSpatialDiffusionSteady(model.ntElectron, model);
					else
						distributionSpatialDiffusion(model.ntElectron, model);
				}
				else {
					readEandRParamSpace("electronInjection_vol", model.ntElectron.injection, 0, 1);
					readEandRParamSpace("electronDistribution_vol", model.ntElectron.distribution, 0, 1);
				}
				
				if (calculateNTprotons) {
					if (accMethod == 0)
						injection(model.ntProton,model);
					else
						injection(model.ntProton,model);
						//injectionFokkerPlanckOneZone(model.ntProton,model);
					if (accMethod == 0)
						distributionSpatialDiffusionSteady(model.ntProton, model);
						//distributionMultiZoneRadial(model.ntProton, model);
					else
						distributionSpatialDiffusionSteady(model.ntProton, model);
						//distributionFokkerPlanckCompleteSteadyState(model.ntProton, model);
						//distributionFokkerPlanckSpatialDiffusionTimeDependent(model.ntProton, model);
						//distributionFokkerPlanckRadial(model.ntProton,model);
				}
				else {
					readEandRParamSpace("protonInjection_vol", model.ntProton.injection,0,1);
					readEandRParamSpace("protonDistribution_vol", model.ntProton.distribution,0,1);
				}
			}
			else {
				if (calculateNTelectrons) {
					readEandRParamSpace("electronInjection_vol",model.ntElectron.injection,0,1);
					readEandRParamSpace("electronDistribution_vol",model.ntElectron.distribution,0,1);
				}
				if (calculateNTprotons) {
					readEandRParamSpace("protonInjection_vol",model.ntProton.injection,0,1);
					readEandRParamSpace("protonDistribution_vol",model.ntProton.distribution,0,1);
				}
			}
			
			if (calculateNonThermalLum)
				nonThermalRadiation(model, "lumNonThermal.dat");
			else {
				readEandRParamSpace("photonDensity", model.photon.distribution, 0, 0);
				readEandRParamSpace("NTphotonDensity", model.ntPhoton.distribution, 0, 0);
				readEandRParamSpace("NTphotonDensity", model.ntPhoton.injection, 0, 0);
				readEandRParamSpace("opticalDepth_gg", model.tau_gg, 0, 0);
			}
			
			if (calculateSecondaries)
				secondariesProcesses(model);
			else {
				readEandRParamSpace("muonDistribution", model.ntMuon.distribution, 0, 0);
				readEandRParamSpace("pionDistribution", model.ntChargedPion.distribution, 0, 0);
			}
			
			// NEUTRINO TRANSPORT
			if (calculateNeutrinos)
				injectionNeutrino(model.neutrino, model);
				
			// NEUTRON TRANPOSRT
			if (calculateNeutronInj) {
				injectionNeutrons(model);
				radiativeLossesNeutron(model.ntNeutron,model,"neutronTimescales.dat");
				if (calculateNeutronDis)
					distributionNeutronsAGN(model);
				if (calculateJetDecay)
					jetNeutronDecay(model);
			}
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
	}
	return 0;
}*/