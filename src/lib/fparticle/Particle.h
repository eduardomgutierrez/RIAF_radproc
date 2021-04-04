#pragma once

// External libraries
#include <boost/property_tree/ptree.hpp>

// Local libraries
#include <fmath/physics.h>
#include <fparameters/ParamSpace.h>
#include <fparameters/ParamSpaceValues.h>

/*
	Particle State
	
	Among other things, it has the current distribution of particles of a given type within the model.
*/
class Particle {
public:

	const std::string id;
	double mass;
	double logEmin;
	double logEmax;

	const double emin() const { return EV_TO_ERG*pow(10.0, this->logEmin); };
	const double emax() const { return EV_TO_ERG*pow(10.0, this->logEmax); };

	template<typename T> T getpar(boost::property_tree::ptree& cfg,
									const std::string &path, const T& def = T{});

	Particle(const std::string& id);	
	Particle(Particle&& other);
	
	/* Creates the vectors for injection and distribution 
	   according to the registered dimensions. */
	void initialize();
	void configure(boost::property_tree::ptree& cfg);

	ParamSpace ps;
	ParamSpaceValues injection;
	ParamSpaceValues distribution;
	Dimension* eDim() const;
};

