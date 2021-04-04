// External libraries
#include <boost/property_tree/ptree.hpp>

// Local libraries
#include <fparameters/Dimension.h>

// Local headers
#include "Particle.h"

template<typename T> T Particle::getpar(boost::property_tree::ptree& cfg,
				const std::string &path, const T& def )
				{
					return cfg.get<T>("particle." + id + "." + path,
										cfg.get<T>("particle.default." + path, def));
				}

template int Particle::getpar<int>(boost::property_tree::ptree& cfg, const std::string &path, const int& def );
template double Particle::getpar<double>(boost::property_tree::ptree& cfg, const std::string &path, const double& def );

namespace {
	void zeroToN(Vector& v) {
		for (int i = 0; i < (int)v.size(); ++i) {
			v[i] = i;
		}
	}
}

Particle::Particle(const std::string& id):
	id(id),
 	mass{ 0.0 },
	injection{ ps, false }, // these PSVs are not initialized immediately, only after this PS has been constructed
 	distribution{ ps, false }
{
}

Particle::Particle(Particle&& other):
	ps(other.ps),
	injection(this->ps, false),
	distribution(this->ps, false),
	id(other.id),
	mass(other.mass),
	logEmax(other.logEmax),
	logEmin(other.logEmin)
{
	injection.values = std::move(other.injection.values);
	distribution.values = std::move(other.distribution.values);
}

void Particle::configure(boost::property_tree::ptree& cfg) {
	mass = cfg.get<double>("mass", mass);
	logEmin = cfg.get<double>("dim.energy.min", logEmin);
	logEmax= cfg.get<double>("dim.energy.max", logEmax);
}

void Particle::initialize() {
	injection.initialize();
	distribution.initialize();
}

Dimension* Particle::eDim() const { return ps.dimensions[0]; }