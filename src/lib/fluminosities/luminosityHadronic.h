#pragma once

#include <fparticle/Particle.h>

/*[erg s^-1 cm^-3 ]*/
double luminosityHadronic(double E,const double density, double temp);

//double luminosityHadronic(double E, const Particle& creator,
//	const ParamSpaceValues& denf, const SpaceCoord& psc);