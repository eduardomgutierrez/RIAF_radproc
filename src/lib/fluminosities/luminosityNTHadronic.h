#pragma once
#include <fparticle/Particle.h>


/* Gamma ray emissivity for pp colissions, for a power-law htadron distribution
  From Kelner, Aharonian & Bugayov, 2006.*/

double luminosityNTHadronic(double E, const Particle& creator,
	const double density, const SpaceCoord& psc);
	
double luminosityThHadronic(double E, const double temp,
	const double density, const SpaceCoord& psc);