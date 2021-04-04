#include "luminosityHadronic.h"
#include <gsl/gsl_sf_bessel.h>
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <nrMath/integrators.h>
#include <fmath/interpolation.h>
#include <flosses/crossSectionInel.h>
#include <fmath/physics.h>
#include <algorithm>

#include <gsl/gsl_math.h>

#include <boost/math/special_functions/bessel.hpp>



double fntHadron(double x, const Particle& p, const double density, const SpaceCoord& psc) //funcion a integrar   x=Ecreator; L=L(Ega)
{	
	double Kpi = 0.17;
	double eval = p.mass*cLight2+x/Kpi;
	
	//double Ekin = Ep/Kpi;
	
	double distCreator = 0.0;
	if (eval < p.emax() && eval > p.emin()) {
		distCreator = p.distribution.interpolate({ { 0, eval } }, &psc);
	}
	
	//double thr = 0.0016; //1GeV
	//double sigma = 30e-27*(0.95+0.06*log(Ekin/thr));
	
	double Eth = 1.22e9 * EV_TO_ERG;
	double l = log10((protonMass*cLight2+x/Kpi)/1.6); //evaluada en eval
	double sigma = 1.e-27 * (34.3+1.88*l+0.25*l*l) * P2(1.0 - pow(Eth/eval,4));
	double pionEmiss = cLight*density*sigma*distCreator/Kpi;  //sigma = crossSectionHadronicDelta(Ekin)
															  //lo saco asi pongo la condicion Ekin > Ethr en el limite de la int
	
	double result = pionEmiss/sqrt(P2(x)-P2(chargedPionMass*cLight2));
	return (x > chargedPionMass*cLight2 ? result : 0.0);
}

double luminosityNTHadronic(double E, const Particle& creator,
	const double density, const SpaceCoord& psc)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	double Max  = creator.emax();   //esto es un infinito 
	double Min  = std::max(E+P2(chargedPionMass*cLight2)/(4*E),thr*Kpi); //== Ekin > Ethr
	Min = E+0.25*P2(chargedPionMass*cLight2)/E;
	double integral = integSimpsonLog(Min, Max,[&](double x)
				{
					return fntHadron(x,creator,density,psc);
				},100);
	//double integral = RungeKuttaSimple(Min, Max, 
	//	[&](double x) {return fntHadron(x, creator, density, psc); }
	//);    //integra entre Emin y Emax

	double luminosity = 2.0*integral*P2(E); // [erg s^-1 cm^-3 ]

	return luminosity; 
}

double fntHadronTh(double x, const double temp, const double density, const SpaceCoord& psc) //funcion a integrar   x=Ecreator; L=L(Ega)
{	
	double Kpi = 0.17;
	double eval = protonMass*cLight2+x/Kpi;
	
	//double Ekin = Ep/Kpi;
	double g = eval / (protonMass*cLight2);
	double beta = sqrt(1.0-1.0/(g*g));
	double theta = boltzmann*temp/(protonMass*cLight2);
	double bessel = gsl_sf_bessel_Kn(2, 1.0/theta);
	double distCreator = (bessel > 0.0 ? density * g*g*beta / (theta*bessel) * exp(-g/theta) / (protonMass*cLight2) : 0.0);
	
	//double thr = 0.0016; //1GeV
	//double sigma = 30e-27*(0.95+0.06*log(Ekin/thr));
	
	double Eth = 1.22e9 * EV_TO_ERG;
	double l = log10((protonMass*cLight2+x/Kpi)/1.6); //evaluada en eval
	double sigma = 1.e-27 * (34.3+1.88*l+0.25*l*l) * P2(1.0 - pow(Eth/eval,4));
	double pionEmiss = cLight*density*sigma*distCreator/Kpi;  //sigma = crossSectionHadronicDelta(Ekin)
															  //lo saco asi pongo la condicion Ekin > Ethr en el limite de la int
	
	double result = pionEmiss/sqrt(P2(x)-P2(chargedPionMass*cLight2));
	return (x > chargedPionMass*cLight2 ? result : 0.0);
}

double luminosityThHadronic(double E, const double temp,
	const double density, const SpaceCoord& psc)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	double Max  = pow(10,1.5)*protonMass*cLight2;   //esto es un infinito 
	double Min  = std::max(E+P2(chargedPionMass*cLight2)/(4*E),thr*Kpi); //== Ekin > Ethr
	Min = E+0.25*P2(chargedPionMass*cLight2)/E;
	double integral = (Min < Max ? integSimpsonLog(Min, Max,[&](double x)
				{
					return fntHadronTh(x,temp, density,psc);
				},100) : 0.0);
	//double integral = RungeKuttaSimple(Min, Max, 
	//	[&](double x) {return fntHadron(x, creator, density, psc); }
	//);    //integra entre Emin y Emax

	double luminosity = 2.0*integral*P2(E); // [erg s^-1 cm^-3 ]

	return 2.0*luminosity; 
}





/*
class luminosityHadronic2;

double heaviside(double x,double a,double b)
{
	return (a <= x && x <= b ? 1.0 : 0.0);
}

double auxf3(double dGeV, void *params)
{
	struct four_d_params *p = (struct four_d_params *) params;
	double sGeV = p->p1;
	double gx = p->p2;
	double GammaGeV = p->p3;
	double isoGeV = p->p4;
	
	double pGeV = protonMass*cLight2/1.6e-3;
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double gd = (sGeV + dGeV*dGeV - pGeV*pGeV)/(2.0*dGeV*sqrt(sGeV));
	double betad = sqrt(1.0-1.0/(gd*gd));
	double epi = (dGeV*dGeV+piGeV*piGeV-pGeV*pGeV)/(2.0*dGeV);
	double ppi = sqrt(epi*epi-piGeV*piGeV);
	double aux1 = 0.5/(betad*gd*ppi);   // REVISAR ESTE p_pi
	double aux2 = P2(dGeV-isoGeV)+GammaGeV*GammaGeV;
	double h = heaviside(gx*piGeV,gd*(epi-betad*ppi),gd*(epi+betad*ppi));
	
	return aux1*h/aux2;    // [GeV^-3]
	
}

double inclusiveSigma(double sGeV)
{
	double pGeV = protonMass*cLight2/1.6e-3;
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double eta=sqrt(P2(sGeV-P2(piGeV)-P2(2.0*pGeV))-4.0*P2(piGeV*2.0*pGeV))/(2.0*piGeV*sqrt(sGeV));
	double pthr = 0.78; //[GeV]
	double sigma=0.0;;

	double p = sqrt(P2(0.5*sGeV/pGeV-pGeV)-P2(pGeV));
	
	if (p >= pthr && p <= 0.96) {
		sigma = 0.032*eta*eta+0.04*pow(eta,6)+0.047*pow(eta,8);
	} else if (p > 0.96 && p <= 1.27) {
		sigma = 32.6*pow(p-0.8,3.21);
	} else if (p > 1.27 && p <= 8.0) {
		sigma = 5.4*pow(p-0.8,0.81);
	} else if (p > 8.0) {
		sigma = 32.0 * log(p) + 48.5 / sqrt(p) - 59.5;
	}
	return 1.0e-27 * sigma;     // [cm^2]
}

double dsigma(double gx, double gr, double sGeV)
{
	double error;
	int status;
	
	double GammaGeV = 0.0575;
	double pGeV = protonMass*cLight2/1.6e-3;
	double piGeV = neutralPionMass * cLight2 / 1.6e-3;
	double isoGeV = 1.236;
	double atan1 = atan((sqrt(sGeV)-pGeV-isoGeV)/GammaGeV);
	double atan2 = atan((pGeV+piGeV-isoGeV)/GammaGeV);
	double aux1 = GammaGeV/(atan1-atan2);
	
	double Min = pGeV+piGeV;
	double Max = sqrt(sGeV)-pGeV;
	
	struct four_d_params auxf3_params = {sGeV,gx,GammaGeV,isoGeV};
	gsl_function gsl_auxf3;
		gsl_auxf3.function = &auxf3;
		gsl_auxf3.params = &auxf3_params;
	
	double integ = integrator_qags(&gsl_auxf3,Min,Max,0,1.0e-2,100,&error,&status);
	
	return inclusiveSigma(sGeV)*aux1*integ*piGeV;  // [cm^2]
}

double auxf2(double gx, void *params)
{
	struct four_d_params *p = (struct four_d_params *) params;
	double gr = p->p1;
	double epi = p->p2;
	double normtemp = p->p3;
	double sGeV = p->p4;
	
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double g = epi/piGeV;
	double beta = sqrt(1.0-1.0/(g*g));
	double betax = sqrt(1.0-1.0/(gx*gx));
	double q = sqrt(2.0*(gr+1.0))/normtemp;
	
	double f1 = exp(-q*g*gx*(1.0-beta*betax))-exp(-q*g*gx*(1.0+beta*betax));
	
	return f1/(betax*gx) * dsigma(gx,gr,sGeV);
}

double auxf(double gr, void *params)
{
	double error;
	int status;
	
	struct two_d_params *p = (struct two_d_params *) params;
	double epi = p->p1;
	double normtemp = p->p2;
	
	double piGeV = neutralPionMass*cLight2 / 1.6e-3;
	double pGeV = protonMass*cLight2/1.6e-3;
	double sGeV = 2.0*P2(pGeV)*(gr+1.0);
	double ji = (sGeV-4.0*P2(pGeV)+P2(piGeV))/(2.0*sqrt(sGeV));
	double Max = ji/piGeV;
	
	struct four_d_params auxf2_params = {gr,epi,normtemp,sGeV};
	gsl_function gsl_auxf2;
		gsl_auxf2.function = &auxf2;
		gsl_auxf2.params = &auxf2_params;
	
	double integral = integrator_qags(&gsl_auxf2,1.0,Max,0,1.0e-2,100,&error,&status);
	double result = (gr*gr-1.0) / sqrt(2.0*(gr+1.0)) * integral;
}

double fPion(double epi, double density, double temp)
{
	double error;
	int status;
	
	double normtemp = boltzmann* temp / (protonMass*cLight2);
	struct two_d_params auxf_params = {epi,normtemp};
	gsl_function gsl_auxf;
		gsl_auxf.function = &auxf;
		gsl_auxf.params = &auxf_params;

	double k2 = boost::math::cyl_bessel_k(2, 1.0/normtemp);
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double constant = cLight*density*density/(4.0*piGeV*normtemp*k2*k2);

	return constant * integrator_qags(&gsl_auxf,1.0,1.0e3,0,1.0e-2,100,&error,&status);
}

double fHadron(double epi, void *params)
{
	struct two_d_params *p = (struct two_d_params *) params;
	double density = p->p1;
	double temp= p->p2;
	
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double qpi = fPion(epi,density,temp);
	return qpi / sqrt(P2(epi)-P2(piGeV));   // [cm-3 s-1 GeV-2]
}


double luminosityHadronic2(double E, const double density, double temp)
{
	double error;
	int status;
	
	struct two_d_params fHadron_params = {density,temp};
	gsl_function gsl_fHadron;
		gsl_fHadron.function = &fHadron;
		gsl_fHadron.params = &fHadron_params;
	
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	double Min  = cHadron(E);
	Min = Min / 1.6e-3;
	double Max = 1.0e3;  // [en GeV]
	
	double integral = 2.0 * integrator_qags(&gsl_fHadron,Min,Max,0,1.0e-2,100,&error,&status);
		
	double jpp = integral * E*planck * 0.25/pi; // [erg s^-1 Hz^-1 cm^-3]
	jpp = jpp / (1.6e-3);
	return jpp;
}*/