#include "luminositySynchrotron.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_synchrotron.h>
#include "opticalDepthSSA.h"
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>


double fSyn(double x, double E, const Particle& creator, const ParamSpaceValues& magf, const SpaceCoord& psc)         //funcion a integrar   x=Ee; L=L(Ega)
{

	const double magneticField = magf.get(psc);

	double distCreator;
	if (x < creator.emin() || x > creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, x } }, &psc); 
	}
	double Erest = creator.mass*cLight2;
	double cte = sqrt(3.0)*P3(electronCharge)*magneticField / (planck*Erest);
	double Echar = 3.0*electronCharge*planck*magneticField*P2(x/Erest) / 
					(4.0*pi*creator.mass*cLight);
	
	double aux = E/Echar;  //aca el aux es el x real

	double result = cte*1.85*distCreator*pow(aux,1.0/3.0)*exp(-aux);  

	return result;
}


double luminositySynchrotron(double E, const Particle& c, const SpaceCoord& psc, const ParamSpaceValues& magf)
{
	
	double integralS = RungeKuttaSimple(c.emin(), c.emax(), [&](double e){
		return fSyn(e, E, c, magf, psc);
	});

	double luminosityS = integralS*E; //multiplico por E asi obtengo luminosidad

	if (luminosityS > 0.0){ return luminosityS; }
	else { return 0.0; }

}



double fSyn2(double x, double E, const Particle& creator, double magf, const SpaceCoord& psc)         //funcion a integrar   x=Ee; L=L(Ega)
{
	double distCreator;
	if (x < creator.emin() || x > creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, x } }, &psc); 
	}
	
	double Erest = creator.mass*cLight2;
	double cte = sqrt(3.0)*P3(electronCharge)*magf / (planck*Erest);
	double Echar = 3.0*electronCharge*planck*magf*P2(x/Erest) / 
					(4.0*pi*creator.mass*cLight);
	
	double aux = E/Echar;  //aca el aux es el x real

	double result = cte*1.85*distCreator*pow(aux,1.0/3.0)*exp(-aux);

	return result;
}

double luminositySynchrotron2(double E, const Particle& c, const SpaceCoord& psc, double magf)
{
	double integralS = RungeKuttaSimple(c.emin(), c.emax(), [&](double e){
		return fSyn2(e, E, c, magf, psc);
	});

	double luminosityS = integralS*E; //multiplico por E asi obtengo luminosidad

	if (luminosityS > 0.0){ return luminosityS; }
	else { return 0.0; }
	
}

double fSyn3(double Ee, double E, const Particle& creator, double magf, const SpaceCoord& psc)
{
	double distCreator = (Ee > creator.emin() && Ee < creator.emax()) ? 
				creator.distribution.interpolate({{0,Ee}},&psc) : 0.0;
	double gamma_e = Ee/(electronMass*cLight2);
	double muMin = -0.999;
	double muMax = 0.999;
	size_t nMu = 30;
	double dMu = (muMax-muMin)/nMu;
	double sum = 0.0;
	double mu = muMin;
	double Constant = 3.0*planck*electronCharge*magf*gamma_e*gamma_e / (4.0*pi*electronMass*cLight);
	for (size_t jMu=0;jMu<nMu;jMu++) {
		double Ec = Constant * sqrt(1.0-mu*mu);
		double x = E/Ec;
		double integ = (x < 10.0) ? integSimpson(log(x),log(10.0),[&](double logxp)
						{
							double xp = exp(logxp);
							return gsl_sf_bessel_Knu(5.0/3.0,xp)*xp;
						},50) : 0.0;
		sum += dMu*integ;
		mu += dMu;
	}
	return distCreator/(gamma_e*gamma_e) * sum;
}

double phiSy(double x, double si)
{
	return (x < 800) ? si*gsl_sf_synchrotron_1(x) : 0.0;
}


double luminositySynchrotronExact(double Eph, const Particle& c, const SpaceCoord& psc, double magf)
{
	double integralS = integSimpsonLog(c.emin(),c.emax(),[Eph,magf,&c,&psc](double Ee)
			{
				double g = Ee/(c.mass*cLight2);
				double integAng = integSimpson(0.0,0.9999,[&](double mu)
				//double integAng = qromo(0.0,0.9999,[&](double mu)
						{
							double si = sqrt(1.0-mu*mu);
							double Ec = 3.0*electronCharge*planck*magf*si*g*g/(4*pi*c.mass*cLight);
							double x = Eph/Ec;
							return phiSy(x,si);
						//});
						},20);
				return c.distribution.interpolate({{0,Ee}},&psc)*integAng;
			},80);
	double constant = sqrt(3.0)*gsl_pow_3(electronCharge)*magf/(planck*c.mass*cLight2);
	return (integralS > 0.0) ? constant*integralS*Eph : 0.0;
}

double luminositySynchrotronExactSec(double Eph, const Particle& c, const SpaceCoord& psc, double magf)
{
	double integralS = integSimpsonLog(c.emin(),c.emax(),[Eph,magf,&c,&psc](double Ee)
			{
				double g = Ee/(c.mass*cLight2);
				double tau = Eph/Ee;
				double integAng = integSimpson(0.0,0.9999,[Eph,g,magf,&c,tau](double mu)
						{
							double si = sqrt(1.0-mu*mu);
							double Ec = 3.0*electronCharge*planck*magf*si*g*g/(4*pi*c.mass*cLight);
							double x = Eph/Ec * (1.0/(1.0-tau));
							if (tau < 1.0)
								return phiSy(x,si)*(1.0-tau);
							else return 0.0;
						},20);
				return c.distribution.interpolate({{0,Ee}},&psc)*integAng;
			},80);
	double constant = sqrt(3.0)*gsl_pow_3(electronCharge)*magf/(planck*c.mass*cLight2);
	return (integralS > 0.0) ? constant*integralS*Eph : 0.0;
}

double Sfunc(double z)
{
	double a = integSimpsonLog(z,50,[&](double q)
				{
					return (pow(q,1.5) - pow(z,1.5)) * gsl_sf_bessel_Knu(5.0/3.0,q);
				},40);
	return (a > 0.0 ? 2.0/(3.0*sqrt(z)) * a / (6*pi) : 0.0);
}

double luminositySynchrotronBackReaction(double Eph, const Particle& c, const SpaceCoord& psc, double magf)
{
	double eps = sqrt(3.0*electronCharge/(thomson*magf));
	double integralS = integSimpsonLog(c.emin(),c.emax(),[Eph,magf,&c,&psc,eps](double Ee)
			{
				double tau = Eph/Ee;
				double g = Ee/(c.mass*cLight2);
				if (g > eps) {
					double Ec = 1.5*g*g*electronCharge*magf/(c.mass*cLight) * planck;
					double x = (Eph/Ec) / (1.0-tau);
					return (tau < 1.0) ? c.distribution.interpolate({{0,Ee}},&psc)*Sfunc(x)*(1.0-tau) : 0.0;
				} else {
					double integAng = integSimpson(0.0,0.9999,[&](double mu)
						{
							double si = sqrt(1.0-mu*mu);
							double Ec = 3.0*electronCharge*planck*magf*si*g*g/(4*pi*c.mass*cLight);
							double x = Eph/Ec;
							return phiSy(x,si);
						},20);
					return c.distribution.interpolate({{0,Ee}},&psc)*integAng;
				}
			},80);
	double constant = sqrt(3.0)*gsl_pow_3(electronCharge)*magf/(planck*c.mass*cLight2);
	return (integralS > 0.0) ? constant*integralS*Eph : 0.0;
}
/*
double luminositySynchrotron_conSSA(double E, const Particle& creator)
{

	double Emax = creator.emax();
	double Emin = creator.emin();

	double tau = opticalDepthSSA(E, creator.mass, Emin, Emax, creator);  //E=Eph

	double factorSSA = (1.0-exp(-tau))/tau;

	if (factorSSA > 1.0 || factorSSA == 0.0) //1.0e-3)  //lo cambio solo en el caso que interesa
	{ factorSSA = 1.0; }	
	
	//double integralS = RungeKuttaSimple(Emin, Emax, bind(fSyn,_1,E,creator));

	double integralS = RungeKuttaSimple(creator.emin(), creator.emax(),   //RungeKuttaSimple(double a, double b, fun1 f)
		[E, &creator](double x){
		return fSyn(x, E, creator);  //double fSyn(double x, double E, const Particle& creator)
	});


	double luminosityS = factorSSA*integralS*E; //multiplico por E asi obtengo luminosidad

	if (luminosityS > 0.0){return luminosityS ; }
	else {return 0.0;} 
	 
}*/


////////////////////////////////////////////////////////////////////////////////////////////////////

double fSynS(double Ee, double Eg, const Particle& creator, const SpaceCoord& psc, double magf)         //funcion a integrar   x=Ee; L=L(Ega)
{
	double distCreator = creator.distribution.interpolate({{0,Ee}},&psc);//interpol(x,Ecreator,Ncreator,Ncreator.size()-1);

	double cte	= pow(3.0,0.5)*P3(electronCharge)*magf/(planck*creator.mass*cLight2);

	double Echar = 3*electronCharge*planck*magf*P2(Ee)/(4*pi*P3(creator.mass)*cLight*P2(cLight2));
	
	double tau = Eg/Ee;
	double aux = Eg/(Echar*(1-tau));  //aca el aux es el x real

	double result = cte*1.85*(1-tau)*distCreator*pow(aux,(1.0/3.0))*exp(-aux);  

	return result;    // esta condicion la puse en el limite inferior tau<1 ? result : 0.0;
}



double luminositySynchrotronSec(double Eg, const Particle& c, const SpaceCoord& psc, double magf)
{
	double Emax = c.emax();
	double Emin = c.emin();
	double inf  = std::max(Emin,Eg);   //esto lo agrego asi le saco la condicion sobre tau < 1

	double integralS = RungeKuttaSimple(Emin, Emax, [&Eg,&c,&psc,magf](double Ee){
		return fSynS(Ee,Eg,c,psc,magf);
	});

	double luminosityS = integralS*Eg; //multiplico por E asi obtengo luminosidad
	                                 //divido por E asi obtengo emisividad y no luminosidad

	return luminosityS ; 

}