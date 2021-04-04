#pragma once

#include "State.h"

//void distributionFast(Particle& p, State& st);
//void distributionDetailed(Particle& p, State& st);
void distributionOneZone_analytical(Particle& p, State& st);
void distributionMultiZone_analytical(Particle& p, State& st);
//void distributionNew(Particle& p, State& st);
//void distributionNewE(Particle& p, State& st);
//void distributionGAMERA(Particle& p, State& st);
//void distributionFokkerPlanck(Particle& p, State& st);
void distributionSecondaries(Particle& p, State& st);
void distributionMultiZone(Particle& p, State& st);
void distributionMultiZonePrueba(Particle& p, State& st);
void distributionFokkerPlanckOneZone(Particle& p, State& st);
void distributionFokkerPlanckMultiZone(Particle& p, State& st);
void distributionFokkerPlanckRadial(Particle& p, State& st);
void distributionFokkerPlanckComplete(Particle& p, State& st);
void distributionFokkerPlanckCompleteSteadyState(Particle& p, State& st);
void distributionSpatialDiffusion(Particle& p, State& st);
void distributionSpatialDiffusionSteady(Particle& p, State& st);
void distributionSpatialDiffusionSteady2(Particle& p, State& st);

//void distributionFokkerPlanckMultiZoneTimeDependent(Particle& p, State& st);
