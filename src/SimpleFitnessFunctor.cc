#include "SimpleFitnessFunctor.h"
#include "Phenotype.h"

SimpleFitnessFunctor::SimpleFitnessFunctor(size_t in, size_t out, double reqFitness) : FitnessFunctor(in, out), reqFitness(reqFitness) {}

void SimpleFitnessFunctor::operator()(std::vector<Phenotype> &population, FitnessRecorder &fr) {
    for(Phenotype &p : population) {
        double fitness = operator()(p);
        fr(p, fitness);
        if(fitness >= reqFitness) {
            throw FitEnough(p);
        }
    }
}

double SimpleFitnessFunctor::requiredFitness() const {
    return reqFitness;
}
