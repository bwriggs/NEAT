#pragma once

//Function to calculate a fitness given
//an input and an output (higher fitness is better)

#include <vector>
#include <cstddef>

#include "BaseFunctor.h"
#include "Phenotype.h"

class FitnessRecorder {
public:
    virtual void operator()(Phenotype&, double) = 0;
};

class FitnessFunctor : public BaseFunctor {
public:
    FitnessFunctor(size_t, size_t);
    //get minimum fitness output to stop training
    double requiredFitness() const;
    virtual void operator()(std::vector<Phenotype>&, FitnessRecorder&) = 0;

    struct FitEnough {
        Phenotype phenotype;
        FitEnough(Phenotype);
    };
};
