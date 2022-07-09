#pragma once

// a class to encapsulate all random logic

#include "LibNeat.h"

#include <random>

class Parameters;

class Random {
public:
    //use specific seed (for testing)
    Random(RandomDeps*, size_t);
    //use seed from std::random_device
    Random(RandomDeps*);
    //get a random weight
    double randomWeight();
    //get a perturbed weight
    double perturbedWeight(double);
    //return true if the mutation should be done
    Mutation getMutation();
    bool crossInterSpecies();
    bool stillDisabled();
    std::mt19937 &getGen();

    //get a random int
    size_t randomInt(size_t, size_t);
    size_t randomInt(size_t);

    //get a random real
    double randomReal(double, double);
    
    //return true arg fraction of the time, false the rest
    bool evalProb(double);
private:
    RandomDeps *deps;

    //std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<size_t> intDist;
    std::uniform_real_distribution<double> realDist;
    std::discrete_distribution<size_t> mutationDist; 

    double randRange(double, double);
};
