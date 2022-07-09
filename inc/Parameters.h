#pragma once

#include <cstdlib>

struct Parameters {
    //c in sigmoid(x)=1/(1+e^(cx))
    //typically tuned to be near linear on [-0.5,0.5]
    //and achieve almost all of [0,1] for "small weights"
    double sigmoidConst;
    static constexpr double defaultSigmoidConst = -4.9;

    //similarity
    //coefficients weighted sum
    double excessWeight;
    double disjointWeight;
    double weightWeight;
    //maximum difference within species
    double deltaThresh;

    //mutation probabilities
    double probWeight;
    double probUW;
    double probNode;
    double probConn;

    //cross over
    //rate at which children are made from
    //parents of different species
    double interSpeciesRate;
    //rate at which inidividuals are simply
    //mutated into the next gen
    double noCrossRate;
    //the chance a gene is disabled
    //is the child if it is disabled in either
    //parent
    double chanceDisabled;
    //minimum size of a species for
    //its top performer to be cloned
    //into the next generation
    size_t minSizeForChamp;
    //maximum number of generations a species
    //can stagnate without being terminated
    size_t maxStagnantSpeciesGens;
    //maximum number of generations the population
    //can stagnate before being pruned
    size_t maxStagnantGens;
    //number of species left after a pruning
    size_t sizeAfterPrune;
    double survivalThresh;

    //other
    //number of networks in each generation
    size_t populationSize;
    //maximum number of generations to test before
    //giving up
    int maxGenerations;
    size_t minGenesToNormalize;

    //my own
    //minimum random weight
    double minRW;
    //minimum random weight
    double maxRW;
    //size of uniform weight perturbations
    double deltaUW;
};

//function to get default parameters as
//specified in the original NEAT paper
Parameters defaultParams();

double sigmoid(double);