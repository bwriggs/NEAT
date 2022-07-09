#include "Parameters.h"
#include <cmath>

Parameters defaultParams() {
    Parameters params;
    
    //similarity
    params.excessWeight = 1.0;
    params.disjointWeight = 1.0;
    params.weightWeight = 0.4;
    params.deltaThresh = 3.0;

    //mutation
    params.probWeight = 0.8;
    params.probUW = 0.9;
    params.probNode = 0.03;
    params.probConn = 0.05;

    //cross over
    params.interSpeciesRate = 0.001;
    params.noCrossRate = 0.25;
    params.chanceDisabled = 0.75;
    params.minSizeForChamp = 5;
    params.maxStagnantSpeciesGens = 15;
    params.maxStagnantGens = 20;
    params.sizeAfterPrune = 2;
    params.survivalThresh = 0.20;

    //other
    params.populationSize = 150;
    params.maxGenerations = 300;
    params.minGenesToNormalize = 20;

    //my own
    params.minRW = -8;
    params.maxRW = 8;
    params.deltaUW = 2.5;

    return params;
}

//fast enough?
//small networks may be able to tolerate a
//significantly truncated taylor series in place
//of std::exp
double sigmoid(double x) {
    return 1 / (1 + std::exp(Parameters::defaultSigmoidConst*x));
}
