#include "Random.h"
#include "Parameters.h"
#include <iostream>

Random::Random(RandomDeps *deps, size_t seed) : deps(deps), gen(seed) {
    Parameters *params = deps->parameters();
    mutationDist = std::discrete_distribution<size_t>{params->probWeight * params->probUW, params->probWeight * (1 - params->probUW),
        params->probConn, params->probNode, 1 - params->probWeight - params->probConn - params->probNode};
}

Random::Random(RandomDeps *deps) : Random::Random(deps, std::random_device()()) {}

double Random::randomWeight() {
    return randRange(deps->parameters()->minRW, deps->parameters()->maxRW);
}

double Random::perturbedWeight(double weight) {
    double newWeight = randRange(weight - deps->parameters()->deltaUW, weight + deps->parameters()->deltaUW);
    if (newWeight < deps->parameters()->minRW) {
        return deps->parameters()->minRW;
    }
    if(newWeight > deps->parameters()->maxRW) {
        return deps->parameters()->maxRW;
    }
    return newWeight;
}

Mutation Random::getMutation() {
    size_t mutation = mutationDist(gen);
    return static_cast<Mutation>(mutation);
}

bool Random::crossInterSpecies() {
    return evalProb(deps->parameters()->interSpeciesRate);
}

size_t Random::randomInt(size_t end) {
    return intDist(gen) % end;
}

size_t Random::randomInt(size_t min, size_t end) {
    return min + randomInt(end - min);
}

double Random::randomReal(double min, double max) {
    return min + realDist(gen) * (max - min);
}

bool Random::evalProb(double prob) {
    return realDist(gen) < prob;
}

double Random::randRange(double min, double max) {
    return min + (max - min) * realDist(gen);
}

bool Random::stillDisabled() {
    return evalProb(deps->parameters()->chanceDisabled);
}

std::mt19937 &Random::getGen() {
    return gen;
}
