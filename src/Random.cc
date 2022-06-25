#include "Random.h"
#include "Parameters.h"

Random::Random(RandomDeps *deps, size_t seed) : deps(deps), gen(seed) {}

Random::Random(RandomDeps *deps) : Random(deps, std::random_device{}()) {}

double Random::randomWeight() {
    return randRange(deps->parameters()->minRW, deps->parameters()->maxRW);
}

double Random::perturbedWeight(double weight) {
    return randRange(weight - deps->parameters()->deltaUW, weight + deps->parameters()->deltaUW);
}

bool Random::changeWeight() {
    return evalProb(deps->parameters()->probWeight);
}

bool Random::addConn() {
    return evalProb(deps->parameters()->probConn);
}

bool Random::addNode() {
    return evalProb(deps->parameters()->probNode);
}

bool Random::useRandomWeight() {
    return !evalProb(deps->parameters()->probUW);
}

bool Random::usePerturbedWeight() {
    return evalProb(deps->parameters()->probUW);
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
