#include "FitnessFunctor.h"

FitnessFunctor::FitnessFunctor(size_t in, size_t out) : BaseFunctor(in, out) {}

FitnessFunctor::FitEnough::FitEnough(Phenotype phenotype) : phenotype(phenotype) {}
