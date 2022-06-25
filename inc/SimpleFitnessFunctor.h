#include "FitnessFunctor.h"

class SimpleFitnessFunctor : public FitnessFunctor {
public:
    SimpleFitnessFunctor(size_t, size_t, double);
    void operator()(std::vector<Phenotype>&, FitnessRecorder&) override;
    virtual double operator()(Phenotype&) = 0;
    double requiredFitness() const;
private:
    double reqFitness;
};