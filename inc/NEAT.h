#pragma once

// main class
// used to evolve a maximizer
// to a given fitness function

#include "Parameters.h"
#include "Genome.h"
#include "Random.h"
#include "FitnessFunctor.h"
#include "TaskFunctor.h"
#include "LibNeat.h"

#include <memory>
#include <utility>
#include <cstddef>
#include <limits>
#include <vector>
#include <unordered_set>

class NEAT {
public:
    NEAT();
    NEAT(size_t);
    //create a network that does the task
    //evaluated by the fitness function
    Genome train(FitnessFunctor&, bool useBias = true);
    //get a functor implementing a genome
    std::unique_ptr<TaskFunctor> getFunctionFromGenome(const Genome&) const;

    //all [hyper-]parameters that control training
    Parameters params;
private:
    struct NEATDeps : public AllDeps {
    public:
        NEATDeps(Parameters*);
        Random *random() override;
        InnovationTracker *innovationTracker() override;
        Parameters *parameters() override;
        Random *_random;
        InnovationTracker *_innovationTracker;
        Parameters *_parameters;
    };
    NEATDeps allDeps;
    //in NEAT class (not in train method) so seeding happens "once"
    Random random;

    struct FitnessResult {
        size_t id;
        double fitness;
        FitnessResult(size_t = 0, double = std::numeric_limits<double>::min());
    };

    struct NEATFR : public FitnessRecorder  {
        NEATFR(std::unordered_map<size_t, double>&);
        void operator()(Phenotype&, double) override;
        std::unordered_map<size_t, double> &fitnesses;
    };

    struct Species {
        Species() = default;
        Species(Genome);

        std::unordered_set<size_t> ids;
        
        Genome representative;
        FitnessResult bestResult;
        FitnessResult bestCurrent;
    };

    struct FitnessComp {
        bool operator()(const FitnessResult&, const FitnessResult&) const;
    };

    struct CrossTicket {
        size_t genomeId;
        size_t speciesIdx;
        bool avail;
        CrossTicket() = default;
        CrossTicket(size_t, size_t, bool = true);
    };

    std::vector<size_t> getShares(std::vector<double> realShares, size_t total) const;
};

