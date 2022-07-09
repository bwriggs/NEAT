#include "NEAT.h"
#include "FitnessFunctor.h"
#include "InnovationTracker.h"
#include "Random.h"
#include "LibNeat.h"

#include <unordered_map>
#include <queue>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>

NEAT::NEAT() : params(defaultParams()), allDeps(&params), random(&allDeps) {
    allDeps._random = &random;
}
NEAT::NEAT(size_t seed) : params(defaultParams()), allDeps(&params), random(&allDeps, seed) {
    allDeps._random = &random;
}

Genome NEAT::train(FitnessFunctor &ff, bool useBias) {
    InnovationTracker innovTracker;
    allDeps._innovationTracker = &innovTracker;

    //mutation counts for debugging
    // size_t mCounts[5] = {0};

    //make initial population
    Genome base(ff.inputs(), ff.outputs(), 0, &allDeps, useBias);
    std::vector<Genome> population;
    for(size_t i = 0; i < params.populationSize; ++i) {
        population.push_back(base.keepStructure());
        population.back().setId(i+1);
    }
    std::list<Species> speciesList;
    std::unordered_map<size_t, GenomeInfo> id2Info;
    //std::unordered_map<size_t, std::list<Species>::iterator> genome2Species;
    FitnessResult bestFitness;
    //train
    for (size_t g = 0; params.maxGenerations < 0 || g < (size_t)params.maxGenerations; ++g) {
        //speciate
        //clear species
        id2Info.clear();
        for(Species &species : speciesList){
            species.ids.clear();
            species.bestCurrent = FitnessResult();
        }
        for(auto git = population.begin(); git != population.end(); ++git) {
            for(auto sit = speciesList.begin(); sit != speciesList.end(); ++sit) {
                if (git->delta(sit->representative) <= params.deltaThresh) {
                    id2Info[git->id()] = GenomeInfo(&*git, sit);
                    sit->ids.insert(git->id());
                    //deal with break / fallthrough ambiguity
                    goto FOUND_SPECIES;
                }
            }
            //make new species
            id2Info[git->id()] = GenomeInfo(&*git,
                speciesList.insert(speciesList.end(), Species(*git)));
            FOUND_SPECIES:
            while(false);
        }
        // std::cout << "species size: " << speciesList.size() << std::endl;

        // for(Species &species : speciesList) {
        //     std::cout << species.representative;
        //     for(size_t id : species.ids) {
        //         std::cout << ' ' << id;
        //     }
        //     std::cout << std::endl;
        // }

        //make phenotypes
        std::vector<Phenotype> phenotypes;
        for(Genome &g : population) {
            phenotypes.emplace_back(g);
        }

        //calculate fitness
        //std::unordered_map<size_t, double> fitnesses;
        NEATFR fr(id2Info);
        try {
            ff(phenotypes, fr);
        } catch (FitnessFunctor::FitEnough &fffe) {
            //return if fit enough
            return *id2Info[fffe.phenotype.id()].genome;
        }

        //update progress stats ? (should this be after species adjustment)
        for (auto &kvp : id2Info) {
            GenomeInfo &gi = kvp.second;
            if (gi.fitness > gi.species->bestResult.fitness) {
                gi.species->bestResult.fitness = gi.fitness;
                gi.species->bestResult.id = g;
            }
            if (gi.fitness > gi.species->bestCurrent.fitness) {
                gi.species->bestCurrent.fitness = gi.fitness;
                gi.species->bestCurrent.id = gi.genome->id();
            }
        }

        //remove stagnant species and update bestFitness
        for (auto sit = speciesList.begin(); sit != speciesList.end();) {
            if (sit->bestResult.fitness > bestFitness.fitness) {
                bestFitness.fitness = sit->bestResult.fitness;
                bestFitness.id = g;
                std::cout << "g=" << g << ", " << bestFitness.fitness << std::endl;
            }
            if ((g - sit->bestResult.id) > params.maxStagnantSpeciesGens) {
                // std::cout << "removing species" << std::endl;
                sit = speciesList.erase(sit);
            } else {
                ++sit;
            }
        }

        //prune if needed
        if ((g - bestFitness.id) > params.maxStagnantGens && speciesList.size() > params.sizeAfterPrune) {
            // std::cout << "pruning" << std::endl;
            //reproduce from only the top params.sizeAfterPrune species
            auto sit = speciesList.begin();
            for(size_t i = 0; i < params.sizeAfterPrune; ++i, ++sit);
            speciesList.erase(sit, speciesList.end());
            bestFitness = FitnessResult();
        }

        //need to eliminate the bottom X of the population first
        //first adjust by iterating over remaining species
        //keep track of lowest fitnesses and their ids
        std::priority_queue<FitnessResult, std::vector<FitnessResult>, FitnessComp> fitMinHeap;
        size_t totalRem = 0;
        for (Species &species : speciesList) {
            for(size_t id : species.ids) {
                //explicit fitness sharing
                id2Info[id].fitness /= species.ids.size();
                fitMinHeap.emplace(id, id2Info[id].fitness);
            }
            totalRem += species.ids.size();
        }
        // std::cout << "survivalFitness = " << survivalFitness << std::endl;
        //actually remove 
        size_t toRemove = (1 - params.survivalThresh) * fitMinHeap.size();
        for (size_t i = 0; i < toRemove; ++i) {
            if(fitMinHeap.empty()) {
                std::cout << "no more to rm" << std::endl;
                break;
            }
            FitnessResult toRemove = fitMinHeap.top(); fitMinHeap.pop();
            if(id2Info.find(toRemove.id) == id2Info.end()) {
                std::cout << "no such genome (" << toRemove.id << ')' << std::endl;
                continue;
            } 
            Species &species = *id2Info[toRemove.id].species;
            auto idt = species.ids.find(toRemove.id);
            if (idt != species.ids.end()) {
                species.ids.erase(idt);
            }
        }

        //maybe assign a new random representative now
        //also record number of champions
        size_t numChamps = 0;
        for (auto sit = speciesList.begin(); sit != speciesList.end();) {
            if (sit->ids.size() > params.minSizeForChamp) {
                ++numChamps;
            }
            if(sit->ids.size() > 0) {
                sit->representative = *id2Info[*sit->ids.begin()].genome;
                ++sit;
            } else {
                sit = speciesList.erase(sit);
            }
        }
                
        if(speciesList.size() == 0) {
            //default pop
            bestFitness = FitnessResult();
            for(size_t i = 0; i < params.populationSize; ++i) {
                population[i] = base.keepStructure();
                population[i].setId(i + 1);
            }
            continue;
        }
        
        //For normal, intraspecies crossover:
        //*use total(adjusted)/average(unadjusted) fitnesses of species
        //  to deterministically assign a number of offspring to each species
        //*use fitnesses within each species to calculate number of offSpring each genome will have
        //*produce a randomally ordered list of 2 * (total normal offspring) genome ids
        //  (2 for each offspring each genome will have)
        //*for each genome (particular order?), for each offspring to have,
        //  cross with the first ticket with a non-matching genome id

        //For interspecies crossover
        //just pick from a random other species tickets
        std::vector<Genome> newPopulation;

        std::vector<double> totalFitnesses;
        for (Species &species : speciesList) {
            double totalFitness = 0;
            for (size_t id : species.ids) {
                totalFitness += id2Info[id].fitness;
            }
            totalFitnesses.push_back(totalFitness);
        }

        //No cross, but mutate
        size_t numNoCross = params.noCrossRate * (params.populationSize - numChamps);
        size_t numCross = params.populationSize - numNoCross - numChamps;
        // std::cout << numChamps << " + " << numNoCross << " + " << numCross << " = " << (numChamps + numNoCross + numCross) << std::endl;
        std::vector<size_t> noCrossShares = getShares(totalFitnesses, numNoCross);
        std::vector<size_t> crossShares = getShares(totalFitnesses, numCross);

        //add no cross genomes and set up tickets for cross genomes
        auto ncit = noCrossShares.begin();
        auto cit = crossShares.begin();
        std::vector<CrossTicket> tickets;
        for (auto sit = speciesList.begin(); sit != speciesList.end(); ++sit, ++ncit, ++cit) {
            if(sit->ids.size() == 0) {
                assertrc(*ncit == 0 && *cit == 0);
                continue;
            }
            std::vector<size_t> ids(sit->ids.begin(), sit->ids.end());
            std::vector<double> specFits(sit->ids.size());
            for(auto it = ids.begin(); it != ids.end(); ++it) {
                specFits[it - ids.begin()] = id2Info[*it].fitness;
            }
            std::vector<size_t> specNoCrossShares = getShares(specFits, *ncit);
            std::vector<size_t> specCrossShares = getShares(specFits, *cit);
            assert(ids.size() == specNoCrossShares.size() && ids.size() == specCrossShares.size());
            //add noCross genomes
            for(size_t i = 0; i < specNoCrossShares.size(); ++i) {
                // std::cout << "no cross X " << specNoCrossShares[i] << " of " << *ncit << std::endl;
                newPopulation.insert(newPopulation.end(), specNoCrossShares[i], *id2Info[ids[i]].genome);
            }
            //add tickets for cross genomes
            for(size_t i = 0; i < specCrossShares.size(); ++i) {
                //make tickets (genome id, species id(x), taken bit) in a global vector
                //shuffle vector
                //add pointers to species specific lists
                //start crossing within species, marking taken bit and deleting pointer from list
                //when we randomally cross interspecies, search global pointer for first 
                tickets.insert(tickets.cend(), 2 * specCrossShares[i], CrossTicket(ids[i], cit - crossShares.begin()));
            }
        }

        std::shuffle(tickets.begin(), tickets.end(), random.getGen());
        std::vector<std::vector<CrossTicket*>> specTickets(speciesList.size());
        for(auto it = tickets.begin(); it != tickets.end(); ++it) {
            specTickets[it->speciesIdx].push_back(&*it);
        }
        for(size_t i = 0; i < speciesList.size(); ++i) {
            auto it = specTickets[i].begin();
            for(;;) {
                //advance to next available ticket
                for(; it != specTickets[i].end() && !(*it)->avail; ++it);
                if (it == specTickets[i].end()) {
                    break;
                }
                //take ticket
                Genome &g1 = *id2Info[(*it)->genomeId].genome;
                Genome *g2 = nullptr;
                (*it)->avail = false;
                auto it2(it);
                ++it2;
                size_t g2Id = 0;
                auto lastAvail(it2);

                //cross interspecies if can and should
                auto isit = tickets.begin();
                if(random.crossInterSpecies()) {
                    for(; isit != tickets.end() && (isit->speciesIdx == i || !isit->avail); ++isit);
                    if(isit != tickets.end()) {
                        // std::cout << "interspecies" << std::endl;
                        goto TakeInterSpecies;
                    }
                }
                
                //cross intraspecies if can and should
                for(; it2 != specTickets[i].end(); ++it2) {
                    if((*it2)->avail) {
                        lastAvail = it2;
                        if((*it2)->genomeId != g1.id()) {
                            break;
                        }
                    }
                }
                if (it2 == specTickets[i].end()) {
                    //not picky
                    it2 = lastAvail;
                    if(it2 == specTickets[i].end() || !(*it2)->avail) {
                        //not picky at all
                        for(isit = tickets.begin(); isit != tickets.end() && !isit->avail; ++isit);
                        if (isit == tickets.end()) {
                            it2 = it;
                            // std::cout << "cloning" << std::endl;
                        } else {
                            // std::cout << "not picky at all" << std::endl;
                            goto TakeInterSpecies;
                        }
                    } else {
                        // std::cout << "not picky" << std::endl;
                    }
                } else {
                    // std::cout << "normal" << std::endl;
                }
                //take ticket
                g2Id = (*it2)->genomeId;
                (*it2)->avail = false;
                goto AddGenome;

TakeInterSpecies:
                //take ticket
                g2Id = isit->genomeId;
                isit->avail = false;
AddGenome:
                g2 = id2Info[g2Id].genome;
                newPopulation.push_back(g1.cross(*g2, id2Info[g1.id()].fitness, id2Info[g2->id()].fitness, 0));
            }
        }

        //forget about last generations innovations
        innovTracker.clear();
        //mutate
        for(Genome &genome : newPopulation) {
            // Mutation mutation = random.getMutation();
            // // ++mCounts[static_cast<size_t>(mutation)];
            // switch(mutation) {
            //     case Mutation::PerturbWeight:
            //         // std::cout << 'p' << std::endl;
            //         genome.perturbWeight();
            //         break;
            //     case Mutation::NewWeight:
            //         // std::cout << 'r' << std::endl;
            //         genome.randomWeight();
            //         break;
            //     case Mutation::AddConnection:
            //         // std::cout << 'c' << std::endl;
            //         genome.addConnection();
            //         break;
            //     case Mutation::AddNode:
            //         // std::cout << 'n' << std::endl;
            //         genome.addNode();
            //         break;
            //     case Mutation::None:
            //         // std::cout << 'x' << std::endl;
            //         break;
            //     default:
            //         throw Exception("No such mutation");
            // }
            if (random.evalProb(params.probWeight)) {
                if(random.evalProb(params.probUW)) {
                    genome.perturbWeight();
                } else {
                    genome.randomWeight();
                }
            }
            if(random.evalProb(params.probConn)) {
                genome.addConnection();
            }
            if(random.evalProb(params.probNode)) {
                genome.addNode();
            }
        }

        //add champions
        for (Species &species : speciesList) {
            if (species.ids.size() > params.minSizeForChamp) {
                newPopulation.push_back(*id2Info[species.bestCurrent.id].genome);
            }
        }

        //assign new ids
        for(size_t i = 0; i < newPopulation.size(); ++i) {
            newPopulation[i].setId(i+1);
        }
        // std::cout << "new pop size " << newPopulation.size() << std::endl; 
        population = newPopulation;
    }

    // std::cout << "Counts (perturb, new, conn, node, none): " << std::endl;
    // size_t mSum = 0;
    // for(size_t i = 0; i < 5; ++i) {
    //     mSum += mCounts[i];
    // }
    // for(size_t i = 0; i < 5; ++i) {
    //     std::cout << mCounts[i] << " (" << (mCounts[i] / double(mSum)) << ')' << std::endl;
    // }

    // too many generations:
    throw Exception("Too many generations");
}

std::vector<size_t> NEAT::getShares(std::vector<double> realShares, size_t total) const {
    assert(realShares.size() > 0);
    double sum = 0;
    for (double d : realShares) {
        sum += d;
    }
    std::vector<size_t> intShares;
    size_t rem = total;
    auto rit = realShares.begin();
    for(; rit != realShares.end(); ++rit) {
        size_t take = total * (*rit / sum);
        if(take > rem) {
            rem = 0;
        } else {
            rem -= take;
        }
        intShares.insert(intShares.end(), take);
    }
    auto iit = intShares.begin();
    rit = realShares.begin();
    bool looped = false;
    for(size_t i = 0; i < rem; ++iit, ++rit) {
        if (iit == intShares.end()) {
            assert(rit == realShares.end());
            if (looped) {
                throw Exception("cannot partition");
            }
            looped = true;
            iit = intShares.begin();
            rit = realShares.begin();
        }
        if (abs(*rit) != 0.0) {
            ++(*iit);
            ++i;
        }
    }

    return intShares;
}

std::unique_ptr<TaskFunctor> NEAT::getFunctionFromGenome(const Genome &genome) const {
    return std::unique_ptr<TaskFunctor>(new Phenotype(genome));
}

NEAT::NEATFR::NEATFR(std::unordered_map<size_t, GenomeInfo> &fitnesses) : fitnesses(fitnesses) {}

void NEAT::NEATFR::operator()(Phenotype &p, double fitness) {
    fitnesses[p.id()].fitness = fitness;
}

NEAT::FitnessResult::FitnessResult(size_t id, double fitness) : id(id), fitness(fitness) {}

NEAT::Species::Species(Genome genome) : ids({genome.id()}), representative(genome) {}

bool NEAT::FitnessComp::operator()(const FitnessResult &lhs, const FitnessResult &rhs) const {
    return lhs.fitness > rhs.fitness;
}

NEAT::CrossTicket::CrossTicket(size_t genomeId, size_t speciesIdx, bool avail) 
    : genomeId(genomeId), speciesIdx(speciesIdx), avail(avail) {}

NEAT::NEATDeps::NEATDeps(Parameters *params) : _parameters(params) {}

Parameters *NEAT::NEATDeps::parameters() { return _parameters; }

Random *NEAT::NEATDeps::random() { return _random; }

InnovationTracker *NEAT::NEATDeps::innovationTracker() { return _innovationTracker; }

NEAT::GenomeInfo::GenomeInfo(Genome *genome, std::list<Species>::iterator species): genome(genome), species(species) {}
