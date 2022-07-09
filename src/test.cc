#include <iostream>
#include <vector>
#include <cmath>

#define __NEAT_UNIDIRECTIONAL__

#include "NEAT.h"
#include "SimpleFitnessFunctor.h"

class XORFitnessFunctor : public SimpleFitnessFunctor {
public:
    XORFitnessFunctor(double reqFitness) : SimpleFitnessFunctor(2, 1, reqFitness), gen(std::random_device()()), dist(0,3) {
        tests.emplace_back(0, 0);
        tests.emplace_back(0, 1);
        tests.emplace_back(1, 0);
        tests.emplace_back(1, 1);
    }
    double err(size_t lhs, size_t rhs, double out) {
        return abs(out - double(lhs^rhs));
    }
    double operator()(Phenotype &p) {
        size_t iters = 8;
        double fitness = iters;
        for(size_t i = 0; i < iters; ++i) {
            auto tc = tests[i%4]; //tests[dist(gen)];
            size_t lhs = tc.first;
            size_t rhs = tc.second;
            fitness -= err(lhs, rhs, p(std::vector<double>{(double)lhs, (double)rhs})[0]);
        }
        fitness = (fitness * 4) / iters;
        return fitness * (fitness < 0 ? -fitness : fitness);
    }
    std::vector<std::pair<size_t, size_t>> tests;
    std::mt19937 gen;
    std::uniform_int_distribution<size_t> dist;
};

int main(int argc, char *argv[]) {
    int seed = -1;
    double reqFitness = 15.995;

    try {
        switch(argc) {
            case 3:
                seed = std::stoi(argv[2]);
            [[ fallthrough ]];
            case 2:
                reqFitness = std::stod(argv[1]);
            [[ fallthrough ]];
            case 1:
                break;
            default:
                throw Exception("usage");
    }
    } catch(...) {
        std::cerr << "Usage: " << argv[0] << " [reqFitness = " << reqFitness << " [seed = random]]" << std::endl;
        exit(1);
    }


    NEAT neat = seed >= 0 ? NEAT(seed) : NEAT(); //1234
    neat.params.maxGenerations = -1;
    XORFitnessFunctor xorff(reqFitness);
    try {
        Genome g = neat.train(xorff);
        auto func = Phenotype(g);//neat.getFunctionFromGenome(g);
        std::cout << "Solution:" << std::endl;
        std::cout << g << std::endl;
        for(size_t i = 0; i < 10; ++i) {
            std::cout << (func)(std::vector<double>{0,0})[0] << std::endl;
            std::cout << (func)(std::vector<double>{0,1})[0] << std::endl;
            std::cout << (func)(std::vector<double>{1,0})[0] << std::endl;
            std::cout << (func)(std::vector<double>{1,1})[0] << std::endl;
        }
        std::cout << func << std::endl;
    } catch (Exception &ex) {
        std::cout << ex.what() << std::endl;
    }

    return 0;
}
