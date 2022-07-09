#include "TestConfig.h"

#include "Genome.h"
#include "LibNeat.h"
#include "Random.h"
#include "Genome.h"
#include "Phenotype.h"

#include <iostream>

TEST_CASE("Basic Example") {
    TestDeps deps;
    Genome genome(4, 4, 0, &deps);
    genome.addConnection(ConnectionGene(4, 0, 1, true));
    std::cout << "Genome:" << std::endl;
    std::cout << genome << std::endl;

    Phenotype phenotype(genome);
    std::cout << "Phenotype:" << std::endl;
    std::cout << phenotype << std::endl;

    std::vector<double> input {0.1, 0.1, 0.1, 0.1};
    for(size_t i = 0; i < 5; ++i) {
        std::vector<double> output = phenotype(input);
        std::cout << output[0];
        for(size_t j = 1; j < output.size(); ++j) {
            std::cout << ", " << output[j];
        }
        std::cout << std::endl;
    }

    CHECK(true);
}