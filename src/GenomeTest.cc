#include "TestConfig.h"

#include "Genome.h"
#include "LibNeat.h"
#include "Random.h"
#include "Genome.h"

#include <iostream>
#include <algorithm>
#include <unordered_set>

struct CrossTest {
    size_t inputs;
    size_t outputs;
    bool useBias;

};

TEST_CASE("Add Connection") {
    TestDeps deps;
    Genome genome(2,2,0,&deps);
    std::list<NodeGene> ogNodes = genome.getNodeGenes();
    std::list<ConnectionGene> ogConns = genome.getConnectionGenes();
    genome.addConnection();
    CHECK(ogNodes == genome.getNodeGenes());
    std::list<ConnectionGene> newConns = genome.getConnectionGenes();
    CHECK((ogConns.size() + 1) == newConns.size());
    newConns.pop_back();
    CHECK(ogConns == newConns);
}

TEST_CASE("Add Node") {
    TestDeps deps;
    Genome genome(5,4,0,&deps);
    std::list<NodeGene> ogNodes = genome.getNodeGenes();
    std::list<ConnectionGene> ogConns = genome.getConnectionGenes();
    genome.addNode();
    std::list<NodeGene> newNodes = genome.getNodeGenes();
    std::list<ConnectionGene> newConns = genome.getConnectionGenes();
    size_t newNodeId = newNodes.back().id;
    newNodes.pop_back();
    CHECK(newNodes == ogNodes);
    auto ogit = ogConns.begin();
    auto newit = newConns.begin();
    for(; ogit != ogConns.end() && newit != newConns.end() && *ogit == *newit; ++ogit, ++newit);
    bool ok = (ogit != ogConns.end()) && (newit != newConns.end());
    CHECK(ok);
    CHECK(!newit->enabled);
    newit->enabled = true;
    CHECK(*ogit == *newit);
    ConnectionGene splitConn = *newit;
    for(; ogit != ogConns.end() && newit != newConns.end() && *ogit == *newit; ++ogit, ++newit);
    ok = ogit == ogConns.end() && newit != newConns.end();
    CHECK(ok);
    ok = newit->from == splitConn.from && newit->to == newNodeId && newit->weight == 1 && newit->enabled;
    CHECK(ok);
    ++newit;
    ok = newit != newConns.end();
    CHECK(ok);
    ok = newit->from == newNodeId && newit->to == splitConn.to && newit->weight == splitConn.weight && newit->enabled;
    CHECK(ok);
    ++newit;
    ok = newit == newConns.end();
    CHECK(ok);
}

TEST_CASE("Clone Test") {
    TestDeps deps;
    Genome base(2,1,0,&deps,true);
    Genome clone = base.cross(base, 0, 0, 1);
    bool same = clone.identicalStructure(base);
    if (!same) {
        std::cout << "Base:" << std::endl;
        std::cout << base << std::endl;
        
        std::cout << "Clone:" << std::endl;
        std::cout << clone << std::endl;
    }
    CHECK(same);
}


//generate aligned genomes specs
//make genomes conforming to them
//cross them
//assert the offspring inherited everything
//it was supposed to and nothing it wasn't
TEST_CASE("Cross Batch Test") {
    //read from a list (/generate) of aligned genome specs
    //i.e. the number of matching, disjoint and excess
    //generate inputs, outputs, bias, additional matching, disjoint, excess,
    //and fitnesses
    //create 2 genomes from these specs 
    //cross them and make sure the offspring has the correct genes from each
    //and nothing else 

    size_t N = 100;
    size_t minInputs = 1, maxInputs = 10;
    size_t minOutputs = 1, maxOutputs = 10;
    size_t minAddM = 0, maxAddM = 10;
    size_t minD = 1, maxD = 10;
    size_t minE = 1, maxE = 10;
    double minFit = 0; double maxFit = 10;

    TestDeps deps;
    for(size_t i = 0; i < N; ++i) {
        //generate parameters for this case
        size_t inputs = deps.random()->randomInt(minInputs, maxInputs + 1);
        size_t outputs = deps.random()->randomInt(minOutputs, maxOutputs + 1);
        bool useBias = deps.random()->randomInt(2) == 0;
        size_t addMatch = deps.random()->randomInt(minAddM, maxAddM + 1);
        size_t addDisjoint1 = deps.random()->randomInt(minD, maxD + 1);
        size_t addDisjoint2 = deps.random()->randomInt(minD, maxD + 1);
        size_t addExcess = deps.random()->randomInt(minE, maxE + 1);
        double fit1 = deps.random()->randomReal(minFit, maxFit);
        double fit2 = deps.random()->randomReal(minFit, maxFit);
        size_t exceptExcess = addMatch + addDisjoint1 + addDisjoint2;
        size_t totalAdd = exceptExcess + addExcess;

        //generate a genome with enough room for the additional genes
        std::list<NodeGene> nodeGenes;
        std::list<ConnectionGene> connGenes;
        size_t nid = 0, nextInnov = 0;
        for(; nid < inputs; ++nid) {
            nodeGenes.emplace_back(NodeGene::NodeType::Input, nid);
        }
        if(useBias) {
            nodeGenes.emplace_back(NodeGene::NodeType::Bias, nid++);
        }
        size_t outputBegin = nid;
        size_t outputEnd = nid + outputs;
        for(; nid < outputEnd; ++nid) {
            nodeGenes.emplace_back(NodeGene::NodeType::Output, nid);
        }
        for(size_t in = 0; in < (inputs + (useBias ? 1 : 0)); ++in) {
            for(size_t out = outputBegin; out < outputEnd; ++out) {
                connGenes.emplace_back(in, out, 0, true, nextInnov++);
            }
        }
        Genome g1(inputs, outputs, 0, &deps, nodeGenes, connGenes, useBias);
        //add in enough nodes to make all the extra genes
        size_t additionalAvailable = inputs * inputs + outputs * outputs; //conns within inputs OR outputs
        for (; additionalAvailable < totalAdd; ) {
            ConnectionGene *connGene = g1.findToSplit();
            if(connGene == nullptr) {
                std::cerr << "Cannot make genome" << std::endl;
                CHECK(false);
            }
            size_t newNodeId = g1.getNodeGenes().back().id + 1;
            g1.addNode(NodeGene(NodeGene::NodeType::Hidden, newNodeId));
            size_t newTo = g1.getNodeGenes().size() - ((connGene->from == connGene->to) ? 1 : 2);
            if(useBias && connGene->from != inputs && connGene->to != inputs) {
                --newTo;
            }
            connGene->enabled = false;
            g1.addConnection(ConnectionGene(connGene->from, newNodeId, 0, true, nextInnov++));
            g1.addConnection(ConnectionGene(newNodeId, connGene->to, 0, true, nextInnov++));

            additionalAvailable += newTo;
        }
        size_t initialMatching = g1.getConnectionGenes().size();

        //make a copy of this base genome
        Genome g2(g1);

        std::vector<ConnId> matchingGenes;
        std::vector<ConnId> disjointGenes1;
        std::vector<ConnId> disjointGenes2;
        std::vector<ConnId> excessGenes;
        size_t from, to;
        //find genes ahead of time on a dummy
        Genome dummy(g1);
        for(size_t m = 0; m < addMatch; ++m) {
            dummy.findNewConn(from, to);
            dummy.addConnection(ConnectionGene(from, to, 0));
            matchingGenes.emplace_back(from, to);
        }
        for(size_t d1 = 0; d1 < addDisjoint1; ++d1) {
            dummy.findNewConn(from, to);
            dummy.addConnection(ConnectionGene(from, to, 0));
            disjointGenes1.emplace_back(from, to);
        }
        for(size_t d2 = 0; d2 < addDisjoint2; ++d2) {
            dummy.findNewConn(from, to);
            dummy.addConnection(ConnectionGene(from, to, 0));
            disjointGenes2.emplace_back(from, to);
        }
        for(size_t e = 0; e < addExcess; ++e) {
            dummy.findNewConn(from, to);
            dummy.addConnection(ConnectionGene(from, to, 0));
            excessGenes.emplace_back(from, to);
        }

        //generate a random order to create matching and disjoint1/2 genes
        char types[exceptExcess];
        size_t ti = 0;
        for(size_t j = 0; j < addMatch; ++j) {
            types[ti++] = 'm';
        }
        for(size_t j = 0; j < addDisjoint1; ++j) {
            types[ti++] = '1';
        }
        for(size_t j = 0; j < addDisjoint2; ++j) {
            types[ti++] = '2';
        }
        std::shuffle(types, types + exceptExcess, deps.random()->getGen());

        size_t mIdx = 0, dIdx1 = 0, dIdx2 = 0;
        for(size_t j = 0; j < exceptExcess; ++j, ++nextInnov) {
            switch(types[j]) {
                case 'm':
                    g2.addConnection(*g1.addConnection(ConnectionGene(matchingGenes[mIdx].from, matchingGenes[mIdx].to, 0, true, nextInnov)));
                    ++mIdx;
                    break;
                case '1':
                    g1.addConnection(ConnectionGene(disjointGenes1[dIdx1].from, disjointGenes1[dIdx1].to, 0.1, true, nextInnov));
                    ++dIdx1;
                    break;
                case '2':
                    g2.addConnection(ConnectionGene(disjointGenes2[dIdx2].from, disjointGenes2[dIdx2].to, 0.2, true, nextInnov));
                    ++dIdx2;               
                    break;
            }
        }

        bool excess1 = deps.random()->evalProb(0.5);
        Genome &excessG = excess1 ? g1 : g2;
        for(size_t j = 0; j < addExcess; ++j, ++nextInnov) {
            excessG.addConnection(ConnectionGene(excessGenes[j].from, excessGenes[j].to, 0.3, true, nextInnov));
        }
        Genome offspring = g1.cross(g2, fit1, fit2, 0);
        bool take1 = fit1 > fit2 || (fit1 == fit2 && g1.getConnectionGenes().size() < g2.getConnectionGenes().size());
        size_t counts[4] = {0};
        for(const ConnectionGene &cg : offspring.getConnectionGenes()) {
            if(cg.weight == 0) {
                ++counts[0];
            }
            if(cg.weight == 0.1) {
                ++counts[1];
            }
            if(cg.weight == 0.2) {
                ++counts[2];
            }
            if(cg.weight == 0.3) {
                ++counts[3];
            }
        }

        bool ok = true;
        if(counts[0] != (initialMatching + addMatch)) ok = false;
        if(take1) {
            if (counts[1] != addDisjoint1 || counts[2] != 0 || counts[3] != (excess1 ? addExcess : 0)) {
                ok = false;
            }
        } else {
            if(counts[2] != addDisjoint2 || counts[1] != 0 || counts[3] != (excess1 ? 0 : addExcess)) {
                ok = false;
            }
        }
        
        if (!ok) {
            std::cout << g1 << 'X' << std::endl << g2 << "=" << std::endl << offspring << std::endl;
        }

        CHECK(ok);
    }
}
