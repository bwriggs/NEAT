#include "Genome.h"
#include "Random.h"
#include "InnovationTracker.h"
#include "Parameters.h"

#include <iostream>
#include <unordered_set>
#include <set>
#include <utility>

//used for initial population
//fully connected network with no hidden nodes
Genome::Genome (size_t in, size_t out, size_t id, GenomeDeps *deps, bool useBias) 
: inputs(in), outputs(out), _id(id), deps(deps), useBias(useBias) {
    //make input nodes
    for (size_t i = 0; i < in; ++i) {
        addNode(NodeGene(NodeGene::NodeType::Input, i));
    }

    if (useBias) {
        addNode(NodeGene(NodeGene::NodeType::Bias, in));
    }

    size_t minOutput = nodeGenes.size();
    for(size_t i = 0; i < out; ++i) {
        addNode(NodeGene(NodeGene::NodeType::Output, nodeGenes.size()));
    }
    size_t endOutput = nodeGenes.size();
    for(size_t i = 0; i < (in + (useBias ? 1 : 0)); ++i) {
        for(size_t j = minOutput; j < endOutput; ++j) {
            connectionGenes.emplace_back(i, j, deps->random()->randomWeight());
            deps->innovationTracker()->assignInnov(connectionGenes.back());
        }
    }
}

Genome::Genome(size_t in, size_t out, size_t id, GenomeDeps *deps, std::list<NodeGene> nodeGenes, 
std::list<ConnectionGene> connectionGenes, bool useBias)
: inputs(in), outputs(out), _id(id), deps(deps), nodeGenes(nodeGenes), connectionGenes(connectionGenes), useBias(useBias) {}

Genome::Genome() : Genome(0,0,0,nullptr) {}

Genome Genome::keepStructure() {
    Genome newGenome = *this;

    for(ConnectionGene &cg : newGenome.connectionGenes) {
        cg.weight = deps->random()->randomWeight();
    }

    return newGenome;
}

void Genome::addNode(NodeGene ng) {
    nodeGenes.emplace_back(ng);
}

std::list<ConnectionGene>::iterator Genome::addConnection(ConnectionGene cg) {
    return connectionGenes.insert(connectionGenes.end(), cg);
}

bool Genome::identicalStructure(const Genome &genome) const {
    if(nodeGenes.size() != genome.nodeGenes.size() || connectionGenes.size() != genome.connectionGenes.size()) {
        return false;
    }
    auto nit = nodeGenes.begin();
    auto onit = genome.nodeGenes.begin();
    for(; nit != nodeGenes.end(); ++nit, ++onit) {
        if (*nit != *onit) {
            return false;
        }
    }

    auto cit = connectionGenes.begin();
    auto ocit = genome.connectionGenes.begin();
    for(; cit != connectionGenes.end(); ++cit, ++ocit) {
        //ignore weight and enable
        if(cit->innovationNumber != ocit->innovationNumber ||
            cit->from != ocit->from || cit->to != ocit->to) {
            return false;
        }
    }

    return true;
}

ConnectionGene *Genome::findToSplit() {
    std::vector<ConnectionGene*> enabled;
    for (auto it = connectionGenes.begin(); it != connectionGenes.end(); ++it) {
        if (it->enabled && !(useBias && it->from == inputs)) {
            enabled.push_back(&*it);
        }
    }
    if(enabled.size() == 0) {
        return nullptr;
    }
    return enabled[deps->random()->randomInt(enabled.size())];
}

bool Genome::findNewConn(size_t &from, size_t &to) {
    std::unordered_set<ConnId> taken;
    for(auto &cg : connectionGenes) {
        taken.emplace(cg.from, cg.to);
    }

    size_t randomFrom = deps->random()->randomInt(nodeGenes.size());
    for (size_t i = 0; i < nodeGenes.size(); ++i) {
        from = (randomFrom + i) % nodeGenes.size();
        size_t randomTo = deps->random()->randomInt(nodeGenes.size());
        for(size_t j = 0; j < nodeGenes.size(); ++j) {
            to = (randomTo + j) % nodeGenes.size();
            //disallow a connection from a bias to an input or bias
            if (useBias && (from == inputs) && (to <= inputs)) {
                continue;
            }
            //disallow a connection to a bias
            if (useBias && (to == inputs)) {
                continue;
            }
            if (taken.find(ConnId(from, to)) == taken.end()) {
                return true;
            }
        }
    }

    return false;
}

void Genome::addNode() {
    ConnectionGene *toSplit = findToSplit();
    if(toSplit == nullptr) {
        std::cerr << "none to split" << std::endl;
        return;
    }
    toSplit->enabled = false;
    size_t newNodeId = nodeGenes.size() == 0 ? 0 : (nodeGenes.back().id + 1);
    nodeGenes.emplace_back(NodeGene::NodeType::Hidden, newNodeId);
    connectionGenes.emplace_back(toSplit->from, newNodeId, 1);
    deps->innovationTracker()->assignInnov(connectionGenes.back());
    connectionGenes.emplace_back(newNodeId, toSplit->to, toSplit->weight);
    deps->innovationTracker()->assignInnov(connectionGenes.back());
}

void Genome::addConnection() {
    size_t from, to;
    if (findNewConn(from, to)) {
        connectionGenes.emplace_back(from, to, deps->random()->randomWeight(), true);
        deps->innovationTracker()->assignInnov(connectionGenes.back());
    } else {
        std::cerr << "fully connected" << std::endl;
        //throw 1;
    }
}

void Genome::randomWeight() {
    size_t take = deps->random()->randomInt(connectionGenes.size());
    auto it = connectionGenes.begin();
    for(size_t i = 0; i < take; ++i, ++it);
    it->weight = deps->random()->randomWeight();
}

void Genome::perturbWeight() {
    size_t take = deps->random()->randomInt(connectionGenes.size());
    auto it = connectionGenes.begin();
    for(size_t i = 0; i < take; ++i, ++it);
    it->weight = deps->random()->perturbedWeight(it->weight);
}

std::pair<Genome::AlignedGenome, Genome::AlignedGenome> Genome::alignGenomes(const Genome &other) const {
    const std::list<ConnectionGene> &cgs1 = connectionGenes;
    const std::list<ConnectionGene> &cgs2 = other.getConnectionGenes();
    size_t min1 = cgs1.front().innovationNumber, max1 = cgs1.back().innovationNumber;
    size_t min2 = cgs2.front().innovationNumber, max2 = cgs2.back().innovationNumber;

    auto it1 = cgs1.begin();
    auto it2 = cgs2.begin();

    AlignedGenome ag1;
    AlignedGenome ag2;

    for (;it1 != cgs1.end() && it2 != cgs2.end();) {
        if(it1->innovationNumber < it2->innovationNumber &&
        (it1->innovationNumber < min2 || it1->innovationNumber > max2)) {
            // excess
            ag1.excess.push_back(&*it1); ++it1;
        } else if(it1->innovationNumber < it2->innovationNumber) {
            // disjoint
            ag1.disjoint.push_back(&*it1); ++it1;
        } else if(it1->innovationNumber == it2->innovationNumber) {
            //match
            ag1.match.push_back(&*it1); ++it1;
            ag2.match.push_back(&*it2); ++it2;
        } else if(it2->innovationNumber < it1->innovationNumber &&
        (it2->innovationNumber < min1 || it2->innovationNumber > max1)) {
            // excess
            ag2.excess.push_back(&*it2); ++it2;
        } else if(it2->innovationNumber < it1->innovationNumber) {
            // disjoint
            ag2.disjoint.push_back(&*it2); ++it2;
        }
    }
    //just excess I think
    for(; it1 != cgs1.end(); ++it1) {
        ag1.excess.push_back(&*it1);
    }
    for(; it2 != cgs2.end(); ++it2) {
        ag2.excess.push_back(&*it2);
    }

    return std::make_pair(ag1, ag2);
}

double Genome::delta(const Genome &other) const {
    std::pair<AlignedGenome, AlignedGenome> alignedGenomes = alignGenomes(other);
    double N = 1;
    if(connectionGenes.size() >= other.getConnectionGenes().size()
    && connectionGenes.size() >= deps->parameters()->minGenesToNormalize) {
        N = connectionGenes.size();
    } else if (other.getConnectionGenes().size() >= connectionGenes.size()
    && other.getConnectionGenes().size() >= deps->parameters()->minGenesToNormalize) {
        N = other.getConnectionGenes().size();
    }
    double avgWeightDiff = 0;
    auto &matches1 = alignedGenomes.first.match;
    auto &matches2 = alignedGenomes.second.match;
    auto mit1 = matches1.begin();
    auto mit2 = matches2.begin();
    for(;mit1 != matches1.end() && mit2 != matches2.end(); ++mit1, ++mit2) {
        if((*mit1)->weight >= (*mit2)->weight) {
            avgWeightDiff += (*mit1)->weight - (*mit2)->weight;
        } else {
            avgWeightDiff += (*mit2)->weight - (*mit1)->weight;
        }
    }
    if (mit1 != matches1.end() || mit2 != matches2.end()) {
        throw Exception("mismatched matches");
    }
    if (matches1.size() > 1) {
        avgWeightDiff /= matches1.size();
    }

    return (deps->parameters()->excessWeight * (alignedGenomes.first.excess.size() + alignedGenomes.second.excess.size())
            + deps->parameters()->disjointWeight * (alignedGenomes.first.disjoint.size() + alignedGenomes.second.disjoint.size())) / N
            + deps->parameters()->weightWeight * avgWeightDiff;
}

Genome Genome::cross(const Genome &other, double f1, double f2, size_t newId) {
    // align genes
    std::pair<AlignedGenome, AlignedGenome> alignedGenomes = alignGenomes(other);

    //keep genes in order
    std::set<NodeGene, NodeCompare> newNGSet;
    std::set<ConnectionGene, ConnCompare> newConnSet;
    std::unordered_map<size_t, const NodeGene*> id2Node1;
    std::unordered_map<size_t, const NodeGene*> id2Node2;
    for (const NodeGene &ng : nodeGenes) {
        id2Node1[ng.id] = &ng;
    }
    for (const NodeGene &ng : other.getNodeGenes()) {
        id2Node2[ng.id] = &ng;
    }

    // take matching genes randomally
    auto &matches1 = alignedGenomes.first.match;
    auto &matches2 = alignedGenomes.second.match;
    auto mit1 = matches1.begin();
    auto mit2 = matches2.begin();
    for (; mit1 != matches1.end() && mit2 != matches2.end(); ++mit1, ++mit2) {
        bool from1 = deps->random()->evalProb(0.5);
        ConnectionGene newConnGene = from1 ? **mit1 : **mit2;
        auto &id2Node = from1 ? id2Node1 : id2Node2;
        if(!((*mit1)->enabled && (*mit2)->enabled) && deps->random()->stillDisabled()) {
            newConnGene.enabled = false;
        } else {
            newConnGene.enabled = true;
        }
        newConnSet.insert(newConnGene);
        newNGSet.emplace(*id2Node[newConnGene.from]);
        newNGSet.emplace(*id2Node[newConnGene.to]);
    }
    if (mit1 != matches1.end() || mit2 != matches2.end()) {
        throw Exception("mismatched matches");
    }

    //determine which parents will contribute their
    //disjoint and excess (take more fit, smaller if equal fitness)
    bool oneMoreFit = f1 > f2 || (f1 == f2 && connectionGenes.size() < other.getConnectionGenes().size());
    AlignedGenome &moreFitAG = oneMoreFit ? alignedGenomes.first : alignedGenomes.second;
    auto &id2Node = oneMoreFit ? id2Node1 : id2Node2;
    for (auto it = moreFitAG.disjoint.begin(); ;++it) {
        if (it == moreFitAG.disjoint.end()) {
            it = moreFitAG.excess.begin();
        }
        if (it == moreFitAG.excess.end()) {
            break;
        }
        ConnectionGene newConnGene = **it;
        if (!newConnGene.enabled && !deps->random()->stillDisabled()) {
            newConnGene.enabled = true;
        }
        newConnSet.insert(newConnGene);
        if(id2Node.find(newConnGene.from) == id2Node.end() || id2Node.find(newConnGene.to) == id2Node.end()) {
            std::cout << "Node not found on" << std::endl;
            std::cout << newConnGene << std::endl;
        }
        newNGSet.emplace(*id2Node[newConnGene.from]);
        newNGSet.emplace(*id2Node[newConnGene.to]);
    }

    bool newUseBias = false;
    size_t newInputs = 0;
    size_t newOutputs = 0;
    for (auto &ng : newNGSet) {
        switch(ng.nodeType) {
            case NodeGene::NodeType::Bias:
                newUseBias = true;
                break;
            case NodeGene::NodeType::Input:
                ++newInputs;
                break;
            case NodeGene::NodeType::Output:
                ++newOutputs;
                break;
            default:
                break;
        }
    }
    std::list<NodeGene> newNodeGenes(newNGSet.begin(), newNGSet.end());
    std::list<ConnectionGene> newConnGenes(newConnSet.begin(), newConnSet.end());
    return Genome(newInputs, newOutputs, newId, deps, std::move(newNodeGenes), std::move(newConnGenes), newUseBias);
}

size_t Genome::id() const {
    return _id;
}

void Genome::setId(size_t id) {
    _id = id;
}

const std::list<ConnectionGene>& Genome::getConnectionGenes() const {
    return connectionGenes;
}

const std::list<NodeGene>& Genome::getNodeGenes() const {
    return nodeGenes;
}

std::ostream &operator<<(std::ostream &out, const Genome &genome) {
    out << "ID: " << genome.id() << std::endl;
    for(auto it = genome.nodeGenes.begin(); it != genome.nodeGenes.end(); ++it) {
        if(it != genome.nodeGenes.begin()) {
            out << ", ";
        }
        out << *it;
    }
    out << std::endl;

    for(auto it = genome.connectionGenes.begin(); it != genome.connectionGenes.end(); ++it) {
        if(it != genome.connectionGenes.begin()) {
            out << ", ";
        }
        out << *it;
    }
    out << std::endl;

    return out;
}

NodeGene::NodeGene(NodeType nodeType, size_t id) : nodeType(nodeType), id(id) {}

bool NodeGene::operator==(const NodeGene &other) const {
    return id == other.id && nodeType == other.nodeType;
}
bool NodeGene::operator!=(const NodeGene &other) const {
    return !operator==(other);
}

std::ostream &operator<<(std::ostream &out, const NodeGene &ng) {
    out << '(' << ng.id << ", ";
    switch (ng.nodeType) {
        case NodeGene::NodeType::Input:
            out << "Input";
            break;
        case NodeGene::NodeType::Hidden:
            out << "Hidden";
            break;
        case NodeGene::NodeType::Output:
            out << "Output";
            break;
        case NodeGene::NodeType::Bias:
            out << "Bias";
            break;
        case NodeGene::NodeType::Default:
            out << "Default";
            break;
        default:
            std::cerr << "Unknown node type" << std::endl;
            break;
    }

    out << ')';

    return out;
}

ConnectionGene::ConnectionGene(size_t from, size_t to, double weight, bool enabled, size_t innovationNumber)
    : from(from), to(to), weight(weight), enabled(enabled), innovationNumber(innovationNumber) {}

bool ConnectionGene::operator==(const ConnectionGene &other) const {
    return from == other.from && to == other.to && innovationNumber == other.innovationNumber
        && weight == other.weight && enabled == other.enabled;
}

bool ConnectionGene::operator!=(const ConnectionGene &other) const {
    return !operator==(other);
}

std::ostream &operator<<(std::ostream &out, const ConnectionGene &cg) {
    out << '(' << cg.from << ", " << cg.to << ", " << cg.weight << ", "
        << cg.enabled << ", " << cg.innovationNumber << ')'; 
    return out;
}