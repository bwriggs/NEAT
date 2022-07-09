#include "Phenotype.h"
#include "Parameters.h"
#include "Genome.h"

#include <unordered_map>
#include <deque>
#include <queue>
#include <iostream>

//the weights in one row of a weight matrix
//along with node associated with the row's output
struct WeightRow {
    std::vector<double> weights;
    const NodeGene *ng;
    WeightRow(size_t size, const NodeGene *ng) : weights(size), ng(ng) {}
};

//make phenotype (efficiently callable recurrent neural network)
Phenotype::Phenotype(const Genome &genome) 
: TaskFunctor(genome.inputs, genome.outputs), useBias(genome.useBias), _id(genome.id()) {
    //general idea is to calculate just like a feed-forward neural net EXCEPT
    //  we add in the recurrent contributions before sigmoiding and multiplying
    //    by each weight matrix
    // // DO WE ? // //  we "shave off" outputs as they are computed (i.e. record them in an
    // // DO WE ? // //    output std::vector and trucate the Eigen::VectorXd)
    //  we zero out out recurrent inputs and add our recurrent outputs
    //    to the recurrent inputs of previous layers (to be used in the next call)

    //make map from nodes' ids to their outputs
    //so we can do a bredth first (i.e. layer by layer) "search"
    std::unordered_map<size_t, std::vector<const ConnectionGene*>> node2Outputs;
    for(const ConnectionGene &cg : genome.connectionGenes) {
        if(!cg.enabled) continue;
        auto it = node2Outputs.find(cg.from);
        if(it == node2Outputs.end()) {
            node2Outputs[cg.from] = std::vector<const ConnectionGene*>{&cg};
        } else {
            it->second.push_back(&cg);
        }
    }

    //make a map from node ids to nodes
    std::unordered_map<size_t, const NodeGene*> id2Node;
    for (const NodeGene &ng : genome.nodeGenes) {
        id2Node[ng.id] = &ng;
    }

    //mapping from the id of a node to its layer
    //and index in its vector
    //also functions as a visited structure
    std::unordered_map<size_t, LSIDX> id2Idx;
    //nodes in current layer
    std::queue<const NodeGene*> layer;
    //initialize with inputs
    for (const NodeGene &ng : genome.nodeGenes) {
        if(ng.nodeType == NodeGene::NodeType::Input || ng.nodeType == NodeGene::NodeType::Bias) {
            layer.push(&ng);
            LSIDX lsidx(0, ng.id);
            id2Idx.emplace(ng.id, lsidx);
        }
    }
    recurrentVals.emplace_back(layer.size());
    for (size_t l = 1; layer.size() > 0; ++l) {
        //want to produce three things per layer
        //1) recurrent outputs to add to the recurrentVals
        //   of previous layers
        //2) a weight matrix to compute feed-forward outputs
        //   from this layer
        //3) a list of outputs to be recorded from the vector 
        //   after adding recurrentVals, sigmoiding
        //   and multiplying by the weight matrix
        //*additionally have to maintain id2Idx

        //columns of the weight matrix
        //stored as deque so inputs/hidden can be pushed to the front
        //and outputs can be pushed to the back (so they can be efficiently removed)
        //in Eigen::VectorXd form
        std::deque<WeightRow> weightMatrix;
        //for constant access to the middle of the deque
        std::unordered_map<size_t, std::vector<double>*> id2Row;

        //store in variable since layer will grow
        size_t layerSize = layer.size();
        
        //fill the weight matrix deque and record recurrent links
        for(size_t i = 0; i < layerSize; ++i) {
            const NodeGene *from = layer.front(); layer.pop();
            for(const ConnectionGene *cg : node2Outputs[from->id]) {
                const NodeGene *to = id2Node[cg->to];
                auto prevRes = id2Idx.find(to->id);
                if (prevRes != id2Idx.end()) {
                    //recurrent
                    //want to register connection
                    LSIDX fromLSIDX = id2Idx[from->id];
                    LSIDX toLSIDX = prevRes->second;
                    RecurrentLink rl(fromLSIDX.idx, cg->weight, toLSIDX);
                    auto links = recurrentLinks.find(fromLSIDX.layerIdx);
                    if (links == recurrentLinks.end()) {
                        recurrentLinks[fromLSIDX.layerIdx] = std::vector<RecurrentLink>{rl};
                    } else {
                        links->second.push_back(rl);
                    }
                } else {
                    if (id2Row.find(to->id) == id2Row.end()) {
                        if(to->nodeType == NodeGene::NodeType::Output) {
                            weightMatrix.emplace_back(layerSize, to);
                            id2Row[to->id] = &weightMatrix.back().weights;
                        } else {
                            weightMatrix.emplace_front(layerSize, to);
                            id2Row[to->id] = &weightMatrix.front().weights;
                        }
                        //not sure if outputs should be pushed or not, lets see
                        layer.push(to);
                    }
                    id2Row[to->id]->operator[](i) = cg->weight;
                }  
            }
        }

        Eigen::MatrixXd newWeights(weightMatrix.size(), layerSize);
        size_t idx = 0;
        //update id2Idx and translate deque to MatrixXd
        for(auto rowIt = weightMatrix.begin(); rowIt != weightMatrix.end(); ++rowIt, ++idx) {
            newWeights.row(idx) = Eigen::VectorXd::Map(rowIt->weights.data(), rowIt->weights.size());
            id2Idx.emplace(rowIt->ng->id, LSIDX(l, idx));
        }
        //produce outputs
        for(auto rowIt = weightMatrix.rbegin(); rowIt != weightMatrix.rend() 
        && rowIt->ng->nodeType == NodeGene::NodeType::Output; ++rowIt) {
            OutputInfo oi(--idx, rowIt->ng->id);
            auto outputRes = layer2Outputs.find(l);
            if(outputRes == layer2Outputs.end()) {
                layer2Outputs[l] = std::vector<OutputInfo>{oi};
            } else {
                outputRes->second.push_back(oi);
            }
        }

        if(newWeights.size() > 0) {
            weights.push_back(newWeights);
            recurrentVals.emplace_back(weightMatrix.size());
        }
    }
}

std::vector<double> Phenotype::operator()(std::vector<double> input) {
    //L layers, final layer is output only.*
    //*is there a way to disable connections such that this isn't true?
    //initialize Eigen vec with input vals
    //At first layer
    //  initialize input Eigen vector with inputs (optionally w bias)
    //  add first layer recurrent vals
    //  sigmoid (to get first layer outputs)
    //For each layer (including first)
    //  update our and previous recurrent inputs with
    //    our recurrent outputs
    //  record network outputs
    //  break if last layer
    //  generate next layer by applying weights, adding
    //    recurrent vals and sigmoiding
    //Return all outputs
    //In total
    //  L recurrent inputs
    //  L-1 weight matrices
    //  up to L-1 output lists

    if (useBias) input.push_back(1);
    Eigen::VectorXd vec = Eigen::VectorXd::Map(input.data(), input.size());
    // std::cout << "Before Recurr: " << vec << std::endl;
    vec += recurrentVals[0];
    // std::cout << "After Recurr: " << vec << std::endl;
    vec = vec.unaryExpr(&sigmoid);
    // std::cout << "After sigmoid: " << vec << std::endl;

    std::vector<double> output(outputs());
    for(size_t i = 0; ; ++i) {
        //update recurrent vals
        recurrentVals[i].setZero();
        auto recurrRes = recurrentLinks.find(i);
        if (recurrRes != recurrentLinks.end()) {
            for(const RecurrentLink& rl : recurrRes->second) {
                // std::cout << "Adding recurr: (" << rl.lsidx.layerIdx << ", " 
                    // << rl.lsidx.idx << ") --" << rl.weight << ">> " << rl.fromIdx << std::endl;
                recurrentVals[rl.lsidx.layerIdx][rl.lsidx.idx] += vec(rl.fromIdx) * rl.weight;
            }
        }

        //record outputs
        auto outputRes = layer2Outputs.find(i);
        if (outputRes != layer2Outputs.end()) {
            size_t trueNumInputs = inputs() + (useBias ? 1 : 0);
            for (const OutputInfo &oi : outputRes->second) {
                size_t outputIdx = oi.nodeIdx - trueNumInputs;
                // std::cout << "Recording output: " << outputIdx << "->" << oi.vecIdx << std::endl;
                output[outputIdx] = vec(oi.vecIdx);
            }
        }

        if(i >= weights.size()) {
            break;
        }

        vec.applyOnTheLeft(weights[i]);
        // std::cout << "After weights: " << vec << std::endl;
        vec += recurrentVals[i + 1];
        // std::cout << "After recurr: " << vec << std::endl;
        vec = vec.unaryExpr(&sigmoid);
        // std::cout << "After sigmoid: " << vec << std::endl;
    }

    return output;
}

size_t Phenotype::Phenotype::id() {
    return _id;
}

//layer space index, this is the data we are going to have
//while evaluating an input vector
//used to record outputs and recurrent links
LSIDX::LSIDX(size_t layerIdx, size_t idx) : layerIdx(layerIdx), idx(idx) {}
LSIDX::LSIDX() : LSIDX(0,0) {}

RecurrentLink::RecurrentLink(size_t fromIdx, double weight, LSIDX lsidx)
    : fromIdx(fromIdx), weight(weight), lsidx(lsidx) {}

OutputInfo::OutputInfo(size_t vecIdx, size_t nodeIdx) : vecIdx(vecIdx), nodeIdx(nodeIdx) {}

std::ostream &operator<<(std::ostream &out, const Phenotype &phenotype) {
    out << "Matrices (" << phenotype.weights.size() << ')' << std::endl;
    auto it = phenotype.weights.begin();
    if (it != phenotype.weights.end()) {
        out << *it << std::endl;
        ++it; 
    }
    for(; it != phenotype.weights.end(); ++it) {
        out << "------------" << std::endl;
        out << *it << std::endl;
    }
    out << std::endl;

    out << "Reccurent links (" << phenotype.recurrentLinks.size() << ')' << std::endl;
    for(auto kvp : phenotype.recurrentLinks) {
        out << kvp.first << " (" << kvp.second.size() << "): ";
        auto it = kvp.second.begin();
        if (it != kvp.second.end()) {
            out << it->fromIdx << ' ' << it->weight << "-> (" << it->lsidx.layerIdx << ',' << it->lsidx.idx << ')';
            ++it;
        }
        for(; it != kvp.second.end(); ++it) {
            out << ',' << it->fromIdx << ' ' << it->weight << "-> (" << it->lsidx.layerIdx << ',' << it->lsidx.idx << ')';
        }
        out << std::endl;
    }

    out << "Outputs (" << phenotype.layer2Outputs.size() << ')' << std::endl;
    for(auto kvp : phenotype.layer2Outputs) {
        out << kvp.first << ": ";
        auto it = kvp.second.begin();
        if (it != kvp.second.end()) {
            out << '(' << it->vecIdx << ',' << it->nodeIdx << ')';
            ++it;
        }
        for(; it != kvp.second.end(); ++it) {
            out << '(' << it->vecIdx << ',' << it->nodeIdx << ')';
        }
        out << std::endl;
    }

    return out;
}
