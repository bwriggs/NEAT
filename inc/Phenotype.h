#pragma once

// the actual function described by a genome

#include "TaskFunctor.h"

#include <Eigen/Dense>
#include <unordered_map>

class Genome;
class Parameters;

//layer space index, this is the data we are going to have
//while evaluating an input vector
//used to record outputs and recurrent links
struct LSIDX {
    size_t layerIdx;
    size_t idx;
    LSIDX(size_t, size_t);
    LSIDX();
};

//connection from one layer
//to a previous layer, used to get
//all "extra" (i.e. those from future layers)
//inputs to a node
struct RecurrentLink {
    size_t fromIdx;
    double weight;
    LSIDX lsidx;
    RecurrentLink(size_t, double, LSIDX lsidx);
};


//used to record outputs and remove them
//from further processing as their layer's are
//reached
struct OutputInfo {
    size_t vecIdx;
    size_t nodeIdx;
    OutputInfo(size_t, size_t);
};

class Phenotype : public TaskFunctor {
public:
    Phenotype(const Genome&);
    //evaluate at inputs
    std::vector<double> operator()(std::vector<double>) override;

    size_t id();

    friend std::ostream &operator<<(std::ostream&, const Phenotype&);

private:
    //matrices of weights
    std::vector<Eigen::MatrixXd> weights;
    //the values added by recurrent connections
    std::vector<Eigen::VectorXd> recurrentVals;
    //map from layer id to recurrent contributions originating there
    std::unordered_map<size_t, std::vector<RecurrentLink>> recurrentLinks;
    //map from layer id to outputs (so they can be recorded and "shaved off")
    std::unordered_map<size_t, std::vector<OutputInfo>> layer2Outputs;

    bool useBias;
    size_t _id;
};
