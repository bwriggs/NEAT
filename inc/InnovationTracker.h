#pragma once

//the logic to ensure the same
//structural innovations get the
//same historical origin

#include <cstdlib>
#include <unordered_map>

//only connections have innovation
//numbers, they can be identified
//by their from and to fields
struct InnovId {
    size_t from, to;
    InnovId(size_t,size_t);

    bool operator==(const InnovId &) const;
};

//define a hash for this structure (used by unordered_map)
namespace std {
    template<>
    struct hash<InnovId> {
        size_t operator()(const InnovId &id) const {
            //good enough? .num_buckets() or something to check
            return (std::hash<size_t>{}(id.from) << 5) ^ std::hash<size_t>{}(id.to);
        }
    };
};

class ConnectionGene;
class InnovationTracker {
public:
    InnovationTracker();
    //get the innovation number of a structural change
    void assignInnov(ConnectionGene &);
    //forget innovations (called at the start
    //of a new generation)
    void clear();

private:
    //reference to a global innovation number
    size_t nextInnovation;
    //structural innovations from this generation
    //and their innovationNumbers (cleared by clear())
    std::unordered_map<InnovId, size_t> innovs;
};
