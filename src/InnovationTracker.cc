#include "InnovationTracker.h"
#include "Genome.h"

InnovationTracker::InnovationTracker() : nextInnovation(0) {}

void InnovationTracker::assignInnov(ConnectionGene &cg) {
    //try insert, returns pair<iterator, bool> = <in or inserted,inserted or not> 
    auto innov = innovs.emplace(InnovId(cg.from, cg.to), nextInnovation);
    if(innov.second) {
        ++nextInnovation;
    }

    cg.innovationNumber = innov.first->second;
    return;
}

void InnovationTracker::clear() {
    innovs.clear();
}

void InnovationTracker::reset() {
    nextInnovation = 0;
}

InnovId::InnovId(size_t from, size_t to) : from(from), to(to) {}

bool InnovId::operator==(const InnovId &other) const {
    return from == other.from && to == other.to;
}
