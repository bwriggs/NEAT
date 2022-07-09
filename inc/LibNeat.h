#pragma once

#include <string>

class Random;
class InnovationTracker;
class Parameters;

struct Exception {
    std::string msg;
    std::string what() const;
    Exception(std::string = "");
};

class RandomDeps {
public:
    virtual Parameters *parameters() = 0;
};

class GenomeDeps {
public:
    virtual Random *random() = 0;
    virtual InnovationTracker *innovationTracker() = 0;
    virtual Parameters *parameters() = 0;
};

enum class Mutation { PerturbWeight, NewWeight, AddConnection, AddNode, None };

class AllDeps : virtual public RandomDeps, virtual public GenomeDeps {
public:
    virtual Random *random() = 0;
    virtual InnovationTracker *innovationTracker() = 0;
    virtual Parameters *parameters() = 0;
};

class TestDeps : public AllDeps {
public:
    TestDeps(int = -1);
    ~TestDeps();
    Random *random() override;
    InnovationTracker *innovationTracker() override;
    Parameters *parameters() override;
private:
    Random *rand;
    InnovationTracker *innovTrack;
    Parameters *params;
};

void stackTrace();

struct ConnId {
    size_t from, to;
    ConnId(size_t from, size_t to) : from(from), to(to) {}
    bool operator==(const ConnId &other) const {
        return (from == other.from && to == other.to)
        #ifdef __NEAT_UNIDIRECTIONAL__
            || (to == other.from && from == other.to)
        #endif
        ;
    }
};

namespace std {
    template<>
    struct hash<ConnId> {
        size_t operator()(const ConnId &connId) const {
            return (hash<size_t>{}(connId.from) << 5) ^ hash<size_t>{}(connId.to);
        }
    };
}

#define assertrc(expr) {if(expr){}else{stackTrace(); assert(expr);}}

