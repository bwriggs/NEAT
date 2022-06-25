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

class AllDeps : virtual public RandomDeps, virtual public GenomeDeps {
public:
    virtual Random *random() = 0;
    virtual InnovationTracker *innovationTracker() = 0;
    virtual Parameters *parameters() = 0;
};

void stackTrace();

#define assertrc(expr) {if(expr){}else{stackTrace(); assert(expr);}}

