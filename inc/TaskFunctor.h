#pragma once

// just a custom Functor

#include <vector>

#include "BaseFunctor.h"

class TaskFunctor : public BaseFunctor {
public:
    TaskFunctor(size_t, size_t);
    virtual std::vector<double> operator()(std::vector<double>) = 0;
};
