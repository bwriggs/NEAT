#pragma once

//Base class for function types

#include <cstdlib>

class BaseFunctor {
public:
    size_t inputs() const;
    size_t outputs() const ;
    BaseFunctor(size_t inputs, size_t ouputs);
private:
    size_t in, out;
};
