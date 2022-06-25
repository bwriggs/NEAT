#include "BaseFunctor.h"

BaseFunctor::BaseFunctor(size_t in, size_t out) : in(in), out(out) {}

size_t BaseFunctor::inputs() const { return in; }

size_t BaseFunctor::outputs() const { return out; }
