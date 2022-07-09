#include "LibNeat.h"
#include "Parameters.h"
#include "Random.h"
#include "InnovationTracker.h"

#include <execinfo.h>
#include <iostream>

Exception::Exception(std::string msg) : msg(msg) {}

std::string Exception::what() const {
    return msg;
}

void stackTrace () {
  void* callstack[1024];
  int i, frames = backtrace(callstack, 1024);
  char** strs = backtrace_symbols(callstack, frames);
  for (i = 0; i < frames; ++i) {
    std::cout << strs[i] << std::endl;
  }
  free(strs);
}


TestDeps::TestDeps(int seed) {
    innovTrack = new InnovationTracker();
    params = new Parameters();
    *params = defaultParams();
    rand = seed < 0 ? (new Random(this)) : (new Random(this, seed));
}

TestDeps::~TestDeps() {
    delete params;
    delete rand;
    delete innovTrack;
}

Random *TestDeps::random() { return rand; }

Parameters *TestDeps::parameters() { return params; }

InnovationTracker *TestDeps::innovationTracker() { return innovTrack; }

