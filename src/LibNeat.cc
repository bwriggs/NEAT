#include "LibNeat.h"
#include "execinfo.h"
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


