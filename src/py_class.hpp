/*******************************************************************************
 *  - FILE:  py_class.cpp
 *  - DESC:  C++ class wrappers around C code for use by pybind.
 *******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "macros/_macros.h"

class BuildInfo {
 public:
  static std::string GetVersion() {
    std::stringstream os;
    os << BUILD_VERSION;
    return os.str();
  }
};
