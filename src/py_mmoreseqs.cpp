/*******************************************************************************
 *  - FILE:      py_mmoreseqs.cpp
 *  - DESC:    Python bindings via pybind11 for MMOREseqs
 *******************************************************************************/

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "py_class.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
  return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(py_mmoreseqs, m) {
  m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: mmoreseqs
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

  m.def("add", &add, R"pbdoc(
        Add two numbers
        Some other explanation about the add function.
    )pbdoc");

  m.def(
      "subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
        Some other explanation about the subtract function.
    )pbdoc");

  py::class_<BuildInfo> build_info(m, "BUILD_INFO", "Information about current build of MMOREseqs.");
  build_info
      .def("get_version", &BuildInfo::GetVersion);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
