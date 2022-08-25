/*******************************************************************************
 *  - FILE:  py_mmoreseqs.cpp
 *  - DESC:  Python bindings via pybind11 for MMOREseqs
 *******************************************************************************/

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// #include "pybind_classes.hpp"
// #include "easel.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
  return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(py_mmoreseqs, m) {
  m.doc() = R"pbdoc(
        MMOREseqs 
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

  // py::class_<ArgsCpp> args(m, "ARGS", "List of pipeline arguments.");

  // py::class_<WorkerCpp> worker(m, "WORKER", "Worker object to be processed by pipeline.");

  // py::class_<PipelineCpp> pipeline(m, "PIPELINE", "Pipeline object for conducting searches and other routines.");

  // py::class_<BuildInfo> build_info(m, "BUILD_INFO", "Information about current build of MMOREseqs.");
  // build_info
  //     .def("get_version", &BuildInfo::GetVersion);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
