/*******************************************************************************
 *  - FILE:  mmoreseqs_pylib.cpp
 *  - DESC:  Python bindings via pybind11 for MMOREseqs
 *******************************************************************************/

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "easel.h"
#include "pybind_classes.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(mmoreseqs_pylib, m) {
  m.doc() = R"pbdoc(
        MMOREseqs 
        -----------------------
        .. currentmodule:: mmoreseqs
        .. autosummary::
           :toctree: _generate
           run
    )pbdoc";

  m.def(
      "run", [](py::object& pyargv11) {
        PipelineCpp::RunAppFromCommandLine(pyargv11);
      },
      R"pbdoc(
          Parses commandline and runs pipeline with respect to arguments.
      )pbdoc");

  // py::class_<ArgsCpp>
  //     args(m, "ARGS", "List of arguments for determining pipline..");

  // py::class_<WorkerCpp> worker(m, "WORKER", "Worker object to for data to be modified by pipeline.");

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
