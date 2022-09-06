/*******************************************************************************
 *  - FILE:  py_class.hpp
 *  - DESC:  C++ class wrappers around MMORESEQS C codebase for use by pybind.
 *******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "macros/_macros.h"
#include "objects/_objects.h"
#include "parsers/_parsers.h"
#include "pipelines/_pipelines.h"
#include "application.h"

// PyText cast depending on version of Python.
#if PY_VERSION_HEX < 0x03000000
#define MyPyText_AsString PyString_AsString
#else
#define MyPyText_AsString PyUnicode_AsUTF8
#endif

namespace py = pybind11;

// Store C-style commandline args, built from Python-style commandline args.
class CommandLineArgsCpp {
 public:
  CommandLineArgsCpp(py::object& pyargv11) {
    PyObject* pyargv = pyargv11.ptr();
    if (PySequence_Check(pyargv)) {
      Py_ssize_t sz = PySequence_Size(pyargv);
      argc_ = (int)sz;

      argv_ = (char**)malloc(sz * sizeof(char*));
      for (Py_ssize_t i = 0; i < sz; ++i) {
        PyObject* item = PySequence_GetItem(pyargv, i);
        argv_[i] = (char*)MyPyText_AsString(item);
        Py_DECREF(item);
        if (!argv_[i] || PyErr_Occurred()) {
          free(argv_);
          argv_ = nullptr;
          break;
        }
      }
    }

    if (!argv_) {
      // fprintf(stderr,  "argument is not a sequence of strings\n");
      // return;
      if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "could not convert input to argv");
      }
      throw py::error_already_set();
    }
  }

  ~CommandLineArgsCpp() {
    free(argv_);
  }

  void Print() {
    fprintf(stderr, "COMMANDLINE ARGS:\n");
    for (int i = 0; i < argc_; ++i) {
      fprintf(stderr, "%s\n", argv_[i]);
    }
  }

  int argc_ = 0;
  char** argv_ = NULL;
};

class ArgsCpp {
 public:
  ArgsCpp()
      : is_owned_by_worker_(true) {
    args_ = ARGS_Create();
  }

  ArgsCpp(ARGS* args)
      : is_owned_by_worker_(false) {
    args_ = args;
  }

  ~ArgsCpp() {
    if (is_owned_by_worker_) {
      ARGS_Destroy(args_);
    }
  }

  void ParseCommandLine(int argc, char* argv[]) {
  }

  bool is_owned_by_worker_;
  ARGS* args_ = nullptr;

  // static std::unordered_map<std::string, std::pair<type, void*>> flag_map = {
  //     {"--help", {void, nullptr}},
  //     {""}};
};

class WorkerCpp {
 public:
  WorkerCpp() {
    worker_ = WORKER_Create();
    WORKER_Init(worker_);
    args_ = ArgsCpp(worker_->args);
  }

  ~WorkerCpp() {
    WORKER_Destroy(worker_);
  }

  ArgsCpp& GetArgs() {
    return args_;
  }

  WORKER* worker_ = nullptr;
  ArgsCpp args_ = nullptr;
};

class PipelineCpp {
 public:
  PipelineCpp() {
    worker_ = WorkerCpp();
  }

  WorkerCpp& GetWorker() {
    return worker_;
  }

  void SetMode(const std::string& pipeline_mode) {
  }

  static void RunAppFromCommandLine(py::object& pyargv11) {
    CommandLineArgsCpp cmd_args = CommandLineArgsCpp(pyargv11);
    STATUS_FLAG status = APPLICATION_Run(cmd_args.argc_, cmd_args.argv_);
  }

  WorkerCpp worker_;
};

class BuildInfo {
 public:
  static std::string GetVersion() {
    std::stringstream os;
    os << BUILD_VERSION;
    return os.str();
  }
};
