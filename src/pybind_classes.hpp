/*******************************************************************************
 *  - FILE:  py_class.hpp
 *  - DESC:  C++ class wrappers around MMORESEQS C codebase for use by pybind.
 *******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "macros/_macros.h"
#include "objects/_objects.h"
#include "parsers/_parsers.h"
#include "pipelines/_pipelines.h"

// class ArgsCpp {
//   //  public:
//   //   ArgsCpp()
//   //       : is_owned_by_worker_(true) {
//   //     args_ = ARGS_Create();
//   //   }

//   //   ArgsCpp(ARGS* args)
//   //       : is_owned_by_worker_(false) {
//   //     args_ = args;
//   //   }

//   //   ~ArgsCpp() {
//   //     if (is_owned_by_worker_) {
//   //       ARGS_Destroy(args_);
//   //     }
//   //   }

//   //   void ParseCommandLine(int argc, char* argv[]) {
//   //   }

//   //   void Set(const std::string& arg_name, const std::string& value) {
//   //   }

//   //   std::string Get(const std::string& arg) {
//   //     return arg;
//   //   }

//   //   ARGS* args_ = nullptr;
//   //   bool is_owned_by_worker_;

//   //   // static std::unordered_map<std::string, std::pair<type, void*>> flag_map = {
//   //   //     {"--help", {void, nullptr}},
//   //   //     {""}};
// };

// class WorkerCpp {
//   //  public:
//   //   WorkerCpp() {
//   //     worker_ = WORKER_Create();
//   //     WORKER_Init(worker_);
//   //     args_ = ArgsCpp(worker_->args);
//   //   }

//   //   ~WorkerCpp() {
//   //     WORKER_Destroy(worker_);
//   //   }

//   //   ArgsCpp& GetArgs() {
//   //     return args_;
//   //   }

//   //   void ToString() {
//   //     // WORKER_Dump(worker, stdout);
//   //   }

//   //   WORKER* worker_ = nullptr;
//   //   ArgsCpp args_ = nullptr;
// };

// class PipelineCpp {
//   //  public:
//   //   PipelineCpp() {
//   //     worker_ = WorkerCpp();
//   //   }

//   //   WorkerCpp& GetWorker() {
//   //     return worker_;
//   //   }

//   //   void SetMode(const std::string& pipeline_mode) {
//   //   }

//   //   void Run() {}

//   //   WorkerCpp worker_;
// };

class BuildInfo {
 public:
  static std::string GetVersion() {
    std::stringstream os;
    os << BUILD_VERSION;
    return os.str();
  }
};
