/*******************************************************************************
 *  - FILE:  cli_application.hpp
 *  - DESC:  Entry Point to Application. Parses commandline args and starts engine.
 *******************************************************************************/

#include "CLI11.hpp"

#include "objects/_objects.h"
#include "pipelines/_pipelines.h"
#include "parsers/_parsers.h"

// ** Subcommands

//
void EasySearchCommand(CLI::App& mmoreseqs_app);
//
void PrepSearchCommand(CLI::App& mmoreseqs_app);
//
void SearchCommand(CLI::App& mmoreseqs_app);
//
void MMORESearchCommand(CLI::App& mmoreseqs_app);
//
void MMSeqsSearchCommand(CLI::App& mmoreseqs_app);

// ** Options

//
void Add_GeneralOptions(CLI::App& cmd);
//
void Add_InternalOptions(CLI::App& cmd);
//
void Add_DebugOptions(CLI::App& cmd);
//
void Add_AllOptions(CLI::App& cmd);

// ** Helper Functions

//
void StoreFlagArgAsString(STR& flag_dest, std::string& flag_src);
//
void StoreFlagArgAsBool(bool& flag_dest, std::string& flag_src);
//
void StoreFlagArgAsInt(int& flag_dest, std::string& flag_src);
//
void StoreFlagArgAsDouble(int& flag_dest, std::string& flag_src);
