/*
  GNU General Public License v3 (GPLv3)
  Copyright (c) 2016-2020 Gido Schoenmacker
  Complete License at https://github.com/machine2learn/QGS/blob/master/LICENSE.md
*/

#ifndef INCL_COMMAND_LINE_ARGUMENTS_H
#define INCL_COMMAND_LINE_ARGUMENTS_H

#define QGS_VERSION "1.0+ Devel (Dec 2020)"

#include <unordered_map>
#include <vector>
#include <string>
#include "log.h"

namespace QGS {

std::unordered_map<std::string, std::vector<std::string>> command_line_arguments(int, char **);

std::unordered_map<std::string, std::string> create_gtf_filter(std::unordered_map<std::string, std::vector<std::string>> const &);

std::unordered_map<std::string, bool> create_snp_filter(std::unordered_map<std::string, std::vector<std::string>> const &);

}

#endif
