#ifndef INCL_COMMAND_LINE_ARGUMENTS_H
#define INCL_COMMAND_LINE_ARGUMENTS_H

#define QGS_VERSION "1.4 plink-workaround"

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
