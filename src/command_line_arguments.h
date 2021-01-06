/* 
 * This file is part of the QGS distribution https://github.com/machine2learn/QGS/.
 * Copyright (c) 2016-2020 Gido Schoenmacker
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
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
