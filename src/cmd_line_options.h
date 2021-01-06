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

#ifndef INCL_QGS_CMD_LINE_OPTIONS_H
#define INCL_QGS_CMD_LINE_OPTIONS_H

#include "log.h"

namespace QGS {

inline std::unordered_map<std::string, std::vector<std::string>> 
  parse_cmd_line_options(size_t argc, char ** argv) {

  std::unordered_map<std::string, std::vector<std::string>> options;

  std::string current_option;
  for (std::size_t idx = 1; idx != argc; ++idx) {
    std::string token(argv[idx]);
    if (token.substr(0, 2) == "--") {
      current_option = token.substr(2);
      if (options.find(current_option) != options.end())
        LOG(QGS::Log::WARNING) << "Command line option `" << token
          << "` specified more than once. Ignoring all but last.\n";
      options[current_option] = {};
    }
    else {
      options[current_option].push_back(token);
    }
  }

  return options;
}

}

#endif
