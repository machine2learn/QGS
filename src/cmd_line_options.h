#ifndef INCL_QGS_CMD_LINE_OPTIONS_H
#define INCL_CMD_LINE_OPTIONS_H

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
