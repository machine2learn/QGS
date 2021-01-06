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

#ifndef INCL_QGS_LOG_H
#define INCL_QGS_LOG_H

#include <ostream>
#include <iostream>
#include <sstream>

namespace QGS {

class Log {

  public:

    enum Loglevel { TRACE, DEBUG, VERBOSE, INFO, WARNING, FATAL };
    static Log & instance()
    {
        static Log instance;
        return instance;
    }
    
    void level(Log::Loglevel l) {
      s_level = l;
    }
    
    std::ostream & operator()(Loglevel l) {
      return l < s_level ? 
        d_oss : 
        (std::cout << "[" << tostring(l) << "] ");
    }

  private:
    Loglevel s_level;
    std::ostringstream d_oss;
    
    std::string tostring(Loglevel l) const {
      switch (l) {
        case TRACE : return "TRACE";
        case DEBUG : return "DEBUG";
        case VERBOSE: return "VERBOSE";
        case INFO : return "INFO";
        case WARNING : return "WARNING";
        case FATAL : return "FATAL";
      }
      return ""; // unreachable
    }
  
    Log() : s_level{INFO} { d_oss.setstate(std::ios_base::failbit); };
    Log(Log const &); // Not implemented = delete
    Log & operator=(Log const &); // Not implemented = delete
};

inline std::ostream & LOG(Log::Loglevel l) {
  return Log::instance()(l);
}

inline void LOGLVL(Log::Loglevel l) {
  QGS::Log::instance().level(l);
}

}

#endif
