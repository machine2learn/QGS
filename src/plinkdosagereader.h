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

#ifndef INCL_PLINKDOSAGEREADER_H
#define INCL_PLINKDOSAGEREADER_H

#include "snpreader.h"
#include <string>
#include <vector>

class Plinkdosagereader : public SNPreader {
  
  private:
    std::string d_map_fname;
    std::vector<std::string> d_fnames;
    GZfile d_map_file;
    //bool d_force_hardcalls; always read dosages
    std::size_t d_linenr, d_filenr;
    float d_max;

  public:
    Plinkdosagereader(std::string const &);
    Plinkdosagereader(std::vector<std::string> const &);
    bool deep_read(Locus &) final;
    void parse_line(Locus &) final;

  protected:
    bool parse_header() final;
    
  private:
    std::string find_map_file() const;
    bool open_next();
};

#endif
