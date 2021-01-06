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

#ifndef INCL_PLINKBEDREADER_H
#define INCL_PLINKBEDREADER_H

#include "snpreader.h"
#include <string>
#include <vector>

class Plinkbedreader : public SNPreader {
  
  private:
    std::string d_bim_fname, d_fam_fname;
    GZfile d_bim_file, d_fam_file;
    std::size_t d_num_pos_read, d_bytes_per_locus;

  public:
    Plinkbedreader(std::string const &, bool = false);
    bool deep_read(Locus &) final;
    void parse_line(Locus &) final;

  protected:
    bool parse_header() final;
    
  private:
    std::string find_file(std::string const &) const;
};

#endif
