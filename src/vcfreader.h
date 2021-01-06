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

#ifndef INCL_VCFREADER_H
#define INCL_VCFREADER_H

#include "snpreader.h"

class VCFreader : public SNPreader {
  
  private:
    std::string d_format;

  public:
    VCFreader(std::string const &, bool = false, bool = false);
    bool deep_read(Locus &) final;
    void parse_line(Locus &) final;

  protected:
    bool parse_header() final;
    
  private:
    bool read_gt(Locus &);
    bool read_ds(Locus &);
    bool read_plink(Locus &);
};

#endif
