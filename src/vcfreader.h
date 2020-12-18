/*
  GNU General Public License v3 (GPLv3)
  Copyright (c) 2016-2020 Gido Schoenmacker
  Complete License at https://github.com/machine2learn/QGS/blob/master/LICENSE.md
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
