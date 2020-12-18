/*
  GNU General Public License v3 (GPLv3)
  Copyright (c) 2016-2020 Gido Schoenmacker
  Complete License at https://github.com/machine2learn/QGS/blob/master/LICENSE.md
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
