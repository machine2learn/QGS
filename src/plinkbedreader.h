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
    Plinkbedreader(std::string const &);
    bool deep_read(Locus &) final;
    void parse_line(Locus &) final;

  protected:
    bool parse_header() final;
    
  private:
    std::string find_file(std::string const &) const;
};

#endif
