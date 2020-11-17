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
