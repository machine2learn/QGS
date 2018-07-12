#ifndef INCL_VCFREADER_H
#define INCL_VCFREADER_H

#include "snpreader.h"

class VCFreader : public SNPreader {
  
  private:
    bool d_force_hardcalls;
    std::string d_format;

  public:
    VCFreader(std::string const &, bool = false);
    bool deep_read(Locus &) final;
    void parse_line(Locus &) final;

  protected:
    bool parse_header() final;
    
  private:
    bool read_gt(Locus &);
    bool read_ds(Locus &);
};

#endif
