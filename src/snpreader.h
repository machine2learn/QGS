#ifndef INCL_SNPREADER_H
#define INCL_SNPREADER_H

#include "gzfile.h"
#include <string>
#include <iosfwd>
#include <vector>
#include <sstream>
#include <utility> // swap
#include <unordered_map>

class SNPreader {
  
  public:
  
  struct Locus {
    long long chr = 0;
    std::size_t pos = 0;
    std::string id, ref, alt, qual, filter, info_str, format, chr_str;
    std::vector<std::string> palt;
    std::unordered_map<std::string, std::string> pinfo;
    double maf = 0;
    std::vector<float> data_ds;
    bool switch_ar = false;

    void clear() {
      maf = 0;
      switch_ar = false;
      palt.clear();
      data_ds.clear();
    }
    
    void parse_alt() {
      palt.clear();
      std::istringstream iss(alt);
      std::string tmp;
      while (std::getline(iss, tmp, ','))
        palt.push_back(std::move(tmp));
    }
    
    void parse_info() {
      std::istringstream iss(info_str);
      std::string tmp;
      while (std::getline(iss, tmp, ';')) {
        std::size_t const p = tmp.find('=');
        if (p == std::string::npos || p == tmp.size())
          continue;
        pinfo[tmp.substr(0, p)] = tmp.substr(p + 1);
      }
    }
    
    void switch_alt_ref() {
      switch_ar = true;
      if (!palt.empty())
        std::swap(palt[0], ref);
    }
  };
  
  protected:
    std::string d_fname;
    std::istringstream d_buffer;
    GZfile d_file;
    std::vector<std::string> d_header, d_sample;
    std::size_t d_num_samples;
    bool d_hard_calls, d_allow_missings;
  
  public:
    SNPreader(std::string const &);
    virtual ~SNPreader() = default;
    operator bool() const;
    std::size_t num_samples() const;
    std::string const & sample_id(std::size_t) const;
    
    virtual bool deep_read(Locus &) = 0;
    virtual void parse_line(Locus &) = 0;

  protected:
    virtual bool parse_header() = 0;
};

inline SNPreader::operator bool() const {
  return d_file;
}

inline std::size_t SNPreader::num_samples() const {
  return d_num_samples;
}

inline std::string const & SNPreader::sample_id(std::size_t idx) const {
  return d_sample[idx];
}

SNPreader & operator>>(SNPreader &, SNPreader::Locus &);

std::ostream & operator<<(std::ostream &, SNPreader::Locus const &);

inline bool operator<=(SNPreader::Locus const & lhs, SNPreader::Locus const & rhs) {
  return lhs.chr < rhs.chr || (lhs.chr == rhs.chr && lhs.pos <= rhs.pos);
}
    
  
#endif
