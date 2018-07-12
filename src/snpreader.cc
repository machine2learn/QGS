#include "snpreader.h"
#include <sstream>
#include <ostream>
#include <iostream>
#include "log.h"

SNPreader::SNPreader(std::string const & fname)
 :
  d_fname{fname},
  d_file(fname),
  d_num_samples{0}
{
}

SNPreader & operator>>(SNPreader & in, SNPreader::Locus & locus) {
  in.parse_line(locus);
  return in;
}

std::ostream & operator<<(std::ostream & out, SNPreader::Locus const & l) {
  return out << l.chr << ':' << l.pos << " (" << l.id << "): "
             << "[ref=" << l.ref << "; alt=" << l.alt << "; qual="
             << l.qual << "; filter=" << l.filter << "; info="
             << l.info_str << "; format=" << l.format << "; flip="
             << l.switch_ar << ";]";
}
