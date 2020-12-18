/*
  GNU General Public License v3 (GPLv3)
  Copyright (c) 2016-2020 Gido Schoenmacker
  Complete License at https://github.com/machine2learn/QGS/blob/master/LICENSE.md
*/

#include "snpreader.h"
#include <sstream>
#include <ostream>
#include <iostream>
#include "log.h"

SNPreader::SNPreader(std::string const & fname)
 :
  d_fname{fname},
  d_file(fname),
  d_num_samples{0},
  d_hard_calls{false},
  d_allow_missings{false}
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
