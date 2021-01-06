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
