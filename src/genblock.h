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

#ifndef INCL_GENBLOCK_H
#define INCL_GENBLOCK_H

#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>

struct Genblock {
  std::string source, type;
  long long chr = 0;
  std::size_t start = 0, stop = 0;
  std::size_t in_sample = 0, in_ref = 0, overlap = 0;
  char score = '?', strand = '?', phase = '?';

  std::unordered_map<std::string, std::string> attr = {};
};

inline std::string str(Genblock const & gene) {
  std::string s =
    std::to_string(gene.chr) + ":" + 
    std::to_string(gene.start) + "-" + 
    std::to_string(gene.stop) + 
    " [score=" + gene.score +
    ";strand=" + gene.strand +
    ";type=" + gene.type +
    ";source=" + gene.source +
    ";phase=" + gene.phase;
    
    for (auto const & attr : gene.attr)
      s += ";" + attr.first + "=" + attr.second;
      
    s += "]";
    
    return s;
}

inline std::ostream & operator<<(std::ostream & out, Genblock const & g) {
  return out << str(g);
}

inline std::istream & operator>>(std::istream & in, Genblock & g) {
  
  g = {};
  
  std::string line;
  if (!std::getline(in, line)) {
    return in;
  }
  
  if (line.substr(0, 2) == "##") {
    g.type = "comment";
    return in;
  }
  
  std::istringstream iss(line);

  std::string chr;
  if (!(iss >> chr >> g.source >> g.type)) {
    in.setstate(std::ios::failbit);
    return in;
  }

  try {
    if (chr == "chrX")
      g.chr = 23;
    else if (chr == "chrY")
      g.chr = 24;
    else if (chr == "chrM")
      g.chr = 25;
    else
      g.chr = chr.length() > 2
        ? std::stoull(chr.substr(3))
        : std::stoull(chr);
  }
  catch (std::invalid_argument const & e) {
    in.setstate(std::ios::failbit);
    return in;
  }

  if (!(iss >> g.start >> g.stop >> g.score >> g.strand >> g.phase)) {
    in.setstate(std::ios::failbit);
    return in;
  }
  
  g.attr["chr"] = chr;
  g.attr["source"] = g.source;
  g.attr["type"] = g.type;

  std::string key, val;
  while(iss >> key >> val)
    g.attr[key] = 
    (val.length() > 3 && val[0] == '"' &&
     val[val.length()-2] == '"' && val[val.length()-1] == ';')
     ?
       val.substr(1, val.length() - 3)
     :
       val;
  
  // set a default name in case none is provided
  if (g.attr.find("gene_name") == g.attr.end()) {
    g.attr["gene_name"] = std::to_string(g.chr) + ":" + 
      std::to_string(g.start) + "-" +
      std::to_string(g.stop);
  }

  // set a default id in case none is provided
  if (g.attr.find("gene_id") == g.attr.end()) {
    g.attr["gene_id"] = g.attr["gene_name"];
  }

  return in;
}

#endif
