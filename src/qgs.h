#ifndef INCL_QGS_H
#define INCL_QGS_H

#include "command_line_arguments.h"
#include "gzfile.h"
#include "genblock.h"
#include "plinkdosagereader.h"
#include "vcfreader.h"
#include "gene_score.h"
#include <string>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp> // ends_with
#include "log.h"

inline bool pass_gene_filter(Genblock const & gb, std::unordered_map<std::string, std::string> const & gtf_filter) {
  
  if (!gb.chr || !gb.start || !gb.stop) {
    LOG(QGS::Log::TRACE) << "Gene: incomplete information on region: skipping.\n";
    return false;
  }

  for (auto const & f : gtf_filter) {
    auto itt = gb.attr.find(f.first);
    if (itt == gb.attr.end()) {
      LOG(QGS::Log::TRACE) << "Region does not have info for gtf_filter `" << f.first << "`: skipping\n";
      return false;
    }
    if (itt->second != f.second) {
      LOG(QGS::Log::TRACE) << "Region fails gtf_filter `" << f.first << '=' << f.second << "` (has value `" << itt->second << "`): skipping\n";
      return false;
    }
  }
  
  return true;
}

inline bool operator<(Genblock const & gb, SNPreader::Locus const & l) {
  return gb.chr < l.chr || (gb.chr == l.chr && gb.stop < l.pos);
}

#endif
