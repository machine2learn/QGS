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

#ifndef INCL_QGS_H
#define INCL_QGS_H

#include "command_line_arguments.h"
#include "gzfile.h"
#include "genblock.h"
#include "plinkdosagereader.h"
#include "plinkbedreader.h"
#include "vcfreader.h"
#include "gene_score.h"
#include <string>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp> // ends_with
#include <iomanip>
#include <iostream>
#include <cmath> // NAN
#include "log.h"

inline std::unique_ptr<SNPreader> open_genetics_file(std::vector<std::string> const & files, bool hard_calls = 0, bool allow_missings = 0) {

  if (files.empty()) {
    LOG(QGS::Log::FATAL) << "File name not provided.\n";
    std::exit(-1);
  }

  if (files.size() > 1) {
    LOG(QGS::Log::VERBOSE) << "Assuming files `" << files.at(0) << "` etc are plink dosage format.\n";
    return std::unique_ptr<SNPreader>{std::unique_ptr<Plinkdosagereader>(new Plinkdosagereader(files))};
  }

  if (boost::algorithm::ends_with(files.at(0), ".dosage") || boost::algorithm::ends_with(files.at(0), ".dosage.gz")) {
    LOG(QGS::Log::VERBOSE) << "Assuming files `" << files.at(0) << "` etc are plink dosage format.\n";
    return std::unique_ptr<SNPreader>{std::unique_ptr<Plinkdosagereader>(new Plinkdosagereader(files.at(0)))};
  }

  if (boost::algorithm::ends_with(files.at(0), ".bed") || boost::algorithm::ends_with(files.at(0), ".bed.gz")) {
    LOG(QGS::Log::VERBOSE) << "Assuming `" << files.at(0) << "` input file is PLINK BED format.\n";
    return std::unique_ptr<SNPreader>{std::unique_ptr<Plinkbedreader>(new Plinkbedreader(files.at(0), allow_missings))};
  }
  
  LOG(QGS::Log::VERBOSE) << "Assuming `" << files.at(0) << "` input file is VCF format (default).\n";
  return std::unique_ptr<SNPreader>{std::unique_ptr<VCFreader>(new VCFreader(files.at(0), hard_calls, allow_missings))};
}

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
