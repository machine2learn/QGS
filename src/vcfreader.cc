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

#include "vcfreader.h"
#include <sstream>
#include <ostream>
#include <iostream>
#include <cmath> // NAN
#include "log.h"


VCFreader::VCFreader(std::string const & fname, bool hard, bool miss)
 :
  SNPreader(fname)
{
  d_hard_calls = hard;
  d_allow_missings = miss;
  parse_header();
}

bool VCFreader::deep_read(VCFreader::Locus & l) {
  if (d_format == "PLINK")
    return read_plink(l);
  if (d_format == "GT")
    return read_gt(l);
  if (d_format == "DS")
    return read_ds(l);
  return false;
}

bool VCFreader::parse_header() {
  
  if (!d_file) {
    LOG(QGS::Log::FATAL) << "Cannot open input file `" << d_fname 
      << "`: for reading. Aborting.\n";
    std::exit(EXIT_FAILURE);
  }
  
  std::string line;
  bool found_format_tag = false;
  while (std::getline(d_file.handle(), line)) {
    
    if (line.substr(0, 14) == "##source=PLINK") {
      LOG(QGS::Log::WARNING)
        << "File `" << d_fname << "` was created by PLINK\n"
        << "PLINK implementation of VCF files is broken, attempting work-around\n"
        << "As a result, the options --hard-calls and --allow-missings are ignored\n"
        << "To remove this warning and disable the work-around, remove any lines "
        << "starting with `##source=PLINK` from VCF file.\n";
      d_format = "PLINK";
    }

    // find formats in file
    if (line.substr(0, 9) == "##FORMAT=") {
      found_format_tag = true;
      if (d_format != "PLINK" && !d_hard_calls && line.find("ID=DS") != std::string::npos)
        d_format = "DS"; // preferred format
      else if (d_format.empty() && line.find("ID=GT") != std::string::npos)
        d_format = "GT"; // least preferred format
      continue;
    }

    if (line.substr(0, 6) != "#CHROM")
      continue; // ignore everything up until header line

    std::istringstream iss(line);
    std::string token;
    while (iss >> token)
      d_header.push_back(std::move(token));
      
    break;
  }
  
  if (d_header.size() < 10) {
    LOG(QGS::Log::FATAL) << "No samples "
      << "found in file " << d_fname << ": Can't use input file.\n";
    std::exit(EXIT_FAILURE);
  }
  
  d_num_samples = d_header.size() - 9;
  
  std::copy(std::begin(d_header) + 9, std::end(d_header), std::back_inserter(d_sample));

  LOG(QGS::Log::VERBOSE) << "Opened file `" << d_fname << ". "
                             "Read mode: " << d_format << ". "
                             "Found " << d_num_samples << " subjects.\n";
  
  if (d_format.empty()) {
    if (!found_format_tag) {
      LOG(QGS::Log::FATAL) << "VCF file `" << d_fname << "` "
        << "does not contain a FORMAT tag in the header. Add one.\n";
    }
    LOG(QGS::Log::FATAL) << "No supported data format ("
      << (d_hard_calls ? "GT" : "GT, DS")
      << ") found in file " << d_fname << ": Can't use input file.\n";
    std::exit(EXIT_FAILURE);
  }
  
  return true;
}

void VCFreader::parse_line(VCFreader::Locus & l) {
  
  std::string line;
  while (std::getline(d_file.handle(), line)) {

    if (line.substr(0, 1) == "#")
      continue; // comment line, ignore

    d_buffer.ignore(std::numeric_limits<std::streamsize>::max());
    d_buffer.clear();
    d_buffer.str(line);
    
    l.switch_ar = false; // reset switch
    if (!(d_buffer >> l.chr_str >> l.pos >> l.id >> l.ref >> l.alt >> l.qual >> l.filter >> l.info_str >> l.format))
      continue;

    // parse chromosome
    if (l.chr_str == "X" || l.chr_str == "x")
      l.chr = 23;
    else if (l.chr_str == "Y" || l.chr_str == "y")
      l.chr = 24;
    else if (l.chr_str == "MT" || l.chr_str == "mt" || l.chr_str == "Mt" || l.chr_str == "mT")
      l.chr = 25;
    else {
    
      // parse "chr1" etc
      if (l.chr_str.size() > 3 && l.chr_str.substr(0,3) == "chr")
        l.chr_str = l.chr_str.substr(3);

      try {
        l.chr = std::stoul(l.chr_str);
        if (l.chr == 0 || l.chr > 26)
          throw "not a chromosome code";
      }
      catch (...) {
       continue;
      }
    }
      
    for (char & c : l.ref)
      c = std::toupper(c);

    for (char & c : l.alt)
      c = std::toupper(c);
      
    break;
  }
}

bool VCFreader::read_gt(VCFreader::Locus & l) {

  if (!d_buffer) {
    LOG(QGS::Log::VERBOSE) << "In file `" << d_fname << "` "
      << "locus " << l << " appears to be duplicated. "
      << "Skipping.\n";
    return false;
  }
 
  // find GT in format
  std::size_t const gt_pos = l.format.find("GT");
  if (gt_pos == std::string::npos)
    return false;
  std::size_t const gt_idx = std::count(l.format.begin(), l.format.begin() + gt_pos, ':');

  l.data_ds.clear();
  l.data_ds.resize(d_num_samples, 0);
  std::size_t idx = -1;
  char c, prev_allele = '!';

  std::size_t colon_count = 0, allele_count = 2;
  long double ds_sum = 0;
  std::size_t male_unilog = 0;
  while (d_buffer.get(c)) {

    if (c == ':') {
      ++colon_count;
      continue;
    }
    if (c == ' ' || c == '\t') {
      colon_count = 0;
      ++idx;
      if (allele_count == 1 && l.chr > 22) {
        ++male_unilog;
        ++idx;
      }
      else if (allele_count != 2) {
        LOG(QGS::Log::WARNING) << "In file `" << d_fname << "` "
          << "subject " << d_sample.at(idx / 2) << " has incomplete "
          << "data for locus " << l << ". Skipping.\n";
        return false;
      }
      allele_count = 0;
      continue;
    }

    if (colon_count != gt_idx) {
      continue;
    }

    if (c == '|' || c == '/') {
      ++idx;
      continue;
    }
    
    ++allele_count;
    prev_allele = c;
    
    if (c == '.') {

      LOG(QGS::Log::TRACE) << "vcf_reader: missing data point "
        << "for subject #" << idx / 2 << '\n';
      
      if (!d_allow_missings)
        return false; // exclude snp

      l.data_ds[idx/2] = NAN; // missing

      continue;
    }

    if (c < '0' || c > '9') {
      LOG(QGS::Log::WARNING) << "vcf_reader: unexpected character in file `" << d_fname
                                 << "`: `" << c << "`\n";
      continue;
    }

    if ( (prev_allele == '1' && !l.switch_ar) || (prev_allele == '0' && l.switch_ar) ) {
      ds_sum += 1;
      l.data_ds[idx/2] += 1;
    }

  }
  
  if (allele_count == 1 && l.chr > 22) {
    ++male_unilog;
    ++idx;
  }
  else if (allele_count != 2) {
    LOG(QGS::Log::WARNING) << "Skipping locus.\n";
    return false;
  }
  
  if (d_num_samples * 2 != idx + 1) {
    LOG(QGS::Log::WARNING) << "Read " << (idx + 1) / 2 << 
      " individuals, expected " << d_num_samples << ": skipping locus "
      << l << '\n';
    return false;
  }

  l.maf = ds_sum / (d_num_samples * 2.0 - male_unilog);
  if (l.maf > 0.5) {
    LOG(QGS::Log::TRACE) << "Flipped MAF FOR " << l << '\n';
    l.maf = 1 - l.maf;
  }

  return true;
}

bool VCFreader::read_ds(VCFreader::Locus & l) {

  if (!d_buffer) {
    LOG(QGS::Log::VERBOSE) << "In file `" << d_fname << "` "
      << "locus " << l << " appears to be duplicated. "
      << "Skipping.\n";
    return false;
  }

  // find DS in format
  std::size_t const gt_pos = l.format.find("DS");
  if (gt_pos == std::string::npos) {
    LOG(QGS::Log::WARNING) << "In file " << d_fname << " for " 
      << l.chr << ":" << l.pos << " no " << d_format << " info found";
    return false;
  }
  std::size_t const gt_idx = std::count(l.format.begin(), l.format.begin() + gt_pos, ':');

  l.data_ds.clear(); // not needed, strictly
  l.data_ds.resize(d_num_samples, 0);
  std::size_t idx = 0;
  char c;

  std::size_t colon_count = 0;
  long double ds_sum = 0;
  while (d_buffer.get(c)) {

    if (c == ':')
      ++colon_count;
      
    if (colon_count == gt_idx) {
      float ds;
      if (!(d_buffer >> ds)) {
        LOG(QGS::Log::WARNING) << "Failed to read " << l << "\n";
        LOG(QGS::Log::WARNING) << "Something is wrong with the VCF file. Skipping locus. \n";
        return false;
      }
      
      if (l.switch_ar)
        ds = 2 - ds;
        
      //LOGQGS::Log::WARNING << "read " << ds;
      l.data_ds[idx] = ds;
      ++idx;
      ds_sum += ds;
    }

    if (c == ' ' || c == '\t') {
      colon_count = 0;
    }

  }
  
  if (d_num_samples != idx) {
    LOG(QGS::Log::WARNING) << "Read " << idx << " individuals, expected " << d_num_samples << ": skipping locus\n";
    return false;
  }
   
  l.maf = ds_sum / (d_num_samples * 2.0);
  if (l.maf > 0.5) {
    LOG(QGS::Log::TRACE) << "Flipped MAF FOR " << l << '\n';
    l.maf = 1 - l.maf;
  }
  return true;
}

bool VCFreader::read_plink(VCFreader::Locus & l) {

  if (!d_buffer) {
    LOG(QGS::Log::VERBOSE) << "In file `" << d_fname << "` "
      << "locus " << l << " appears to be duplicated. "
      << "Skipping.\n";
    return false;
  }

  l.data_ds.clear(); // not needed, strictly
  l.data_ds.resize(d_num_samples, 0);
  
  // parse format
  std::vector<std::string> format;
  std::size_t p1 = 0, p2 = 0;
  for (; (p2 = l.format.find(":", p1)) != std::string::npos; p1 = p2 + 1)
    format.push_back(l.format.substr(p1, p2 - p1));
  format.push_back(l.format.substr(p1));

  std::string genostr;
  std::size_t sample_idx = 0;
  long double ds_sum = 0;
  while (d_buffer >> genostr) {
    
    // parse genostr
    std::unordered_map<std::string, std::string> genotype;
    std::size_t idx = 0;
    for (p1 = 0, p2 = 0; (p2 = genostr.find(":", p1)) != std::string::npos; p1 = p2 + 1)
      genotype[format.at(idx++)] = genostr.substr(p1, p2 - p1);
    genotype[format.at(idx++)] = genostr.substr(p1);
    
    // parse genotype
    bool correct_read = false;
    auto itt = genotype.find("DS");
    if (itt != genotype.end()) {
      // dosage available
      try {
        l.data_ds[sample_idx] = std::stof(itt->second);
        correct_read = true;
      }
      catch (...) {
        LOG(QGS::Log::DEBUG) << "Failed to read " << l << "\n";
        LOG(QGS::Log::DEBUG) << "Something is wrong with the VCF file (ds). Trying GT.\n";
      }
    }
    if (!correct_read && (itt = genotype.find("GT")) != genotype.end()) {
      // genotype available
      try {
        int const gt1 = itt->second.at(0) - '0';
        int const gt2 = itt->second.at(2) - '0';
        if (gt1 < 0 || gt1 > 1 || gt2 < 0 || gt2 > 1)
          throw "Parse error";
        l.data_ds[sample_idx] = gt1 + gt2;
        correct_read = true;
      }
      catch (...) {
        LOG(QGS::Log::WARNING) << "Failed to read " << l << "\n";
        LOG(QGS::Log::WARNING) << "Something is wrong with the VCF file (gt). Skipping locus. \n";
      }
    }
    
    if (!correct_read)
      return false;
      
    if (l.switch_ar)
      l.data_ds[sample_idx] = 2 - l.data_ds[sample_idx];
    ds_sum += l.data_ds[sample_idx];
    ++sample_idx;
  }
  
  if (d_num_samples != sample_idx) {
    LOG(QGS::Log::WARNING) << "Read " << sample_idx << " individuals, expected " << d_num_samples << ": skipping locus\n";
    return false;
  }
   
  l.maf = ds_sum / (d_num_samples * 2.0);
  if (l.maf > 0.5) {
    LOG(QGS::Log::TRACE) << "Flipped MAF (" << l.maf << ") for " << l << '\n';
    l.maf = 1 - l.maf;
  }
  
  return true;
  
}

