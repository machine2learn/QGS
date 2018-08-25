#include "plinkbedreader.h"
#include <sstream>
#include <ostream>
#include <iostream>
#include <numeric>
#include <bitset>
#include "log.h"


Plinkbedreader::Plinkbedreader(std::string const & fname)
 :
  SNPreader(fname),
  d_bim_fname{find_file("bim")},
  d_fam_fname{find_file("fam")},
  d_bim_file(d_bim_fname),
  d_fam_file(d_fam_fname),
  d_num_pos_read{0},
  d_bytes_per_locus{0}
{
  
  if (!d_file) {
    LOG(QGS::Log::FATAL) << "Cannot open input file `" << d_fname 
      << "`: for reading. Aborting.\n";
    std::exit(EXIT_FAILURE);
  }
  
  if (d_bim_fname.empty() || !d_bim_file) {
    LOG(QGS::Log::FATAL) << "Cannot find bim file of `" << d_fname 
      << "`. Aborting.\n";
    std::exit(EXIT_FAILURE);
  }
  
  if (d_fam_fname.empty() || !d_fam_file) {
    LOG(QGS::Log::FATAL) << "Cannot find fam file of `" << d_fname 
      << "`. Aborting.\n";
    std::exit(EXIT_FAILURE);
  }
  
  unsigned char const magic_bytes[] = {0x6c, 0x1b, 0x01};
  for (auto const magic_byte : magic_bytes) {
    unsigned char const byte = d_file.handle().get();
    if (byte != magic_byte) {
      LOG(QGS::Log::FATAL) << "Input file `" << d_fname << "not "
        "recognised as PLINK BED format, or in id-major format. Maybe "
        "your PLINK version is out of date (or too new?). Aborting.\n";
      std::exit(EXIT_FAILURE);
    }
  }
  
  if (!parse_header()) {
    LOG(QGS::Log::FATAL) << "Unable to parse `" << d_fname << ". Aborting.\n";
    std::exit(EXIT_FAILURE);
  }
  
  d_bytes_per_locus = d_num_samples / 4;
  if (d_num_samples % 4 != 0)
    ++d_bytes_per_locus;
  
  LOG(QGS::Log::VERBOSE) << "Opened file `" << d_fname << ". "
   "Read mode: plink bed. "
   "Found " << d_num_samples << " subjects.\n";
}

std::string Plinkbedreader::find_file(std::string const & ext) const {
  
  std::string fbase = d_fname;

  LOG(QGS::Log::TRACE) << "Attempting to find plink map file...\n";
  
  std::vector<std::string> files = {fbase + "." + ext, fbase + "." + ext + ".gz"};
  
  for (std::size_t n = 0; n != 2; ++n) {
    std::size_t pos = fbase.find_last_of('.');
    if (pos != std::string::npos) {
      fbase = fbase.substr(0, pos);
      files.push_back(fbase + "." + ext);
      files.push_back(fbase + "." + ext + ".gz");
    }
  }

  // attempting to find file
  for (auto const & file : files) {
    LOG(QGS::Log::TRACE) << "Probe file `" << file  << "`.\n";
    std::ifstream tmp(file);
    if (tmp) {
      LOG(QGS::Log::TRACE) << "File found.\n";
      return file;
    }
  }
  
  LOG(QGS::Log::TRACE) << "File not found.\n";
  return "";

}

bool Plinkbedreader::deep_read(SNPreader::Locus & l) {
  l.data_ds.clear();
  l.data_ds.resize(d_num_samples);
  
  if (d_num_pos_read == 0) {
    LOG(QGS::Log::FATAL) << "Plinkbedreader::deep_read: should never happen.\n";
    std::exit(EXIT_FAILURE);
  }

  // remove skipped loci from read buffer
  if (d_num_pos_read > 1)
    d_file.handle().ignore(d_bytes_per_locus * (d_num_pos_read - 1));
  
  d_num_pos_read = 0;
  
  // do actual reading
  std::vector<char> buffer;
  buffer.resize(d_bytes_per_locus);
  
  d_file.handle().read(&buffer[0], d_bytes_per_locus);
  std::size_t gcount = d_file.handle().gcount();
  
  if (gcount != d_bytes_per_locus) {
    LOG(QGS::Log::WARNING) << "Plinkbedreader::deep_read: read error.\n";
    return false;
  }
  
  std::size_t const dosages[] = {l.switch_ar ? 2U : 0U, 1U, l.switch_ar ? 0U : 2U};
  
  std::size_t sample_idx = 0;
  long double total_dosage = 0;
  for (auto byte : buffer) {
    for (std::size_t offset = 0; offset != 4; ++offset) {
      unsigned char mask = 0x03 << (offset * 2);
      char const val = (mask & byte) >> (offset * 2);
      switch (val) {
        case 0x00: l.data_ds[sample_idx] = dosages[0]; total_dosage += dosages[0]; break;
        case 0x01: l.data_ds[sample_idx] = -999; break;
        case 0x02: l.data_ds[sample_idx] = dosages[1]; total_dosage += dosages[1]; break;
        case 0x03: l.data_ds[sample_idx] = dosages[2]; total_dosage += dosages[2]; break;
        default:
          LOG(QGS::Log::FATAL) << "Plinkbedreader::deep_read: Error in "
            "PLINK BED format. Check file `" << d_fname << "`\n";
          std::exit(EXIT_FAILURE);
      }
      std::cout << std::bitset<8>(byte) << ": " << std::bitset<8>(mask) << ": " << std::bitset<8>(val) << '\n';
      if (++sample_idx == d_num_samples)
        break;
    }
  }
  
  l.maf = total_dosage / (d_num_samples * 2);
  
  return true;
}

bool Plinkbedreader::parse_header() {
  
  // parse the fam file
  std::string line, fid, iid, f_iid, m_iid, sex, pheno;
  std::size_t fam_line_nr = 0;
  while (std::getline(d_fam_file.handle(), line)) {
    ++fam_line_nr;
    std::istringstream iss(line);
    if (!(iss >> fid >> iid >> f_iid >> m_iid >> sex >> pheno)) {
      LOG(QGS::Log::FATAL) << "PLINK fam file `" << d_fam_fname << 
        " read error on line " << fam_line_nr << ": cannot parse `" <<
        line << "` as subject. Aborting.\n";
      return false;
    }
    d_sample.push_back(fid + "_" + iid);
  }
  
  d_num_samples = d_sample.size();
  return true;
}

void Plinkbedreader::parse_line(SNPreader::Locus & l) {
  
  l.chr = 0;
  l.switch_ar = false;
  
  // read next bim entry
  std::string line;
  while (std::getline(d_bim_file.handle(), line)) {
    
    ++d_num_pos_read;
    std::istringstream iss(line);

    std::string cm_pos;
    if (!(iss >> l.chr_str >> l.id >> cm_pos >> l.pos >> l.ref >> l.alt))
      continue;

    // parse chromosome
    if (l.chr_str == "X" || l.chr_str == "x")
      l.chr = 23;
    else if (l.chr_str == "Y" || l.chr_str == "y")
      l.chr = 24;
    else if (l.chr_str == "MT" || l.chr_str == "mt" || l.chr_str == "Mt" || l.chr_str == "mT")
      l.chr = 25;
    else {
      try {
        l.chr = std::stoul(l.chr_str);
        if (l.chr == 0 || l.chr > 26)
          throw "not a chromosome code";
      }
      catch (...) {
       continue;
      }
    }
    
    break;
  }
  
  if (!d_bim_file)
    d_file.handle().setstate(std::ios_base::failbit);
}
