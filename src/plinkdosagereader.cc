#include "plinkdosagereader.h"
#include "naturalsort.h"
#include <sstream>
#include <ostream>
#include <iostream>
#include <numeric>
#include "log.h"


Plinkdosagereader::Plinkdosagereader(std::string const & fname)
 :
  SNPreader(fname),
  d_map_fname{find_map_file()},
  d_map_file(d_map_fname),
  d_force_hardcalls{false},
  d_linenr{0},
  d_filenr{0},
  d_max{-1}
{
  
  if (!d_file) {
		LOG(QGS::Log::FATAL) << "Cannot open input file `" << d_fname 
      << "`: for reading. Aborting.";
		std::exit(EXIT_FAILURE);
  }
  
  if (d_map_fname.empty() || !d_map_file) {
		LOG(QGS::Log::FATAL) << "Cannot find map file of `" << d_fname 
      << "`.Aborting.";
		std::exit(EXIT_FAILURE);
  }
  
  parse_header();
}

Plinkdosagereader::Plinkdosagereader(std::vector<std::string> const & files)
 :
  SNPreader(""),
  d_fnames(files),
  d_map_file(""),
  d_force_hardcalls{false},
  d_linenr{0},
  d_filenr{0},
  d_max{-1}
{
	
	SI::natural::sort(d_fnames);
	
	d_filenr = -1;
  if (!open_next()) {
		LOG(QGS::Log::FATAL) << "Cannot open input files: empty list."
      << "Aborting.";
		std::exit(EXIT_FAILURE);
  }
}

std::string Plinkdosagereader::find_map_file() const {
  
  std::string fbase = d_fname;

  LOG(QGS::Log::TRACE) << "Attempting to find plink map file...";
  
  std::vector<std::string> files = {fbase + ".map", fbase + ".map.gz"};
  
  for (std::size_t n = 0; n != 2; ++n) {
    std::size_t pos = fbase.find_last_of('.');
    if (pos != std::string::npos) {
      fbase = fbase.substr(0, pos);
      files.push_back(fbase + ".map");
      files.push_back(fbase + ".map.gz");
    }
  }

  // attempting to find map file
  for (auto const & file : files) {
    LOG(QGS::Log::VERBOSE) << "Probe map file `" << file  << "`.";
    std::ifstream tmp(file);
    if (tmp) {
      LOG(QGS::Log::VERBOSE) << "Map file found.";
      return file;
    }
  }
  
  LOG(QGS::Log::VERBOSE) << "Map file not found.";
  return "";

}

bool Plinkdosagereader::deep_read(SNPreader::Locus & l) {
	
	l.data_ds.clear();

  std::vector<float> data;
  data.reserve(d_num_samples * 3);

  float f, max = 0;
  while (d_buffer >> f) {
    max = f > max ? f : max;
    data.push_back(f);
  }
  
  if (d_max < 0) {
    d_max = max;
    LOG(QGS::Log::TRACE) << "First dosage line contains max dosage"
      " of " << d_max << ": we assume a 0-" <<
      (d_max > 1 ? 2 : 0) << " dosage scale.\n";
  }
  
  if (d_max <= 1 && max > 1) {
    LOG(QGS::Log::FATAL) << "Reading locus " << l.id << 
       " we discovered our intial guess of 0-1 dosages was incorrect. "
      "Please discard current results and rerun with dosage specified.";
    std::exit(EXIT_FAILURE);
  }
    
  std::size_t const num_read = data.size();
  
  if (num_read == d_num_samples) {
    // data represents dosages themselves: easy
    if (d_max <= 1)
      for (auto & ds : data)
        ds *= 2;
    if (l.switch_ar) { // SWITCH
      //LOG(QGS::Log::WARNING) << "Locus was switched, using 2 - dosage."; // SWITCH
      for (auto & ds : data) // SWITCH
        ds = 2 - ds; // SWITCH
    }
    std::swap(l.data_ds, data);
    long double const sum = std::accumulate(std::begin(l.data_ds), std::end(l.data_ds), 0.0L);
    l.maf = sum / (d_num_samples * 2);
    return true;
  }
  
  if (d_max > 1) {
    LOG(QGS::Log::FATAL) << "Dosage file `" << d_fname << 
       "` contains probabilities but has value " << d_max <<
       " > 1 on line " << d_linenr << ". Aborting run.";
    std::exit(EXIT_FAILURE);
  }
  
  l.data_ds.reserve(num_read);
  
  
  if (num_read == d_num_samples * 2) {

    // data represents prob. of A1/A1 and A1/A2
    for (std::size_t idx = 0; idx != num_read; idx += 2) {
      float a1 = data[idx], a2 = data[idx + 1], a3 = 1 - (data[idx] + data[idx+1]);
      if (a3 < 0) {
        LOG(QGS::Log::FATAL) << "Dosage file `" << d_fname << 
           "` contains probabilities " << data[idx] << " + " <<
           data[idx+1] << " > 1 on line " << d_linenr <<
           " > 1. Aborting run (1).";
        std::exit(EXIT_FAILURE);
      }
			if (l.switch_ar) { // SWITCH
        //LOG(QGS::Log::WARNING) << "Locus was switched, using swap(a1, a3)"; // SWITCH
				std::swap(a1, a3); // SWITCH
      }
      float const val = a2 + 2 * a3;
      if (val < 0 || val > 2) {
        LOG(QGS::Log::FATAL) << "Dosage file `" << d_fname << 
           "` contains probabilities " << data[idx] << " + " <<
           data[idx+1] << " > 1 on line " << d_linenr <<
           " > 1. Aborting run (2).";
        std::exit(EXIT_FAILURE);
      }
      l.data_ds.push_back(val);
    }
    
    if (l.data_ds.size() != d_num_samples) {
			LOG(QGS::Log::FATAL) << "Didn't read enough dosages.";
			std::exit(EXIT_FAILURE);
    }
    
    long double const sum = std::accumulate(std::begin(l.data_ds), std::end(l.data_ds), 0.0L);
    l.maf = sum / (d_num_samples * 2);

    return true;
  }
  
  
  if (num_read == d_num_samples * 3) {
    // data represents prob. of A1/A1, A1/A2, and A2/A2
    for (std::size_t idx = 0; idx != num_read; idx += 3) {
			float a1 = data[idx], a2 = data[idx + 1], a3 = data[idx+2];
			if (l.switch_ar) { // SWITCH
        //LOG(QGS::Log::WARNING) << "Locus was switched, swap(a1, a3) (2)"; // SWITCH
				std::swap(a1, a3); // SWITCH
      }
			float const val = a2 + 2 * a3;
      if (val < 0 || val > 2) {
        LOG(QGS::Log::FATAL) << "Dosage file `" << d_fname << 
           "` contains probabilities " << data[idx] << ", " <<
           data[idx+1] << ", " << data[idx+2] << " on line " << d_linenr
           << " > 1. Aborting run (3).";
        std::exit(EXIT_FAILURE);
      }
      l.data_ds.push_back(val);
    }
    
    long double const sum = std::accumulate(std::begin(l.data_ds), std::end(l.data_ds), 0.0L);
    l.maf = sum / (d_num_samples * 2);

    return true;
  }
  
	return false;
}

bool Plinkdosagereader::parse_header() {
  
  std::string line;
  if (!std::getline(d_file.handle(), line)) {
    return false;
  }

  std::istringstream iss(line);
  std::string snp, a1, a2;
  if (!(iss >> snp >> a1 >> a2)) {
    return false;
  }
  
  if (snp != "SNP" || a1 != "A1" || a2 != "A2") {
    return false;
  }

  // everything is in order, read ids
  
  std::string fid, iid;
  std::vector<std::string> sample;
  sample.reserve(d_num_samples);
  while (iss >> fid) {
    if (!(iss >> iid)) {
      LOG(QGS::Log::WARNING) << "Found sample fid without iid in "
        "file `" << d_fname << "` on line 1. fid=" << fid << 
        ": ignoring individual.";
      break;
    }
    sample.push_back(fid + "_" + iid);
  }
  
  if (!d_num_samples) { // first (only?) file
		d_sample = std::move(sample);
		d_num_samples = d_sample.size();
		if (!d_num_samples) {
      LOG(QGS::Log::FATAL) << "File `" << d_fname << "` does "
        "not have any samples. Aborting.";
      std::exit(EXIT_FAILURE);
    }
	}
  else { // not the first file
		if (sample != d_sample) {
      LOG(QGS::Log::FATAL) << "File `" << d_fname << "` has "
        "different subjects than previous file: can't proceed. Aborting.";
      std::exit(EXIT_FAILURE);
	  }
  }

  LOG(QGS::Log::VERBOSE) << "Opened file `" << d_fname << ". "
                             "Read mode: plink dosage. "
                             "Found " << d_num_samples << " subjects.";
                             
  return true;
}

void Plinkdosagereader::parse_line(SNPreader::Locus & l) {
  
  // here we should do three things
  // 1. read the map file and get the rs-nr & chr & pos
  // 2. read the dosage file line and get rs-nr & a1 & a2
  // 3. make sure rs-nrs match

  ++d_linenr;

  std::string line;
  if (!std::getline(d_map_file.handle(), line)) {
		if (d_map_file.handle().eof()) { // we've reached EOF
			LOG(QGS::Log::TRACE) << "Map EOF: opening next";
			if (open_next())
				return parse_line(l);
		}
    LOG(QGS::Log::TRACE) << "Can't read line from map file.";
    d_file.handle().setstate(std::ios_base::failbit);
    return;
  }
  
  // 9 rs573167194 0 141000084
  std::istringstream mapline(line);
  
  auto old_locus = l;
  
  if (!(mapline >> l.chr >> l.id >> l.info_str >> l.pos)) {
    LOG(QGS::Log::TRACE) << "Can't parse line from map file.";
    d_file.handle().setstate(std::ios_base::failbit);
    return;
  }
  
  if (old_locus.chr > l.chr || (old_locus.chr == l.chr && old_locus.pos > l.pos)) {
		LOG(QGS::Log::FATAL) << "File `" << d_fname << "` line " << d_linenr
		  << " has wrong locus order.\nPrev. locus: " << old_locus
		  << "\nCurrent locus: " << l << "\nAborting.";
		std::exit(EXIT_FAILURE);
  }
    
  // test if empty
  std::string tmp;
  if (mapline >> tmp) {
    LOG(QGS::Log::WARNING) << "Unexpected data in map file `" 
      << d_map_fname << "`. Value=" << tmp << ". Ignoring.";
  }
  
  if (!std::getline(d_file.handle(), line)) {
    LOG(QGS::Log::TRACE) << "Can't read line from dosage file.";
    return;
  }
  
  // rs573167194 T A 1 0......
  d_buffer = std::istringstream(line);
  
  std::string id;
  if (!(d_buffer >> id >> l.ref >> l.alt)) {
    LOG(QGS::Log::TRACE) << "Can't parse line from dosage file.";
    d_file.handle().setstate(std::ios_base::failbit);
    return;
  }
  
  if (id != l.id) {
    LOG(QGS::Log::WARNING) << "Dosage and map file out of sync on "
      "map line " << d_linenr << ": read snps " << id << " (dosage) "
      "and " << l.id << " (map). Stopping read.";
    l.chr = 0;
    l.pos = 0;
  }

}

bool Plinkdosagereader::open_next()
{

	++d_filenr;

  if (d_filenr >= d_fnames.size())
    return false;

  d_linenr = 0;

	d_fname = d_fnames[d_filenr];
	d_file = GZfile(d_fname);
	
  if (!d_file) {
		LOG(QGS::Log::FATAL) << "Cannot open input file `" << d_fname 
      << "`: for reading. Aborting.";
		std::exit(EXIT_FAILURE);
  }
	
  d_map_fname = find_map_file();
  d_map_file = GZfile(d_map_fname);
  
  if (d_map_fname.empty() || !d_map_file) {
		LOG(QGS::Log::FATAL) << "Cannot find map file of `" << d_fname 
      << "`.Aborting.";
		std::exit(EXIT_FAILURE);
  }
  
  if (!parse_header()) {
		// assuming this file is empty
		LOG(QGS::Log::WARNING) << "Input file `" << d_fname 
      << "` does not have proper header: skipping file.";
    return open_next();
  }
  
  return true;
}
