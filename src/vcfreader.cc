#include "vcfreader.h"
#include <sstream>
#include <ostream>
#include <iostream>
#include "log.h"


VCFreader::VCFreader(std::string const & fname, bool force_hardcalls)
 :
  SNPreader(fname),
  d_force_hardcalls{force_hardcalls}
{
  parse_header();
}

bool VCFreader::deep_read(VCFreader::Locus & l) {
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
  while (std::getline(d_file.handle(), line)) {
		
		// find formats in file
		if (line.substr(0, 9) == "##FORMAT=") {
			if (!d_force_hardcalls && line.find("ID=DS") != std::string::npos)
			  d_format = "DS"; // preferred format
//			else if (d_format != "DS" && line.find("ID=GP") != std::string::npos)
//			  d_format = "GP";
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
		LOG(QGS::Log::FATAL) << "No supported data format ("
		  << (d_force_hardcalls ? "GT" : "GT, DS")
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
      
    d_buffer = std::istringstream(line);
    
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
	
	//std::cerr << d_buffer.str() << "\n\n\n";

  if (!d_buffer && !l.data_ds.empty()) {
    LOG(QGS::Log::DEBUG) << "Request to read same position twice: duplicate snp position in sample file?\n";
    // locus already read, may occor if sample file has multiple snps at one position
    // TODO currently first snp result is discarded because of storing map<chr, map<pos, Score>>
    // think about different solution?
    return true;
  }
  
  //std::cerr << "Doing deep read of " << d_fname << " locus " << l.pos << '\n';
      
  // find GT in format
  std::size_t const gt_pos = l.format.find("GT");
  if (gt_pos == std::string::npos)
    return false;
  std::size_t const gt_idx = std::count(l.format.begin(), l.format.begin() + gt_pos, ':');

  l.data_ds.clear();
  l.data_ds.resize(d_num_samples, 0);
  std::size_t idx = -1;
  char c, prev_allele = '!';
  
  //std::cerr << "Expecting " << d_num_samples << " samples\n";

  std::size_t colon_count = 0, allele_count = 2;
  long double ds_sum = 0;
  std::size_t male_unilog = 0;
  while (d_buffer.get(c)) {
		
		//std::cerr << '\n' << ++c_count << ":\t`" << c << "`\t" << idx;

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
				  << "data for locus " << l;
				LOG(QGS::Log::WARNING) << "Skipping locus.";
				break;
			}
			allele_count = 0;
      //std::cerr << " => " << idx;
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
      // missing data point:
      LOG(QGS::Log::TRACE) << "vcf_reader: missing data point "
        << "for subject " << idx / 2 << '\n';
      
      // for MAF, count as zero
			//ds_sum += 0;

			l.data_ds[idx/2] = -99; // missing, might change to -98 later

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
		LOG(QGS::Log::WARNING) << "In file `" << d_fname << "` "
			<< "subject " << d_sample.at(idx / 2) << " has incomplete "
			<< "data for locus " << l << ":" << ". Skipping locus.\n";
		return false;
	}
  
  if (d_num_samples * 2 != idx + 1) {
    LOG(QGS::Log::WARNING) << "Read " << (idx + 1) / 2 << 
      " individuals, expected " << d_num_samples << ": skipping locus "
      << l << '\n';
    return false;
  }

  l.maf = ds_sum / (d_num_samples * 2.0 - male_unilog);

  return true;
}

bool VCFreader::read_ds(VCFreader::Locus & l) {
	
  //LOGQGS::Log::WARNING << "reading " << d_fname << " => " << l.chr << ":" << l.pos;

  if (!d_buffer && !l.data_ds.empty()) {
    LOG(QGS::Log::DEBUG) << "Request to read same position twice: duplicate snp position in sample file?\n";
    // locus already read, may occor if sample file has multiple snps at one position
    // TODO currently first snp result is discarded because of storing map<chr, map<pos, Score>>
    // think about different solution?
    return true;
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

  std::unordered_map<char, std::size_t> maf_count;
  std::size_t colon_count = 0;
  long double ds_sum = 0;
  while (d_buffer.get(c)) {

    if (c == ':')
      ++colon_count;
      
    if (colon_count == gt_idx) {
      float ds;
      if (!(d_buffer >> ds)) {
        LOG(QGS::Log::WARNING) << "failed to read ds float\n";
        std::exit(EXIT_FAILURE);
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

  return true;
}

