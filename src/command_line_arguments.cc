#include "command_line_arguments.h"
#include "gzfile.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include "cmd_line_options.h"

namespace QGS {

std::unordered_map<std::string, std::vector<std::string>>  command_line_arguments(int argc, char ** argv) {
  
  auto options = QGS::parse_cmd_line_options(argc, argv);
  
  struct Option {
    std::string desc;
    bool required;
    std::string type;
    std::string defaul;
  };
  
  std::unordered_map<std::string, Option> qgs_arguments;
  
  // misc arguments
  qgs_arguments["help"] = {"print this message", 0, ""};
  qgs_arguments["version"] = {"print version information", 0, ""};
  
  // required arguments
  qgs_arguments["sample"] = {"path to input file(s) in VCF or Plink dosage format. may be gz-compressed.", 1, "multi"};
  qgs_arguments["reference"] = {"path to reference file in VCF format. may be gz-compressed.", 1, "string"};
  qgs_arguments["genes"] = {"path to gene database in Gene Transfer Format. my be gz-compressed.", 1, "string"};
  qgs_arguments["out"] = {"filename for writeable output file. any existing file will be overwritten without warning.", 1, "string"};

  // defaulted arguments
  qgs_arguments["flank"] = {"symmetrical flanking region in kb", 0, "uint", "0"};
  qgs_arguments["pre-flank"] = {"asymmetrical pre-region flanking region in kb", 0, "uint"};
  qgs_arguments["post-flank"] = {"asymmetrical post-region flanking region in kb", 0, "uint"};
  qgs_arguments["maf"] = {"minor allele frequency threshold", 0, "double", "0.01"};
  
  // filter options
  qgs_arguments["gtf-filter"] = {"key-value pair(s) for filtering the gene file, e.g. type=gene gene_type=protein_coding", 0, "multi"};
  qgs_arguments["include-snps"] = {"file name to whitespace-separated snp-names for snps to include. snps not in the file will be excluded.", 0, "string"};
  qgs_arguments["exclude-snps"] = {"file name to whitespace-separated snp-names for snps to exclude. snps in the file will be excluded.", 0, "string"};
  qgs_arguments["chr"] = {"Chromosome to include, numeric value. X = 23, Y = 24, MT = 25", 0, "uint"};
  
  // print output arguments
  qgs_arguments["verbose"] = {"give more verbose output", 0};
  qgs_arguments["debug"] = {"give a lot more verbose output", 0};
  qgs_arguments["trace"] = {"give all possible output (warning: a lot)", 0};
  
  // misc
  qgs_arguments["hard-calls"] = {"force the program to use hard calls instead of DS data", 0};
  qgs_arguments["allow-missings"] = {"include loci with missing data points, printing 'nan' for those subjects/genes", 0};
  qgs_arguments["fill-missings"] = {"include loci with missing data points, assuming 0/0 for any missing variant", 0};
  qgs_arguments["weight-by"] = {"Weighted QGS. Provide VCF INFO field name which contains the weight for each locus, e.g. --weight-by R2", 0, "string"};
  qgs_arguments["delimiter"] = {"replace standard comma output delimiter, must be single character", 0, "string"};
  qgs_arguments["output-variants"] = {"include variants used for the QGS in output file", 0};
  
  // version information
  if (options.find("version") != options.end()) {
    std::exit(EXIT_SUCCESS);
  }
  
  // find longest option name
  std::size_t longest_option = 0;
  for (auto itt = qgs_arguments.cbegin(); itt != qgs_arguments.cend(); ++itt)
    longest_option = longest_option > itt->first.size() ? longest_option : itt->first.size();
  
  // print help
  if (options.find("help") != options.end()) {
    for (auto itt = qgs_arguments.cbegin(); itt != qgs_arguments.cend(); ++itt) {
      std::cout << "--" << std::left << std::setw(longest_option + 4) << itt->first;
      std::cout << (itt->second.required ? "[REQ]" : "[OPT]");
      std::cout << " " << itt->second.desc;
      if (!itt->second.defaul.empty())
        std::cout << " (default value: " << itt->second.defaul << ")";
      std::cout << '\n';
    }
    std::exit(EXIT_SUCCESS);
  }
  
  // add all defaulted arguments
  for (auto itt = qgs_arguments.cbegin(); itt != qgs_arguments.cend(); ++itt) {
    if (itt->second.defaul.empty())
      continue;
    if (options.find(itt->first) == options.end())
      options[itt->first].push_back(itt->second.defaul);
  }
  
  // test whether all required arguments are set
  bool everything_ok = true;
  for (auto itt = qgs_arguments.cbegin(); itt != qgs_arguments.cend(); ++itt) {
    if (!itt->second.required)
      continue;
    if (options.find(itt->first) == options.end()) {
      LOG(QGS::Log::FATAL) << "Required flag `--" << itt->first
        << "` is missing\n";
      everything_ok = false;
    }
  }
  
  // test whether arguments have correct types
  for (auto itt = qgs_arguments.cbegin(); itt != qgs_arguments.cend(); ++itt) {
    if (options.find(itt->first) == options.end())
      continue;

    if (itt->second.type.empty() && !options[itt->first].empty()) {
      everything_ok = false;
      LOG(QGS::Log::FATAL) << "Flag `--" << itt->first
        << "` should be provided without any value\n";
    }
    else if (itt->second.type == "multi" && options[itt->first].empty()) {
      everything_ok = false;
      LOG(QGS::Log::FATAL) << "Flag `--" << itt->first
        << "` should be provided with 1 or more values\n";
    }
    else if ((itt->second.type == "string" || itt->second.type == "uint" || itt->second.type == "double") && options[itt->first].size() != 1) {
      everything_ok = false;
      LOG(QGS::Log::FATAL) << "Flag `--" << itt->first
        << "` should be provided with 1 value\n";
    }
    else if (itt->second.type == "uint") {
      try {
        long long tmp = stoll(options[itt->first][0]);
        if (tmp < 0 || std::to_string(tmp) != options[itt->first][0])
          throw 1;
      }
      catch (...) {
        everything_ok = false;
        LOG(QGS::Log::FATAL) << "Flag `--" << itt->first
          << "` should be provided with a non-negative integer numeric value\n";
      }
    }
    else if (itt->second.type == "double") {
      try {
        stod(options[itt->first][0]);
      }
      catch (...) {
        everything_ok = false;
        LOG(QGS::Log::FATAL) << "Flag `--" << itt->first
          << "` should be provided with a numeric value\n";
      }
    }
  }

  if (!everything_ok)
    std::exit(EXIT_FAILURE);
    
  // test whether all arguments are actually allowed
  for (auto itt = options.cbegin(); itt != options.cend(); ++itt) {
    if (qgs_arguments.find(itt->first) == qgs_arguments.end()) {
      LOG(QGS::Log::WARNING) << "Unknown flag `--" << itt->first
        << "` is ignored\n";
    }
  }
    
  // print all arguments in effect
  std::cout << "Command line arguments:\n";
  for (auto itt = options.cbegin(); itt != options.cend(); ++itt) {
    std::cout << "--" << std::left << std::setw(longest_option + 4) << itt->first;
    for (std::size_t idx = 0; idx != itt->second.size(); ++idx)
      std::cout << std::string((longest_option + 6) * (idx != 0), ' ') + itt->second[idx] + std::string("\n");
    if (itt->second.empty())
      std::cout << '\n';
  }
  std::cout << '\n';
  

  return options;
}

std::unordered_map<std::string, std::string> create_gtf_filter(std::unordered_map<std::string, std::vector<std::string>> const & args) {
  std::unordered_map<std::string, std::string> filter;
  
  std::unordered_map<std::string, std::vector<std::string>>::const_iterator itt = args.find("gtf-filter");
  if (itt != args.end()) {
    std::vector<std::string> const & filter_vec = itt->second;
    
    for (std::size_t idx = 0; idx != filter_vec.size(); ++idx) {
      std::string const & str = filter_vec[idx];
      std::size_t const eq_idx = str.find('=');
      if (eq_idx == std::string::npos)
        continue;
      filter[str.substr(0, eq_idx)] = str.substr(eq_idx + 1);
    }
  }

  return filter;
}

void read_snp_file(std::string const & fname, std::unordered_map<std::string, bool> & out, bool value) {
  GZfile file(fname);
  if (!file)
    LOG(QGS::Log::WARNING) << "Cannot open file `" << fname << "`. No loci read.\n";
  std::string loc;
  while (file.handle() >> loc)
    out[std::move(loc)] = value;
  if (out.empty())
    LOG(QGS::Log::WARNING) << "No loci read from file `" << fname << "`.\n";
  LOG(QGS::Log::VERBOSE) << "Read " << out.size() << " loci from `" << fname << "`.\n";
}

std::unordered_map<std::string, bool> create_snp_filter(std::unordered_map<std::string, std::vector<std::string>> const & args) {
  std::unordered_map<std::string, bool> out;
  
  if (args.find("include-snps") != args.end() && args.find("exclude-snps") != args.end()) {
    LOG(QGS::Log::WARNING) << "Both --include-snps and --exclude-snps have been set. "
      << "This is not possible. --exclude-snps will be ignored.\n";
  }

  std::unordered_map<std::string, std::vector<std::string>>::const_iterator itt;
  if ((itt = args.find("include-snps")) != args.end())
    read_snp_file(itt->second[0], out, true);
  else if ((itt = args.find("exclude-snps")) != args.end())
    read_snp_file(itt->second[0], out, false);
  return out;
}

}

