#include "qgs.h"

int main(int argc, char ** argv) {
  
  std::cout << "QGS " << QGS_VERSION << "\n\n";
  
  // ##### SET COMMAND LINE CONSTANTS
  
  std::unordered_map<std::string, std::vector<std::string>> const args = QGS::command_line_arguments(argc, argv);
  
  // required arguments
  std::vector<std::string> const & fname_sample = args.find("sample")->second;
  std::vector<std::string> const & fname_ref = args.find("reference")->second;
  std::string const & fname_genes = args.find("genes")->second.at(0);
  std::string const & fname_out = args.find("out")->second.at(0);
  
  // default arguments
  std::size_t const flank = std::stoull(args.find("flank")->second.at(0));
  std::size_t const pre_flank = args.find("pre-flank") != args.end() ? std::stoull(args.find("pre-flank")->second.at(0)) : flank;
  std::size_t const post_flank = args.find("post-flank") != args.end() ? std::stoull(args.find("post-flank")->second.at(0)) : flank;
  double const maf_limit = std::stod(args.find("maf")->second.at(0));
  
  // filter arguments
  std::unordered_map<std::string, std::string> const gtf_filter = QGS::create_gtf_filter(args);
  std::unordered_map<std::string, bool> const snp_filter = QGS::create_snp_filter(args);
  bool const snp_include_filter = args.find("include-snps") != args.end();
  long long const chr_filter = args.find("chr") != args.end() ? std::stoll(args.find("chr")->second.at(0)) : 0;

  // optional arguments
   bool const hard_calls = args.find("hard-calls") != args.end();
   bool const allow_missings = args.find("allow-missings") != args.end();
   char const delimiter = args.find("delimiter") == args.end() ? ',' : args.find("delimiter")->second.at(0)[0];
   
   // verbosity
   if (args.find("trace") != args.end())
     LOGLVL(QGS::Log::TRACE);
   else if (args.find("debug") != args.end())
     LOGLVL(QGS::Log::DEBUG);
   else if (args.find("verbose") != args.end())
     LOGLVL(QGS::Log::VERBOSE);

  // ##### OPEN INPUT AND OUTPUT FILES
  
  std::ofstream out_file(fname_out);
  if (!out_file) {
    LOG(QGS::Log::FATAL) << "Cannot open output file `" 
      << fname_out << "` for writing.\n";
    return 1;
  }

  GZfile gene_file(fname_genes);
  if (!gene_file)  {
    LOG(QGS::Log::FATAL) << "Cannot open gene file `" 
      << fname_genes << "` for reading.\n";
    return 1;
  }
  
  // open genetics files
  std::unique_ptr<SNPreader> reference_file = open_genetics_file(fname_ref, hard_calls);
  std::unique_ptr<SNPreader> sample_file = open_genetics_file(fname_sample, hard_calls, allow_missings);

  // ##### START OUTPUT

  LOG(QGS::Log::VERBOSE) << "Outputting subject ids and header.\n";
  out_file << "gene_name" << delimiter << "gene_id" << delimiter << "chr" 
           << delimiter << "start" << delimiter << "stop" << delimiter 
           << "Nsample" << delimiter << "Nref" << delimiter << "num_loci" << delimiter << "total_num_loci";
  for (std::size_t sidx = 0; sidx != sample_file->num_samples(); ++sidx)
    out_file << delimiter << sample_file->sample_id(sidx);
  out_file << '\n';
  
  // ##### MAIN READ BLOCK
  
  Genblock gb;

  SNPreader::Locus sample_locus, reference_locus;
  
  std::map<std::size_t, std::map<std::size_t, std::map<std::string, std::vector<QGS::Gene_score>>>> scores;

  while (*sample_file && *reference_file && gene_file >> gb) {

    LOG(QGS::Log::TRACE) << "Gene: read " << gb << " from gene file.\n";
    
    if (chr_filter && chr_filter > gb.chr) {
      LOG(QGS::Log::TRACE) << "Gene: gene fails chr filter, skipping.\n";
      continue;
    }
    if (chr_filter && chr_filter < gb.chr)
      break;
    
    if (!pass_gene_filter(gb, gtf_filter))
      continue;

    auto const strand_itt = gb.attr.find("strand");
    if (strand_itt == gb.attr.end() || strand_itt->second != "-") {
      gb.start = gb.start < pre_flank ? 0 : gb.start - pre_flank;
      gb.stop += post_flank;
    }
    else { // flip pre-flank and post-flank for reverse-strand genes
      gb.start = gb.start < post_flank ? 0 : gb.start - post_flank;
      gb.stop += pre_flank;
    }
    
    for (;;) {
      
      if (reference_locus.chr && gb < reference_locus) {
        LOG(QGS::Log::TRACE) << "Gene: gene lies in front of current loci. Read next gene.\n";
        break;
      }
      
      if (!reference_locus.chr || reference_locus <= sample_locus) {
        if (!(*reference_file >> reference_locus))
          break;
        LOG(QGS::Log::TRACE) << "Reference: read " << reference_locus << " from reference file.\n";
        scores[reference_locus.chr][reference_locus.pos][reference_locus.ref] = {};
      }
      else {
        if (!(*sample_file >> sample_locus)) {
          // sample file empty, continue reading to get #snps in last gene right
          sample_locus.chr = 99;
          continue;
        }
        LOG(QGS::Log::TRACE) << "Sample: read " << sample_locus << " from sample file.\n";
      }
      
      if (sample_locus.chr != reference_locus.chr || sample_locus.pos != reference_locus.pos)
        continue; // no hit
        
      if (!snp_filter.empty()) {
        // there is an include or exclude filter
        
        auto itt = snp_filter.find(sample_locus.id);
        if (itt == snp_filter.end())
          itt = snp_filter.find(std::to_string(sample_locus.chr) + ":" + std::to_string(sample_locus.pos));
          
        if (snp_include_filter && itt == snp_filter.end()) {
          LOG(QGS::Log::TRACE) << "Filter: Sample locus not in include filter. Skipping.\n";
          continue;
        }

        if (!snp_include_filter && itt != snp_filter.end()) {
          LOG(QGS::Log::TRACE) << "Filter: Sample locus in exclude filter. Skipping.\n";
          continue;
        }
        
        LOG(QGS::Log::TRACE) << "Filter: Sample locus passed filter. Including.\n";
      }
        
      if (sample_locus.chr != gb.chr || sample_locus.pos < gb.start || sample_locus.pos > gb.stop) {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " does not lie within gene. skipping.\n";
        continue;
      }
        
      LOG(QGS::Log::TRACE) << "Match: locus " 
        << sample_locus.chr << ":" << sample_locus.pos 
        << " found in sample and reference. Attempting deep read.\n";
        
      sample_locus.parse_alt();
      reference_locus.parse_alt();
      
      if (sample_locus.palt.empty() || reference_locus.palt.empty()) {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " is missing alt data: skipping locus.\n";
        continue;
      }
      
      if (sample_locus.ref == reference_locus.palt[0] && reference_locus.ref == sample_locus.palt[0]) {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " has flipped ref/alt data. Correcting.\n";
          sample_locus.switch_alt_ref();
      }
      else if (sample_locus.ref != reference_locus.ref || reference_locus.palt[0] != sample_locus.palt[0]) {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " has ref/alt mismatch between sample and reference. skipping locus.\n";
          sample_locus.clear();
          continue;
      }

      if (!sample_file->deep_read(sample_locus) || sample_locus.maf < maf_limit) {
        LOG(QGS::Log::TRACE) << "Match: sample locus excluded. skipping locus.\n";
        sample_locus.clear();
        continue;
      }
      
      if (!reference_file->deep_read(reference_locus) || reference_locus.maf < maf_limit) {
        LOG(QGS::Log::TRACE) << "Match: reference locus excluded. skipping locus.\n";
        sample_locus.clear();
        reference_locus.clear();
        continue;
      }
      
      LOG(QGS::Log::DEBUG) << "QGS: locus " 
        << sample_locus.chr << ":" << sample_locus.pos 
        << " deep read successful. Calculating QGS.\n";
      
      scores[sample_locus.chr][sample_locus.pos][sample_locus.ref] = QGS::score(sample_locus, reference_locus);
    }
    
    std::map<std::size_t, std::map<std::string, std::vector<QGS::Gene_score>>>::iterator lower_bound = scores[gb.chr].lower_bound(gb.start);
    scores[gb.chr].erase(scores[gb.chr].begin(), lower_bound);
    
    std::vector<QGS::Gene_score> total_score;
    total_score.resize(sample_file->num_samples(), 0);
    std::size_t snp_cnt = 0, total_snp_cnt = 0;
    for (auto itt_pos : scores[gb.chr]) {
      if (itt_pos.first > gb.stop)
        break;
      for (auto itt_var : itt_pos.second) {
        ++total_snp_cnt;
        if (itt_var.second.empty())
          continue;
        ++snp_cnt;
        for (std::size_t sample_idx = 0; sample_idx != total_score.size(); ++sample_idx) {
          if (itt_var.second[sample_idx] < 0 || total_score[sample_idx] < 0)
            total_score[sample_idx] = -99;
          else
            total_score[sample_idx] += itt_var.second[sample_idx];
        }
      }
    }
    
    if (!snp_cnt) {
      LOG(QGS::Log::TRACE) << "Gene: No loci in " << gb.attr["gene_id"]  << ": skipping.\n";
      continue;
    }

    LOG(QGS::Log::VERBOSE) << "QGS: Outputting QGS for " << gb.attr["gene_id"] << " (" << gb.chr << ':' << gb.start << '-' << gb.stop << ") based on "  << snp_cnt << "/" << total_snp_cnt << " loci.\n";

    out_file << gb.attr["gene_name"] << delimiter << gb.attr["gene_id"] << delimiter << gb.chr << delimiter 
             << gb.start << delimiter << gb.stop << delimiter << sample_file->num_samples() << delimiter << reference_file->num_samples() << delimiter << snp_cnt << delimiter << total_snp_cnt;
    for (auto const & score : total_score) {
      if (score < 0)
        out_file << delimiter << "NaN"; // missing data point
      else
        out_file << delimiter << score / (2 * snp_cnt * reference_file->num_samples());
    }
    out_file << '\n';

  }


  std::cout << "\nRun completed\n";
}
