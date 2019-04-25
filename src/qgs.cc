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
   bool const fill_missings = args.find("fill-missings") != args.end();
   bool const allow_missings = fill_missings || args.find("allow-missings") != args.end();
   char const delimiter = args.find("delimiter") == args.end() ? ',' : args.find("delimiter")->second.at(0)[0];
   bool const output_variants = args.find("output_variants") != args.end();
   
   bool const weighted_calculation = args.find("weight-by") != args.end();
   std::string const weight_info_field = weighted_calculation ? args.find("weight-by")->second.at(0) : std::string();


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
  
  // Keep some statistics about how many variants are used for what
  struct Snp_counts {
    std::size_t read = 0,
                outside_blocks = 0,
                inside_blocks = 0,
                overlapping = 0,
                non_overlapping = 0,
                skipped = 0,
                filled_missing = 0,
                used = 0;
  } sample_counts, ref_counts;
  
  
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
      
      Snp_counts & data_counts = (!reference_locus.chr || reference_locus <= sample_locus) ? ref_counts : sample_counts;
      
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
        scores[sample_locus.chr][sample_locus.pos][sample_locus.ref] = {};
      }
      
      ++data_counts.read;
      
      if (sample_locus.chr != reference_locus.chr || sample_locus.pos != reference_locus.pos) {
        ++data_counts.non_overlapping;
        continue; // no hit
      }
      
      ++data_counts.overlapping;
        
      if (!snp_filter.empty()) {
        // there is an include or exclude filter
        
        auto itt = snp_filter.find(sample_locus.id);
        if (itt == snp_filter.end())
          itt = snp_filter.find(std::to_string(sample_locus.chr) + ":" + std::to_string(sample_locus.pos));
          
        if (snp_include_filter && itt == snp_filter.end()) {
          LOG(QGS::Log::TRACE) << "Filter: Sample locus not in include filter. Skipping.\n";
          ++sample_counts.skipped;
          ++ref_counts.skipped;
          continue;
        }

        if (!snp_include_filter && itt != snp_filter.end()) {
          LOG(QGS::Log::TRACE) << "Filter: Sample locus in exclude filter. Skipping.\n";
          ++sample_counts.skipped;
          ++ref_counts.skipped;
          continue;
        }
        
        LOG(QGS::Log::TRACE) << "Filter: Sample locus passed filter. Including.\n";
      }
        
      if (sample_locus.chr != gb.chr || sample_locus.pos < gb.start || sample_locus.pos > gb.stop) {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " does not lie within block. skipping.\n";
        ++sample_counts.outside_blocks;
        ++ref_counts.outside_blocks;
        continue;
      }
        
      LOG(QGS::Log::TRACE) << "Match: locus " 
        << sample_locus.chr << ":" << sample_locus.pos 
        << " found in sample and reference. Attempting deep read.\n";

      ++sample_counts.inside_blocks;
      ++ref_counts.inside_blocks;
        
      sample_locus.parse_alt();
      reference_locus.parse_alt();
      
      if (sample_locus.palt.empty() || reference_locus.palt.empty()) {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " is missing alt data: skipping locus.\n";
        ++sample_counts.skipped;
        ++ref_counts.skipped;
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
          ++sample_counts.skipped;
          ++ref_counts.skipped;
          continue;
      }

      if (!sample_file->deep_read(sample_locus) || sample_locus.maf < maf_limit) {
        LOG(QGS::Log::TRACE) << "Match: sample locus excluded. skipping locus.\n";
        sample_locus.clear();
        ++sample_counts.skipped;
        ++ref_counts.skipped;
        continue;
      }
      
      if (!reference_file->deep_read(reference_locus) || reference_locus.maf < maf_limit) {
        LOG(QGS::Log::TRACE) << "Match: reference locus excluded. skipping locus.\n";
        sample_locus.clear();
        reference_locus.clear();
        ++sample_counts.skipped;
        ++ref_counts.skipped;
        continue;
      }
      
      LOG(QGS::Log::DEBUG) << "QGS: locus " 
        << sample_locus.chr << ":" << sample_locus.pos 
        << " deep read successful. Calculating QGS.\n";
        
      if (fill_missings) {
        std::size_t miss_cnt = 0;
        for (auto & dosage : sample_locus.data_ds)
          if (dosage < 0) {
            ++miss_cnt;
            dosage = 0;
          }
        if (miss_cnt) {
          ++sample_counts.filled_missing;
          LOG(QGS::Log::DEBUG) << "Filled " << miss_cnt
            << " missing values.\n";
        }
      }
      
      scores[sample_locus.chr][sample_locus.pos][sample_locus.ref] = QGS::score(sample_locus, reference_locus);
      
      ++sample_counts.used;
      ++ref_counts.used;

      if (weighted_calculation) {
        // add weight, if any
        sample_locus.parse_info();
        auto const weight_itt = sample_locus.pinfo.find(weight_info_field);
        if (weight_itt == sample_locus.pinfo.end()) {
          LOG(QGS::Log::WARNING) << "QGS: weight INFO field `"
            <<  weight_info_field << "` not available for locus "
            << sample_locus << ". Not including weight.\n";
          continue;
        }
        try {
          float const weight = std::stof(weight_itt->second);
          for (auto & score : scores[sample_locus.chr][sample_locus.pos][sample_locus.ref])
            score *= weight;
        }
        catch (...) {
          LOG(QGS::Log::WARNING) << "QGS: weight INFO field `"
            <<  weight_info_field << "` not a number for locus "
            << sample_locus << ". Not including weight.\n";
        }
      }
    }
    
    std::map<std::size_t, std::map<std::string, std::vector<QGS::Gene_score>>>::iterator lower_bound = scores[gb.chr].lower_bound(gb.start);
    scores[gb.chr].erase(scores[gb.chr].begin(), lower_bound);
    
    std::vector<QGS::Gene_score> total_score;
    total_score.resize(sample_file->num_samples(), 0);
    std::ostringstream used_loci, unused_loci;
    std::size_t snp_cnt = 0, total_snp_cnt = 0;
    for (auto itt_pos : scores[gb.chr]) {
      if (itt_pos.first > gb.stop)
        break;
      for (auto itt_var : itt_pos.second) {
        ++total_snp_cnt;
        if (itt_var.second.empty()) {
          unused_loci << "|" << gb.chr << ":" << itt_pos.first << ":" << itt_var.first;
          continue;
        }
        used_loci << "|" << gb.chr << ":" << itt_pos.first << ":" << itt_var.first;
        ++snp_cnt;
        for (std::size_t sample_idx = 0; sample_idx != total_score.size(); ++sample_idx) {
          if (itt_var.second[sample_idx] < 0 || total_score[sample_idx] < 0)
            total_score[sample_idx] = -99;
          else
            total_score[sample_idx] += itt_var.second[sample_idx];
        }
      }
    }
    
    std::string used_loci_s = used_loci.str();
    if (!used_loci_s.empty())
      used_loci_s = used_loci_s.substr(1);
    std::string unused_loci_s = unused_loci.str();
    if (!unused_loci_s.empty())
      unused_loci_s = unused_loci_s.substr(1);
    
    LOG(QGS::Log::DEBUG) << "Loci included in " << gb.attr["gene_id"] << ": " << used_loci_s << "\n";
    LOG(QGS::Log::DEBUG) << "Loci not included in " << gb.attr["gene_id"] << ": " << unused_loci_s << "\n";

    if (!snp_cnt) {
      LOG(QGS::Log::TRACE) << "Gene: No loci in " << gb.attr["gene_id"]  << ": skipping.\n";
      continue;
    }

    LOG(QGS::Log::VERBOSE) << "QGS: Outputting QGS for " << gb.attr["gene_id"] << " (" << gb.chr << ':' << gb.start << '-' << gb.stop << ") based on "  << snp_cnt << "/" << total_snp_cnt << " loci.\n";

    out_file << gb.attr["gene_name"] << delimiter << gb.attr["gene_id"] << delimiter << gb.chr << delimiter 
             << gb.start << delimiter << gb.stop << delimiter << sample_file->num_samples() << delimiter << reference_file->num_samples() << delimiter << (output_variants ? used_loci_s : std::to_string(snp_cnt)) << delimiter << total_snp_cnt;
    for (auto const & score : total_score) {
      if (score < 0)
        out_file << delimiter << "NaN"; // missing data point
      else
        out_file << delimiter << score / (2 * snp_cnt * reference_file->num_samples());
    }
    out_file << '\n';

  }
  
  std::cout << std::setprecision(3);
  std::cout << "Sample statistics:\n" 
    << "  Read: " << sample_counts.read << "\n"
    << "  Used: " << sample_counts.used << " (" 
    << (sample_counts.used * 100.0) / sample_counts.read << "%)\n"
    << "  Overlapping: " << sample_counts.overlapping << " (" 
    << (sample_counts.overlapping * 100.0) / sample_counts.read << "%)\n"
//    << "  Non-overlapping: " << sample_counts.non_overlapping << " (" 
//    << (sample_counts.non_overlapping * 100.0) / sample_counts.read << "%)\n"
    << "  Skipped: " << sample_counts.skipped << " (" 
    << (sample_counts.skipped * 100.0) / sample_counts.read << "%)\n"
    << "  Inside regions: " << sample_counts.inside_blocks << " (" 
    << (sample_counts.inside_blocks * 100.0) / sample_counts.read << "%)\n"
    << "  Outside regions: " << sample_counts.outside_blocks << " (" 
    << (sample_counts.outside_blocks * 100.0) / sample_counts.read << "%)\n"
    << "  With missings: " << sample_counts.filled_missing << " (" 
    << (sample_counts.filled_missing * 100.0) / sample_counts.read << "%)\n\n";
    
  std::cout << "Reference statistics:\n" 
    << "  Read: " << ref_counts.read << "\n"
    << "  Used: " << sample_counts.used << " (" 
    << (ref_counts.used * 100.0) / ref_counts.read << "%)\n"
    << "  Overlapping: " << ref_counts.overlapping << " (" 
    << (ref_counts.overlapping * 100.0) / ref_counts.read << "%)\n"
//    << "  Non-overlapping: " << ref_counts.non_overlapping << " (" 
//    << (ref_counts.non_overlapping * 100.0) / ref_counts.read << "%)\n"
    << "  Skipped: " << ref_counts.skipped << " (" 
    << (ref_counts.skipped * 100.0) / ref_counts.read << "%)\n"
    << "  Inside regions: " << ref_counts.inside_blocks << " (" 
    << (ref_counts.inside_blocks * 100.0) / ref_counts.read << "%)\n"
    << "  Outside regions: " << ref_counts.outside_blocks << " (" 
    << (ref_counts.outside_blocks * 100.0) / ref_counts.read << "%)\n\n";

  std::cout << "Run completed\n";
}
