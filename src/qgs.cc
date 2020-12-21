/*
  GNU General Public License v3 (GPLv3)
  Copyright (c) 2016-2020 Gido Schoenmacker
  Complete License at https://github.com/machine2learn/QGS/blob/master/LICENSE.md
*/

#include "qgs.h"

int main(int argc, char ** argv) {
  

  // GPL COPYRIGHT NOTICE
  
  std::cout << "QGS " << QGS_VERSION << " (C) Gido Schoenmacker.\n"
    "This program comes with ABSOLUTELY NO WARRANTY.\n"
    "This is free software and you are welcome to redistribute it "
    "under GPLv3 conditions.\n\n";
  
  // ##### SET COMMAND LINE CONSTANTS
  std::unordered_map<std::string, std::vector<std::string>> const args = QGS::command_line_arguments(argc, argv);
  
  // required arguments
  std::vector<std::string> const & fname_sample = args.find("sample")->second;
  std::vector<std::string> const & fname_ref = args.find("reference")->second;
  std::string const & fname_genes = args.find("genes")->second.at(0);
  std::string const & fname_out = args.find("out")->second.at(0);
  
  // default arguments
  std::size_t const flank = std::stoull(args.find("flank")->second.at(0));
  std::size_t const pre_flank = 1000 * (args.find("pre-flank") != args.end() ? std::stoull(args.find("pre-flank")->second.at(0)) : flank);
  std::size_t const post_flank = 1000 * (args.find("post-flank") != args.end() ? std::stoull(args.find("post-flank")->second.at(0)) : flank);
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
  bool const output_variants = args.find("output-variants") != args.end();

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
  
  GZofile out_file(fname_out);
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
  
  if (!out_file) {
    LOG(QGS::Log::FATAL) << "Failed to write to `" 
      << fname_out << "`.\n";
    return 1;
  }
  
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

  struct Scores {
    std::vector<QGS::Gene_score> scores;
    float weight = 1;
    std::string id;
  };
  
  std::map<std::size_t, std::map<std::size_t, std::map<std::string, Scores>>> scores;

  for (;;) {
    
    if (!(*sample_file && *reference_file)) {
      std::size_t num_snps_in_mem = 0;
      for (auto const & map_chr : scores)
        num_snps_in_mem += map_chr.second.size();
      LOG(QGS::Log::TRACE) << "No more genetic data left to read. "
        << num_snps_in_mem << " loci still in memory.\n";
      if (scores.empty() || gb.chr != reference_locus.chr || gb.start > reference_locus.pos) {
        LOG(QGS::Log::TRACE) << "Genetic region starts after last read. Stopping.\n";
        break;
      }
    }

    if (!(gene_file >> gb)) {
      if (gene_file.handle().eof()) {
        LOG(QGS::Log::TRACE) << "No more genetic regions left to read. Stopping.\n";
        break;
      }
      gene_file.handle().clear();
      LOG(QGS::Log::TRACE) << "Failed to read line from gene file. Skipping.\n";
      continue;
    }

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
      
      if (!(*sample_file && *reference_file))
        break;
      
      Snp_counts & data_counts = (!reference_locus.chr || reference_locus <= sample_locus) ? ref_counts : sample_counts;
      
      if (reference_locus.chr && gb < reference_locus && gb < sample_locus) {
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
        
      if (sample_locus.chr != gb.chr || sample_locus.pos < gb.start) {
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
      
      if (sample_locus.ref != reference_locus.ref && 
          sample_locus.ref == reference_locus.palt[0] && 
          reference_locus.ref == sample_locus.palt[0])
      {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " has flipped ref/alt data. Correcting.\n";
          
          // erase variant with wrong ref from read history
          scores[sample_locus.chr][sample_locus.pos].erase(sample_locus.ref);

          sample_locus.switch_alt_ref();
      }
      else if (sample_locus.ref != reference_locus.ref || reference_locus.palt[0] != sample_locus.palt[0]) {
        LOG(QGS::Log::TRACE) << "Match: locus " 
          << sample_locus.chr << ":" << sample_locus.pos 
          << " has ref/alt mismatch between sample and reference. skipping locus.\n";
          sample_locus.clear();
          reference_locus.clear();
          ++sample_counts.skipped;
          ++ref_counts.skipped;
          continue;
      }

      if (!sample_file->deep_read(sample_locus) || sample_locus.maf < maf_limit) {
        if (sample_locus.maf < maf_limit)
          LOG(QGS::Log::TRACE) << "Match: sample locus MAF below filter (" << sample_locus.maf << ").\n";
        else
          LOG(QGS::Log::TRACE) << "Match: sample locus excluded. Skipping locus.\n";
        sample_locus.clear();
        ++sample_counts.skipped;
        ++ref_counts.skipped;
        continue;
      }
      
      if (!reference_file->deep_read(reference_locus) || reference_locus.maf < maf_limit) {
        if (reference_locus.maf < maf_limit)
          LOG(QGS::Log::TRACE) << "Match: reference locus MAF below filter (" << reference_locus.maf << ").\n";
        else
          LOG(QGS::Log::TRACE) << "Match: reference locus excluded. Skipping locus.\n";
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
          if (std::isnan(dosage)) {
            ++miss_cnt;
            dosage = sample_locus.switch_ar ? 2 : 0;
          }
        if (miss_cnt) {
          ++sample_counts.filled_missing;
          LOG(QGS::Log::DEBUG) << "Filled " << miss_cnt
            << " missing values.\n";
        }
      }
      
      scores[sample_locus.chr][sample_locus.pos][sample_locus.ref].scores = QGS::score(sample_locus, reference_locus);
      scores[sample_locus.chr][sample_locus.pos][sample_locus.ref].id = sample_locus.id + "/" + reference_locus.id;
      
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
          for (auto & score : scores[sample_locus.chr][sample_locus.pos][sample_locus.ref].scores)
            score *= weight;
          scores[sample_locus.chr][sample_locus.pos][sample_locus.ref].weight = weight;
        }
        catch (...) {
          LOG(QGS::Log::WARNING) << "QGS: weight INFO field `"
            <<  weight_info_field << "` not a number for locus "
            << sample_locus << ". Not including weight.\n";
        }
      }
    }
    
    auto lower_bound = scores[gb.chr].lower_bound(gb.start);
    scores[gb.chr].erase(scores[gb.chr].begin(), lower_bound);
    
    std::vector<QGS::Gene_score> total_score;
    total_score.resize(sample_file->num_samples(), 0);
    std::ostringstream used_loci, unused_loci;
    std::size_t snp_cnt = 0, total_snp_cnt = 0;
    long double correction_factor = 0, addition_factor = 0;
    for (auto itt_pos : scores[gb.chr]) {
      if (itt_pos.first > gb.stop)
        break;
      for (auto itt_var : itt_pos.second) {
        ++total_snp_cnt;
        if (itt_var.second.scores.empty()) {
          unused_loci << "|" << gb.chr << ":" << itt_pos.first << ":" << itt_var.first;
          continue;
        }
        used_loci << "|" << gb.chr << ":" << itt_pos.first << ":" << itt_var.first << "(" << itt_var.second.id << ")";
        ++snp_cnt;
        if (itt_var.second.weight < 0)
          addition_factor += -2 * itt_var.second.weight;
        else
          correction_factor += 2 * itt_var.second.weight;
        for (std::size_t sample_idx = 0; sample_idx != total_score.size(); ++sample_idx) {
          if (std::isnan(itt_var.second.scores[sample_idx]) || std::isnan(total_score[sample_idx]))
            total_score[sample_idx] = NAN;
          else {
            total_score[sample_idx] += itt_var.second.scores[sample_idx];
          }
        }
      }
    }
    
    std::string used_loci_s = used_loci.str();
    if (!used_loci_s.empty())
      used_loci_s = used_loci_s.substr(1);
    std::string unused_loci_s = unused_loci.str();
    if (!unused_loci_s.empty())
      unused_loci_s = unused_loci_s.substr(1);
    
    if (!snp_cnt) {
      LOG(QGS::Log::TRACE) << "Gene: No loci in " << gb.attr["gene_id"]  << ": skipping.\n";
      continue;
    }

    double constexpr eps = 0.00001; // some small value
    if (correction_factor < eps && addition_factor < eps) {
      LOG(QGS::Log::TRACE) << "Gene: No loci with sufficient weights in " << gb.attr["gene_id"]  << ": skipping.\n";
      continue;
    }

    LOG(QGS::Log::DEBUG) << "Loci included in " << gb.attr["gene_id"] << ": " << used_loci_s << "\n";
    LOG(QGS::Log::DEBUG) << "Loci not included in " << gb.attr["gene_id"] << ": " << unused_loci_s << "\n";

    LOG(QGS::Log::VERBOSE) << "QGS: Outputting QGS for " << gb.attr["gene_id"] << " (" << gb.chr << ':' << gb.start << '-' << gb.stop << ") based on "  << snp_cnt << "/" << total_snp_cnt << " loci.\n";

    out_file << gb.attr["gene_name"] << delimiter << gb.attr["gene_id"] << delimiter << gb.chr << delimiter 
             << gb.start << delimiter << gb.stop << delimiter << sample_file->num_samples() << delimiter << reference_file->num_samples() << delimiter << (output_variants ? used_loci_s : std::to_string(snp_cnt)) << delimiter << total_snp_cnt;

    for (auto const & score : total_score) {
      if (std::isnan(score))
        out_file << delimiter << "NaN"; // missing data point
      else
        out_file << delimiter << (score + addition_factor * reference_file->num_samples()) / (correction_factor * reference_file->num_samples() + addition_factor * reference_file->num_samples());
    }
    if (!(out_file << '\n')) { // final write of this QGS score, check for success
      LOG(QGS::Log::FATAL) << "Failed to write to `" 
        << fname_out << "`.\n";
      return 1;
    }
  }
  
  // output some statistics

  sample_counts.overlapping += ref_counts.overlapping;
  
  std::cout << std::setprecision(3);
  LOG(QGS::Log::INFO) << "Sample statistics:\n" 
    << "  Loci read: " << sample_counts.read << "\n"
    << "  Loci used: " << sample_counts.used << " (" 
    << (sample_counts.used * 100.0) / sample_counts.read << "%)\n";

  LOG(QGS::Log::VERBOSE)
    << "Overlapping: " << sample_counts.overlapping << " (" 
    << (sample_counts.overlapping * 100.0) / sample_counts.read << "%)\n"
    << "  Skipped: " << sample_counts.skipped << " (" 
    << (sample_counts.skipped * 100.0) / sample_counts.read << "%)\n"
    << "  Inside regions: " << sample_counts.inside_blocks << " (" 
    << (sample_counts.inside_blocks * 100.0) / sample_counts.read << "%)\n"
    << "  Outside regions: " << sample_counts.outside_blocks << " (" 
    << (sample_counts.outside_blocks * 100.0) / sample_counts.read << "%)\n"
    << "  With missings: " << sample_counts.filled_missing << " (" 
    << (sample_counts.filled_missing * 100.0) / sample_counts.read << "%)\n\n";

  std::cout << "Run completed\n";
}
