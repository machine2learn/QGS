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

#include "gene_score.h"
#include <iostream>
#include <cmath>
#include <iterator>
#include <numeric>
#include "log.h"

namespace QGS {

std::vector<Gene_score> score(
  SNPreader::Locus const & sample_locus, 
  SNPreader::Locus const & reference_locus)
{
  std::vector<Gene_score> out;
  out.resize(sample_locus.data_ds.size(), 0);

  std::unordered_map<float, Gene_score> cache;

  for (std::size_t sample_idx = 0; sample_idx < sample_locus.data_ds.size(); ++sample_idx) {
    
    auto const dosage = sample_locus.data_ds[sample_idx];
    
    if (std::isnan(dosage)) {
      out[sample_idx] = NAN;
      continue; // missing data point, set as NAN
    }

    // check to see if dosage is in cache
    auto itt = cache.find(dosage);
    if (itt != cache.end()) {
      out[sample_idx] = itt->second;
      continue;
    }

    // else, calculate
    Gene_score gscore = 0;
    for (std::size_t reference_idx = 0; reference_idx < reference_locus.data_ds.size(); ++reference_idx) {
      if (reference_locus.data_ds[reference_idx] < 0) {
        LOG(QGS::Log::TRACE) << "QGS: Reference indiviual "
          << reference_idx << " ignored because of missing data point"
          << ". suggest removal of missing data points in reference.\n";
        continue;
      }
      gscore += dosage > reference_locus.data_ds[reference_idx] ? 
        dosage - reference_locus.data_ds[reference_idx]
          :
        reference_locus.data_ds[reference_idx] - dosage;
    }

    cache[dosage] = gscore;
    out[sample_idx] = gscore;
  }

  return out;
}

}
