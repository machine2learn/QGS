/*
  GNU General Public License v3 (GPLv3)
  Copyright (c) 2016-2020 Gido Schoenmacker
  Complete License at https://github.com/machine2learn/QGS/blob/master/LICENSE.md
*/

#ifndef INCL_GENE_SCORE_H
#define INCL_GENE_SCORE_H

#include "snpreader.h"
#include <vector>

namespace QGS {

using Gene_score = float;

std::vector<Gene_score> score(SNPreader::Locus const &, SNPreader::Locus const &);

}
  
#endif
