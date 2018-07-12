#ifndef INCL_GENE_SCORE_H
#define INCL_GENE_SCORE_H

#include "snpreader.h"
#include <vector>

namespace QGS {

using Gene_score = long double;

std::vector<Gene_score> score(SNPreader::Locus const &, SNPreader::Locus const &);

}
  
#endif
