#include <iostream>
#include "gzfile.h"

int main(int argc, char ** argv) {

  // Open output file
  std::ofstream out_file(argv[1]);
  if (!out_file) {
    std::cerr << "Can't open output file for writing. Abort.\n";
    return 1;
  }
}
