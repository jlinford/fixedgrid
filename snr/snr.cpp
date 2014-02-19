/*
 * snr.cpp
 *
 *  Created on: Aug 8, 2011
 *      Author: jlinford
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libgen.h>

#include "metrics.h"

using namespace std;

bool is_ascii(istream & in)
{
  int c;
  while ((c = in.get()) != EOF && c < 128) 
    ; // Intentional
  // Reset and rewind
  in.clear();
  in.seekg(0, in.beg);
  // All bytes in the file < 128?
  return c == EOF;
}


int main(int argc, char ** argv)
{
  // Check arguments
  if(argc < 3) {
    cout << "Usage: " << argv[0] << " file1 file2" << endl;
    return 1;
  }
  char const * fname1 = argv[1];
  char const * fname2 = argv[2];

  // Open files for reading
  ifstream file1(argv[1]);
  if (!file1) {
    cerr << "ERROR: Can't read " << fname1 << endl;
    return 1;
  }

  ifstream file2(argv[2]);
  if (!file2) {
    cerr << "ERROR: Can't read " << fname2 << endl;
    return 1;
  }

  bool ascii = is_ascii(file1);
  if (ascii != is_ascii(file2)) {
    cerr << "ERROR: Different file encodings" << endl;
    return 1;
  }
  if (ascii) {
    cout << "Comparing text files" << endl;
  } else {
    cout << "Comparing binary files" << endl;
  }

  Metrics m(ascii);
  while (file1.good() && file2.good()) {
    m.update(file1, file2);
  }
  m.summary();

  if (m.good()) {
    cout << "PASSED" << endl;
    return 0;
  } else {
    cout << "FAILED" << endl;
    return 1;
  }
}
