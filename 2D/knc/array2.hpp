/*
 *  array2.hpp
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __ARRAY2_HPP__
#define __ARRAY2_HPP__

#include <cstdio>
#include <errno.h>

namespace fixedgrid {

template < typename T >
class Array2
{
public:

  typedef T element_type;

  Array2(size_t _nrows, size_t _ncols, element_type const & init) :
    data(new element_type[_nrows*_ncols]),
    rows(new element_type*[_nrows]),
    nrows(_nrows), ncols(_ncols), count(_nrows*_ncols)
  {
    for (int i=0; i<count; ++i) {
      data[i] = init;
    }
    for (int i=0; i<nrows; ++i) {
      rows[i] = data + i*ncols;
    }
  }

  ~Array2() {
    delete[] rows;
    delete[] data;
  }

  element_type * operator()(size_t const i) {
    return rows[i];
  }
  element_type const * operator()(size_t const i) const {
    return rows[i];
  }

  element_type & operator()(size_t const i, size_t const j) {
    return rows[i][j];
  }
  element_type const & operator()(size_t const i, size_t const j) const {
    return rows[i][j];
  }

  /*
   * From gnuplot/doc/html/gnuplot.html:
   * For binary input data, the first element of the first row must contain
   * the number of data columns. (This number is ignored for ascii input).
   * Both the coordinates and the data values in a binary input are treated
   * as single precision floats. Example commands for plotting non-uniform matrix data:
   *
   *     splot 'file' nonuniform matrix using 1:2:3  # ascii input
   *     splot 'file' binary matrix using 1:2:3      # binary input
   *
   * Thus the data organization for non-uniform matrix input is
   *
   *      <N+1>  <y0>   <y1>   <y2> ...  <yN>
   *      <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
   *      <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
   *       :      :      :      :   ...    :
   *
   * which is then converted into triplets:
   *
   *      <x0> <y0> <z0,0>
   *      <x0> <y1> <z0,1>
   *      <x0> <y2> <z0,2>
   *       :    :     :
   *      <x0> <yN> <z0,N>
   *      <x1> <y0> <z1,0>
   *      <x1> <y1> <z1,1>
   *       :    :     :
   *
   * These triplets are then converted into gnuplot iso-curves and then gnuplot
   * proceeds in the usual manner to do the rest of the plotting.
   */
  int WriteGnuplotBinaryMatrixFile(char const * fname)
  {
    using namespace std;

    // Open file for writing
    FILE * fout = fopen(fname, "wb");
    if (!fout) {
      int err = errno;
      fprintf(stderr, "ERROR: Can't open '%s' for writing\n", fname);
      return err;
    }

    size_t nmemb = ncols + 1;
    float buffer[nmemb];

    // Write ncols and column indices as float values
    buffer[0] = (float)ncols;
    for (int j=0; j<ncols; ++j) {
      buffer[j+1] = (float)j;
    }
    if (fwrite(buffer, sizeof(float), nmemb, fout) != nmemb) {
      fprintf(stderr, "ERROR: Failed to write %ld elements to '%s'\n", nmemb, fname);
      return -1;
    }

    // Write matrix rows
    for (int i=0; i<nrows; ++i) {
      buffer[0] = (float)i;
      for (int j=0; j<ncols; ++j) {
        buffer[j+1] = (float)rows[i][j];
      }
      if (fwrite(buffer, sizeof(float), nmemb, fout) != nmemb) {
        fprintf(stderr, "ERROR: Failed to write %ld elements to '%s'\n", nmemb, fname);
        return -1;
      }
    }

    // Cleanup and return
    fclose(fout);
    return 0;
  }

private:

  element_type * const data;
  element_type ** const rows;
  size_t const nrows;
  size_t const ncols;
  size_t const count;
};


} // namespace fixedgrid

#endif
