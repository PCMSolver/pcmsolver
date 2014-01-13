#ifndef NPY_HPP
#define NPY_HPP
#include <fstream>
#include <iomanip>
#include <complex>

/*
  Functions for writing out a matrix in .npy format, described in
  https://github.com/numpy/numpy/blob/master/doc/neps/npy-format.txt

  fortran_order = true means that the columns of the matrix are close-packed in
  memory, otherwise the rows are close-packed.

  NOT TESTED on big endian machines, probably the '<f8' should be replaced by
  '>f8' there. Some kind of robust endian detection would be needed.

  Ulf Ekstrom <uekstrom@gmail.com>, December 2013.
 */

#ifdef __BIG_ENDIAN__
#warning FIXME: .npy writing will not work on this big endian machine.
#endif

static void npy_write_double_header(std::ostream &dst, size_t rows, size_t cols, bool fortran_order)
{
  const char *header;
  if (not fortran_order)
    header = "\x93NUMPY\x01\x00\x46\x00{'descr': '<f8', 'fortran_order': False, 'shape': (";
  else
    header = "\x93NUMPY\x01\x00\x46\x00{'descr': '<f8', 'fortran_order': True , 'shape': (";
  dst.write(header,61);
  dst << std::setfill('0') << std::setw(5) << rows << "," << std::setfill('0') << std::setw(5) << cols << "),}    \n";
}

static void npy_write_double_complex_header(std::ostream &dst, size_t rows, size_t cols, bool fortran_order)
{
  const char *header;
  if (not fortran_order)
    header = "\x93NUMPY\x01\x00\x46\x00{'descr': '<c16', 'fortran_order': False, 'shape': (";
  else
    header = "\x93NUMPY\x01\x00\x46\x00{'descr': '<c16', 'fortran_order': True , 'shape': (";
  dst.write(header,62);
  dst << std::setfill('0') << std::setw(5) << rows << "," << std::setfill('0') << std::setw(5) << cols << "),}   \n";
}

static void npy_write_matrix(std::ostream &dst, size_t rows, size_t cols, const double *data, bool fortran_order = false)
{
  npy_write_double_header(dst,rows,cols,fortran_order);
  dst.write(reinterpret_cast<const char *>(data),sizeof(*data)*rows*cols);    
}

static void npy_write_matrix(std::ostream &dst, size_t rows, size_t cols, const std::complex<double> *data, bool fortran_order = false)
{
  npy_write_double_complex_header(dst,rows,cols,fortran_order);
  dst.write(reinterpret_cast<const char *>(data),sizeof(*data)*rows*cols);    
}

#endif // NPY_HPP
