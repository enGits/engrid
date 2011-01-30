//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#ifndef UTILITIES_H
#define UTILITIES_H

/** @file utilities.h Contains several utility functions. */

#include <QVector>
#include <QMap>
#include <QObject>
#include <QtDebug>
#include <QFileDialog>

#include <iostream>
using namespace std;

#include "math/mathvector.h"
#include "math/smallsquarematrix.h"
#include "math/linsolve.h"

#include "vtkUnstructuredGrid.h"

#include <complex>

#ifdef _MSC_VER
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) !_finite(x)
#endif

using namespace std;

/** Restricts value to the CYCLIC [min,max[ range.
 * Equivalent to min + modulo(value-min, max-min)
 * @param value value to restrict
 * @param min minimum
 * @param max maximum (if value=max, then min will be returned)
 */
inline double restrict_to(double value, double min, double max)
{
  return value - floor((value - min) / (max - min))*(max - min);
}

/// Converts degrees to radians
inline double RadiansFromDegrees(double a_deg)
{
  return(a_deg / 180.0*M_PI);
}

/// Converts radians to degrees
inline double DegreesFromRadians(double a_rad)
{
  return(a_rad / M_PI*180.0);
}

double toDouble(QString str);///< Equivalent of the QString::toDouble function, but it also replaces "," with "."

/** Converts a double to a string.
 * @param x double value to convert
 * @param separator symbol to use as decimal separator. It is recommended to always use the default: ".".
 */
QString toString(double x, QString separator = QObject::tr("."));

/// a modulo function which also works for negative numbers
inline int modulo(int a, int b)
{
  return((a % b + b) % b);
}

/// a modulo function which also works for negative numbers
inline double modulo(double a, double b)
{
  while (a < 0) a += b;
  while (a >= b) a -= b;
  return a;
}

/** Adds a suffix to a filename if it is missing.
 * @param filename string to which to add the suffix
 * @param suffix suffix to add. suffix should be of the form 'ext' without the dot!!!
 * @param remove_old_suffix If true, any previous suffix of filename will be removed. ex: foo.bar -> foo.png
 * @return the filename with the new suffix
 */
QString addSuffix(QString filename, QString suffix, bool remove_old_suffix);

///////////////////////////////////////////

/** Transposes a vector of vectors (matrix)
 * @param in the matrix to transpose
 * @param nrows the number of rows of the input matrix
 * @param ncolumns the number of columns of the input matrix
 * @return the transposed "matrix"
 */
template <class T>
QVector < QVector <T> > transpose(QVector < QVector <T> > & in, int nrows, int ncolumns)
{
  QVector < QVector <T> > out(ncolumns);
  QVector <T> col(nrows);
  out.fill(col);

  for (int i = 0; i < in.size(); i++) {
    QVector <T> row = in[i];
    for (int j = 0; j < row.size(); j++) {
      out[j][i] = in[i][j];
    }
  }
  return(out);
}

///////////////////////////////////////////

/// ostream operator for QVectors
template <class T>
ostream &operator<<(ostream &out, QVector<T> const & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    out << vector.at(i);
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

/// ostream operator for QSets
template <class T>
ostream &operator<<(ostream &out, QSet<T> const & set)
{
  out << "[ ";
  foreach(T value, set) out << value << " ";
  out << "]";
  return(out);
}

/// ostream operator for QVectors of QSets
template <class T>
ostream &operator<<(ostream &out, QVector<QSet<T> > & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    QSet<T> set = vector.at(i);
    out << set;
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

/// ostream operator for QVectors of QVectors
template <class T>
ostream &operator<<(ostream &out, QVector<QVector<T> > & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    QVector<T> subvector = vector.at(i);
    out << subvector;
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

/// ostream operator for QMaps
template <class T1, class T2>
ostream &operator<<(ostream &out, QMap<T1, T2> & map)
{
  QMapIterator<T1, T2> i(map);
  out << "[";
  while (i.hasNext()) {
    i.next();
    out << " [" << i.key() << ": " << i.value() << "]";
  }
  out << "]";
  return(out);
}

/// ostream operator for a QVector of pairs
template <class T1, class T2>
ostream &operator<<(ostream &out, QVector < pair<T1, T2> > & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    out << "<";
    out << vector.at(i).first;
    out << ",";
    out << vector.at(i).second;
    out << ">";
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

////////////////////////////////////////////////////

template <class T>
QVector <T> Set2Vector(QSet <T> a_set, bool a_sort)
{
  QVector <T> l_vector(a_set.size());
  qCopy(a_set.begin(),a_set.end(),l_vector.begin());
  if(a_sort) qSort(l_vector.begin(),l_vector.end());
  return(l_vector);
}

template <class T>
QSet <T> Vector2Set(QVector <T> a_vector, bool a_sort)
{
  QSet <T> l_set;
  foreach(T element, a_vector) l_set.insert(element);
  if(a_sort) qSort(l_set.begin(),l_set.end());
  return(l_set);
}

template <class T>
bool duplicates(QVector <T> a_vector)
{
  QSet <T> l_set;
  foreach(T element, a_vector) l_set.insert(element);
  return l_set.size()!=a_vector.size();
}

////////////////////////////////////////////////////

Qt::CheckState int2CheckState(int a);///< Converts an integer into a Qt::CheckState: Qt::Unchecked=0, Qt::PartiallyChecked=1, Qt::Checked=2
int CheckState2int(Qt::CheckState a);///< Converts a Qt::CheckState into an integer: Qt::Unchecked=0, Qt::PartiallyChecked=1, Qt::Checked=2

QString vector2csv(QVector <double> vector);///< Converts a vector into a CSV string
QVector <double> csv2vector(QString csv);///< Converts a CSV string into a vector

/// V(rotated base) = rotationMatrix_X * V(original base)
mat3_t rotationMatrix_X(double a_rad);

/// V(rotated base) = rotationMatrix_Y * V(original base)
mat3_t rotationMatrix_Y(double a_rad);

/// V(rotated base) = rotationMatrix_Z * V(original base)
mat3_t rotationMatrix_Z(double a_rad);

/// V(rotated base) = rotateRelativeZXZ * V(original base)
mat3_t rotateRelativeZXZ(double angle_1_rad, double angle_2_rad, double angle_3_rad);

/// V(rotated base) = rotateAbsoluteZXY * V(original base)
mat3_t rotateAbsoluteZXY(double angle_1_rad, double angle_2_rad, double angle_3_rad);

/// returns the polar angle in radians
double getGamma(vec3_t V);

/// returns the  azimuth angle in radians. returns phi from -pi to pi
double getPhi(vec3_t V);

/// A version of QFileDialog::getExistingDirectory which allows preselecting a directory
QString getDirectory(QWidget * parent = 0, const QString & caption = QString(), const QString & selected = QString());
// QString getDirectory( QWidget * parent = 0, const QString & caption = QString(), const QString & dir = QString(), const QString & selected = QString() , Options options = ShowDirsOnly );

/**
  * Utility function that allows printing selected data from an vtkUnstructuredGrid to any ostream (includes ofstream objects)
  * @param stream ostream object to print to
  * @param grid vtkUnstructuredGrid you want to print data from
  * @param npoints print number of points in the grid
  * @param ncells print number of cells in the grid
  * @param points print points in the grid
  * @param cells print cells in the grid
  */
int cout_grid(ostream &stream, vtkUnstructuredGrid *grid, bool npoints=true, bool ncells=true, bool points=false, bool cells=false);

///////////////////////////////////////////
int addCell(vtkUnstructuredGrid* a_grid, vtkIdType A, vtkIdType B, vtkIdType C, int bc);

///get number of the shortest side of the cell
int getShortestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid);
///get number of the longest side of the cell
int getLongestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid);
//sort sides by length
//QVector <vtkIdType> sortSidesByLength(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid);

///get number of the edge corresponding to node1-node2
int getSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid,vtkIdType a_id_node1,vtkIdType a_id_node2);

// QSet <int> complementary_bcs(QSet <int> &bcs, vtkUnstructuredGrid *a_grid, QVector <vtkIdType> &a_cells);
QString cell2str(vtkIdType id_cell,vtkUnstructuredGrid* grid);

///////////////////////////////////////////
///////////////////////////////////////////
pair<vtkIdType,vtkIdType> OrderedPair(vtkIdType a, vtkIdType b);

const char* VertexType2Str(char T);
char Str2VertexType(QString S);
QDebug operator<<(QDebug dbg, const vec3_t &v);
QDebug operator<<(QDebug dbg, const vec2_t &v);

bool checkVector(vec3_t V);
bool checkVector(vec2_t V);

/// returns the index of a node in a structured triangle grid
inline vtkIdType trigrid_idx(vtkIdType N, int i, int j) {
  int offset = -i * (i - 2 * N - 1) / 2;
  return offset + j;
}

/// returns the index of a node in a structured quad grid
inline vtkIdType quadgrid_idx(vtkIdType N, int i, int j) {
  return i*N + j;
}

// solver functions
typedef complex<double> dcmplx;
QDebug operator<<(QDebug dbg, const dcmplx &c);
dcmplx complex_pow(dcmplx base, double power);
// x^3 + a x^2 + b x + c = 0
int poly_solve_cubic(double a, double b, double c, double * x0, double * x1, double * x2);
// a x^2 + b x + c = 0
int poly_solve_quadratic(double a, double b, double c, double * x0, double * x1);

#endif
