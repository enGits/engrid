// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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

#ifndef SmallSquareMatrix_HH
#define SmallSquareMatrix_HH

#include "mathvector.h"

template <class T, uint_t N> class SmallSquareMatrix;

#include <vector>
#include <cstdarg>
#include <cstdlib>

template <class T, uint_t N>
class SmallSquareMatrix : public MathVector<StaticVector<MathVector<StaticVector<double,N> >,N> >
{
  typedef StaticVector<T,N> svec_t;
  typedef MathVector<StaticVector<double,N> > rvec_t;

protected:

  /** Indicator for precision safe mode
   */
  bool prec_safe;
  T prec_limit;

public:

  /** get a component of the matrix.
   *  @param row the row of the component
   *  @param column the column of the component
   *  @return the component M[row,column]
   */
  T comp(uint_t row, uint_t column) const { return (*this)[row][column]; }

  /** Empty Constructor I */
  SmallSquareMatrix<T,N>(T a_limit = 1e-20)
    : MathVector<StaticVector<MathVector<StaticVector<double,N> >,N> >() { prec_safe = false; a_limit = a_limit; }

  /** Copy Constructor */
  SmallSquareMatrix<T,N>(const SmallSquareMatrix<T,N> &M)
    : MathVector<StaticVector<MathVector<StaticVector<double,N> >,N> >(M) {prec_safe = false;}

  /** Constructor upon an open set of vectors
   *  @param row_col_type string indicating "row" or "column" type setting
   *         (only these two key-words allowed)
   *  @param col_vects an stl-type vector of double vectors
   */
  SmallSquareMatrix<T,N>(string row_col_type, vector<svec_t*> col_vects);

  /** Constructor upon an open set of vectors containing doubles [T = double]
   *  @param row_col_type string indicating "row" or "column" type setting
   *         (only these two key-words allowed)
   *  @param col_vects an stl-type vector of double vectors
   */
  SmallSquareMatrix<T,N>(string row_col_type, vector<rvec_t*> col_vects);

  /** Set a row vector
   *  @param row the row number to be set
   *  @param row_vec a double vector to set
   */
  template<class Tvec>
  void row(uint_t row, Tvec row_vec)
  {
    (*this)[row] = row_vec;
  }

  /** Set a column vector
   *  @param column the column number to be set
   *  @param col_vec a double vector to set
   */
  template<class Tvec>
  void column(uint_t column, Tvec col_vec)
  {
    for(uint_t k_row = 0; k_row<N; k_row++) {
      (*this)[k_row][column] = col_vec[k_row];
    }
  }

  /** Set safe mode. Safe mode throws a Precision_error in
   *  case of a precision limit (e.g. floating exception for "double")
   *  @param a_prec_limit precision limit
   */
  void setSafe(T a_prec_limit) {
    prec_limit = a_prec_limit;
    prec_safe = true;
  }

  /** Set unsafe mode. Safe mode throws a Precision_error in
   *  case of a precision limit (e.g. floating exception for "double")
   */
  void setUnSafe() {
    prec_safe = false;
  }

  /** Access safe mode. Safe mode throws a Precision_error in
   *  case of a precision limit (e.g. floating exception for "double")
   */
  bool isSafe() {
    return prec_safe;
  }

  /** Access safe limit. Safe mode throws a Precision_error in
   *  case of a precision limit (e.g. floating exception for "double")
   */
  T safeLimit() {
    return prec_limit;
  }

  /** Analyse Precision and throw error, if needed
   *  @param det determinant
   *  @param ele_max a linear norm
   */
  void precisionHandling(T det, T ele_max) {
    T qdet, compare_det;
    qdet = det * det;
    double expo = 0.5 / N;
    compare_det = pow(qdet, expo);
    //.... Compare the det^(1/N) with ele_max
    if(compare_det < safeLimit() * ele_max) {
      cerr << "matrix appears to be singular within given precision" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /** Get maximum absolute value of all elements
   *  @param det determinant
   *  @param ele_max a linear norm
   */
  T linNorm_0() {
    uint_t i,j;
    T ele_max, qq;
    ele_max = (*this)[0][0] * (*this)[0][0];
    for(i=0;i<N;i++) {
      for(j=0;j<N;j++) {
	qq = (*this)[i][j] * (*this)[i][j];
	if(qq > ele_max) {
	  ele_max = qq;
        }
      }
    }
    ele_max = sqrt(ele_max);
    return ele_max;
  }

  /** get an identity matrix
   *  @return an identity matrix of the dimension N
   */ 
  static SmallSquareMatrix<T,N> identity();
  
  /** get the determinant of the matrix
   *  @return the determinant
   */
  T det();

  /** get the inverse of the matrix
   *  @return the inverse of the matrix
   */ 
  SmallSquareMatrix<T,N> inverse();

  /** get a sub-matrix.
   *  This returns a submatrix where one row and one column
   *  from the original matrix will be missing.
   *  @param row the row to exclude
   *  @param column the column to exclude
   *  @return the sub-matrix
   */
  SmallSquareMatrix<T,N-1> subMatrix(uint_t row, uint_t column);

  /** get the transposed matrix.
   *  @return the transposed matrix
   */
  SmallSquareMatrix<T,N> transp();

  /// matrix vector multiplication operator
  MathVector<StaticVector<T,N> > operator* (const MathVector<StaticVector<T,N> > &vec) const;

  /// matrix matrix multiplication operator
  SmallSquareMatrix<T,N> operator* (const SmallSquareMatrix<T,N> &mat) const;

  /** Initialize all matrix elements with initvalue
   *  @param initvalue a value
   */
  void initAll(double initvalue);

};

template <class T, uint_t N>
inline MathVector<StaticVector<T,N> > 
SmallSquareMatrix<T,N>::operator* (const MathVector<StaticVector<T,N> > &vec) const 
{
  MathVector<StaticVector<T,N> > result_vec;
  for (uint_t i = 0; i < N; i++) {
    result_vec[i] = 0;
    for (uint_t j = 0; j < N; j++)
      result_vec[i] += comp(i,j) * vec[j];
  }
  return result_vec;
}

template <class T, uint_t N>
inline SmallSquareMatrix<T,N> SmallSquareMatrix<T,N>::operator* (const SmallSquareMatrix<T,N> &mat) const 
{
  SmallSquareMatrix<T,N> result_mat;
  for (uint_t i = 0; i < N; ++i) {
    for (uint_t j = 0; j < N; ++j) {
      result_mat[i][j] = 0;
      for (uint_t k = 0; k < N; ++k) result_mat[i][j] += this->comp(i,k)*mat[k][j];
    }
  }
  return result_mat;
}

template <class T, uint_t N>
inline void SmallSquareMatrix<T,N>::initAll(double initvalue)
{
  for (uint_t i = 0; i < N; i++) {
    for (uint_t j = 0; j < N; j++) {
      (*this)[i][j] = initvalue;
    }
  }
}

template <class T, uint_t N>
inline T SmallSquareMatrix<T,N>::det()
{
  // Determinant of the matrix
  // copy yourself to protect matrix entries
  SmallSquareMatrix<T,N> a = *this;

  int n;
  n=N;
  int k,i,j;
  vector<int> p(n);
  T q,s,max,h,det;

  det=1;
  for(k=0;k<n-1;k++)
    {
      max=0.0;
      p[k]=0;
      for(i=k;i<n;i++)
	{
	  s=0.0;
	  for(j=k;j<n;j++)
	    {
	      s=s+fabs(a[i][j]);
	    }
	  q=fabs(a[i][k])/s;
	  if(q>max)
	    {
	      max=q;
	      p[k]=i;
	    }
	}
      if(!(p[k]==k))
	{
	  det=-det;
	  for(j=0;j<n;j++)
	    {
	      h=a[k][j];
	      a[k][j]=a[p[k]][j];
	      a[p[k]][j]=h;
	    }
	}
      det=det*a[k][k];
      for(i=k+1;i<n;i++)
	{
	  a[i][k]=a[i][k]/a[k][k];
	  for(j=k+1;j<n;j++)
	    a[i][j]=a[i][j]-a[i][k]*a[k][j];
	}
    }
  det=det*a[n-1][n-1];
    
  return det;
}

// Rainers inverter
template <class T, uint_t N>
class InvSmallSquareMatrix
{
protected:
  SmallSquareMatrix<T,N> b;

public:
  /** constructor.
   */ 
  InvSmallSquareMatrix<T,N>(SmallSquareMatrix<T,N> a,
			    bool a_prec_safe, T a_prec_limit) {

    int Smalldim = N;
    int n;
    n=N;
    int k,i,j,l;
    vector<int> p(n);
    T q,s,max,h,det;
    T ele_max = 0;
      
    if(a_prec_safe) {a.setSafe(a_prec_limit);}

    //.. Find maximum element to get a relative value
    if(a.isSafe()) {
      ele_max = a.linNorm_0();
    }

    //.. Get in matrix reduction
    for(k=0;k<Smalldim;k++)
      for(j=0;j<Smalldim;j++)
	b[j][k]=0.0;
    for(j=0;j<Smalldim;j++)
      b[j][j]=1.0;
    det=1;
    for(k=0;k<n-1;k++){
      max=0.0;
      p[k]=0;
      for(i=k;i<n;i++){
	s=0.0;
	for(j=k;j<n;j++) s=s+fabs(a[i][j]);
	q=fabs(a[i][k])/s;
	if(q>max){
	  max=q;
	  p[k]=i;
        }
      }
      if(!(p[k]==k)){
	det=-det;
	for(j=0;j<n;j++){
	  h=a[k][j];
	  a[k][j]=a[p[k]][j];
	  a[p[k]][j]=h;
        }
      }
      det=det*a[k][k];
      for(i=k+1;i<n;i++){
	a[i][k]=a[i][k]/a[k][k];
	for(j=k+1;j<n;j++) a[i][j]=a[i][j]-a[i][k]*a[k][j];
      }
    }
    det=det*a[n-1][n-1];

    //.. Proceed with rest of system reduction
    for(k=0;k<n-1;k++)
      if(!(p[k]==k)){
	for(l=0;l<n;l++){
	  h=b[k][l];
	  b[k][l]=b[p[k]][l];
	  b[p[k]][l]=h;
	}
      }
    for(i=0;i<n;i++){
      for(j=0;j<i;j++){
	for(l=0;l<n;l++)
	  b[i][l]=b[i][l]-a[i][j]*b[j][l];
      }
    }
    for(i=n-1;i>=0;i--){
      for(l=0;l<n;l++){
	s=b[i][l];
	for(k=i+1;k<n;k++)
	  s=s-a[i][k]*b[k][l];
	b[i][l]=s/a[i][i];
      }
    }

    //.. Check Determinant and throw error, if needed
    if(a.isSafe()) {
      a.precisionHandling(det, ele_max);
    }

  }

  SmallSquareMatrix<T,N> inverse() { return b; }
};


template <class T>
class InvSmallSquareMatrix<T,2>
{
protected:
  SmallSquareMatrix<T,2> INV;

public:
  // constructor.
    
  InvSmallSquareMatrix<T,2>(SmallSquareMatrix<T,2> SSM,
			    bool a_prec_safe, T a_prec_limit)
  {
    T ele_max;
    if(a_prec_safe) {
      SSM.setSafe(a_prec_limit);
      ele_max = SSM.linNorm_0();
      // sorce: maple
      T det = (SSM[0][0]*SSM[1][1]-SSM[0][1]*SSM[1][0]);
      SSM.precisionHandling(det, ele_max);
      T t4 = 1/det;
      INV[0][0] = SSM[1][1]*t4;
      INV[0][1] = -SSM[0][1]*t4;
      INV[1][0] = -SSM[1][0]*t4;
      INV[1][1] = SSM[0][0]*t4;
    }
    else {
      // sorce: maple
      T t4 = 1/(SSM[0][0]*SSM[1][1]-SSM[0][1]*SSM[1][0]);
      INV[0][0] = SSM[1][1]*t4;
      INV[0][1] = -SSM[0][1]*t4;
      INV[1][0] = -SSM[1][0]*t4;
      INV[1][1] = SSM[0][0]*t4;
    }
  }
    
  SmallSquareMatrix<T,2> inverse() { return INV; }
};

template <class T>
class InvSmallSquareMatrix<T,3>
{
protected:
  SmallSquareMatrix<T,3> INV;

public:
  // constructor.
   
  InvSmallSquareMatrix<T,3>(SmallSquareMatrix<T,3> SSM,
			    bool a_prec_safe, T a_prec_limit)
  {
    if(a_prec_safe) {
      SSM.setSafe(a_prec_limit);
      T ele_max = SSM.linNorm_0();
      // Source (maple)
      T t4 = SSM[0][0]*SSM[1][1];
      T t6 = SSM[0][0]*SSM[1][2];
      T t8 = SSM[0][1]*SSM[1][0];
      T t10 = SSM[0][2]*SSM[1][0];
      T t12 = SSM[0][1]*SSM[2][0];
      T t14 = SSM[0][2]*SSM[2][0];
      T det = (t4*SSM[2][2]-t6*SSM[2][1]-t8*SSM[2][2]+t10*SSM[2][1]+
	       t12*SSM[1][2]-t14*SSM[1][1]);
      SSM.precisionHandling(det, ele_max);
      T t17 = 1/det;
      INV[0][0] = (SSM[1][1]*SSM[2][2]-SSM[1][2]*SSM[2][1])*t17;
      INV[0][1] = -(SSM[0][1]*SSM[2][2]-SSM[0][2]*SSM[2][1])*t17;
      INV[0][2] = -(-SSM[0][1]*SSM[1][2]+SSM[0][2]*SSM[1][1])*t17;
      INV[1][0] = -(SSM[1][0]*SSM[2][2]-SSM[1][2]*SSM[2][0])*t17;
      INV[1][1] = (SSM[0][0]*SSM[2][2]-t14)*t17;
      INV[1][2] = -(t6-t10)*t17;
      INV[2][0] = -(-SSM[1][0]*SSM[2][1]+SSM[1][1]*SSM[2][0])*t17;
      INV[2][1] = -(SSM[0][0]*SSM[2][1]-t12)*t17;
      INV[2][2] = (t4-t8)*t17;
    }
    else {
      // Source (maple)
      T t4 = SSM[0][0]*SSM[1][1];
      T t6 = SSM[0][0]*SSM[1][2];
      T t8 = SSM[0][1]*SSM[1][0];
      T t10 = SSM[0][2]*SSM[1][0];
      T t12 = SSM[0][1]*SSM[2][0];
      T t14 = SSM[0][2]*SSM[2][0];
      T t17 = 1/(t4*SSM[2][2]-t6*SSM[2][1]-t8*SSM[2][2]+t10*SSM[2][1]+
		 t12*SSM[1][2]-t14*SSM[1][1]);
      INV[0][0] = (SSM[1][1]*SSM[2][2]-SSM[1][2]*SSM[2][1])*t17;
      INV[0][1] = -(SSM[0][1]*SSM[2][2]-SSM[0][2]*SSM[2][1])*t17;
      INV[0][2] = -(-SSM[0][1]*SSM[1][2]+SSM[0][2]*SSM[1][1])*t17;
      INV[1][0] = -(SSM[1][0]*SSM[2][2]-SSM[1][2]*SSM[2][0])*t17;
      INV[1][1] = (SSM[0][0]*SSM[2][2]-t14)*t17;
      INV[1][2] = -(t6-t10)*t17;
      INV[2][0] = -(-SSM[1][0]*SSM[2][1]+SSM[1][1]*SSM[2][0])*t17;
      INV[2][1] = -(SSM[0][0]*SSM[2][1]-t12)*t17;
      INV[2][2] = (t4-t8)*t17;
    }
  }
  SmallSquareMatrix<T,3> inverse() { return INV; }
};



typedef SmallSquareMatrix<double,2> mat2_t;
typedef SmallSquareMatrix<double,3> mat3_t;
typedef SmallSquareMatrix<double,4> mat4_t;
typedef SmallSquareMatrix<double,5> mat5_t;

template <class T, uint_t N>
SmallSquareMatrix<T,N>::SmallSquareMatrix(string row_col_type,
					  vector<svec_t*> col_vects)
  : MathVector<StaticVector<double,N> >()
{
  if(N < col_vects.size()) {
    cout
      << "SmallSquareMatrix<T,N>(string row_col_type, vector<SmallVector<T,N>*> col_vects)"
      << "\n"
      << "too many input vectors" << endl;
  }
  uint_t direct_it = 0;
  if(row_col_type == "Column") {
    for(typename vector<svec_t*>::iterator kk = col_vects.begin();
	kk != col_vects.end();
	kk++)
      {
	(*this).Column(direct_it, *(*kk));
	direct_it++;
      }
  } else if(row_col_type == "Row") {
    for(typename vector<svec_t*>::iterator kk = col_vects.begin();
	kk != col_vects.end();
	kk++)
      {
	(*this).Row(direct_it, *(*kk));
	direct_it++;
      }
  } else {
    cout
      << "SmallSquareMatrix<T,N>(string row_col_type, uint_t num_smvects, ...)"
      << "\n"
      << "Only Row or Column allowed as first argument" << endl;
    exit(EXIT_FAILURE);
  }
  prec_safe = false;
}

template <class T, uint_t N>
SmallSquareMatrix<T,N>::SmallSquareMatrix(string row_col_type,
					  vector<rvec_t*> col_vects)
  : MathVector<StaticVector<MathVector<StaticVector<double,N> >,N> >()
{
  if(N < col_vects.size()) {
    cout
      << "SmallSquareMatrix<real,N>(string row_col_type, vector<rvec_t*> col_vects)"
      << "\n"
      << "too many input vectors" << endl;
  }
  uint_t direct_it = 0;
  if(row_col_type == "Column") {
    for(typename vector<rvec_t*>::iterator kk = col_vects.begin();
	kk != col_vects.end();
	kk++)
      {
	(*this).Column(direct_it, *(*kk));
	direct_it++;
      }
  } else if(row_col_type == "Row") {
    for(typename vector<rvec_t*>::iterator kk = col_vects.begin();
	kk != col_vects.end();
	kk++)
      {
	(*this).Row(direct_it, *(*kk));
	direct_it++;
      }
  } else {
    cout
      <<  "SmallSquareMatrix<real,N>(string row_col_type, vector<rvec_t*> col_vects)"
      << "\n"
      << "Only Row or Column allowed as first argument" << endl;
    exit(EXIT_FAILURE);
  }
  prec_safe = false;
}

template <class T, uint_t N>
SmallSquareMatrix<T, N> SmallSquareMatrix<T, N>::identity()
{
  SmallSquareMatrix<T, N> I;
  for (uint_t i = 0; i < N; i++) {
    for (uint_t j = 0; j < N; j++) {
      if (i==j) I[i][j] = 1;
      else I[i][j] = 0;
    }
  }
  return I;
}

/*
template <class T, uint_t N>
T SmallSquareMatrix<T,N>::Det()
{
  // This construction is required, since a prtioal specializations on a
  // method with multiple template arguments does not work without specializing
  // the whole class :-(
  DetSmallSquareMatrix<T,N> DET(*this);
  return DET.Det();
}
*/

template <class T, uint_t N>
SmallSquareMatrix<T,N> SmallSquareMatrix<T,N>::transp()
{
 SmallSquareMatrix<T,N> M_t;
  for (uint_t i = 0; i < N; i++) {
    for (uint_t j = 0; j < N; j++) {
      M_t[i][j] = comp(j,i);
    }
  }
  return M_t;
}

template <class T, uint_t N>
SmallSquareMatrix<T,N> SmallSquareMatrix<T,N>::inverse()
{
  InvSmallSquareMatrix<T,N> INV(*this, isSafe(), prec_limit);
  return INV.inverse();
}

template <class T, uint_t N>
SmallSquareMatrix<T,N-1> SmallSquareMatrix<T,N>::subMatrix(uint_t row, uint_t column)
{
  uint_t i_new, j_new;
  i_new = 0;
  SmallSquareMatrix<T,N-1> M;
  for (uint_t i = 0; i < N; i++) {
    if (i == row) continue;
    j_new = 0;
    for (uint_t j = 0; j < N; j++) {
      if (j == column) continue;
      M[i_new][j_new] = comp(i,j);
      j_new++;
    }
    i_new++;
  }
  return M;
}

#endif
