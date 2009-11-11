//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#ifndef engrid_H
#define engrid_H

#include <QMessageBox>
#include <QtDebug>
#include <QString>

#include <vtkSmartPointer.h>
#include <vtkLongArray.h>
#include <vtkLongLongArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include "error.h"
#include "math/mathvector.h"
#include "math/smallsquarematrix.h"
#include "geometrytools.h"
#include "containertricks.h"

#ifndef ENGRID_VERSION
#define ENGRID_VERSION "git"
#endif

#ifdef WIN32
typedef vtkLongLongArray vtkLongArray_t;
#else

#include <limits.h>
#if ( __WORDSIZE == 64 )
typedef vtkLongArray vtkLongArray_t;
#else
typedef vtkLongLongArray vtkLongArray_t;
#endif

#endif

#define EG_ERR_RETURN(TXT) \
{ \
  QString line; \
  line.setNum(__LINE__); \
  QString txt = QString(TXT) + "\n\nfile: " + __FILE__ + "\nline:" + line + "\n\n"; \
  Error err; \
  err.setText(txt); \
  err.setType(Error::ExitOperation); \
  throw err; \
};

#ifdef QT_DEBUG
#define EG_BUG \
{ \
  QString line; \
  line.setNum(__LINE__); \
  QString txt = "This seems to be a bug in enGrid"; \
  txt += QString("\n\nfile: ") + __FILE__ + "\nline:" + line + "\n\n"; \
  qWarning()<<txt; \
  abort(); \
};
#else
#define EG_BUG  \
{ \
  QString line; \
  line.setNum(__LINE__); \
  QString txt = "This seems to be a bug in enGrid"; \
  txt += QString("\n\nfile: ") + __FILE__ + "\nline:" + line + "\n\n"; \
  Error err; \
  err.setText(txt); \
  err.setType(Error::ExitOperation); \
  throw err; \
};
#endif


#define EG_ERR_STOP(TXT) \
{ \
  QString line; \
  line.setNum(__LINE__); \
  QString txt = QString(TXT) + "\n\nfile: " + __FILE__ + "\nline:" + line + "\n\n"; \
  Error err; \
  err.setText(txt); \
  err.setType(Error::ExitProgram); \
  throw err; \
};

#define EG_VTKSP(TYPE,VAR) vtkSmartPointer<TYPE> VAR = vtkSmartPointer<TYPE>::New();

/**
 * Get an automatically casted pointer to a cell data array of a grid
 * (vtkUnstructuredGrid).
 * @param TYPE the type to cast to; examples are: vtkLongArray_t, vtkDoubleArray, ...
 * @param VAR  the name for the C++ pointer variable
 * @param GRID a pointer to the vtkUnstructuredGrid
 * @param NAME the sybolic name of the data array; for example "cell_index"
 */
#define EG_VTKDCC(TYPE,VAR,GRID,NAME) \
TYPE *VAR = NULL; \
if (GRID->GetCellData()->GetScalars(NAME)) { \
  VAR = dynamic_cast<TYPE*>(GRID->GetCellData()->GetScalars(NAME)); \
  if (!VAR) { \
    QString msg = "type mismatch ("; \
    int t = GRID->GetCellData()->GetScalars(NAME)->GetDataType(); \
    if (t == VTK_VOID)               msg += "VTK_VOID"; \
    if (t == VTK_BIT)                msg += "VTK_BIT"; \
    if (t == VTK_CHAR)               msg += "VTK_CHAR"; \
    if (t == VTK_SIGNED_CHAR)        msg += "VTK_SIGNED_CHAR"; \
    if (t == VTK_UNSIGNED_CHAR)      msg += "VTK_UNSIGNED_CHAR"; \
    if (t == VTK_SHORT)              msg += "VTK_SHORT"; \
    if (t == VTK_UNSIGNED_SHORT)     msg += "VTK_UNSIGNED_SHORT"; \
    if (t == VTK_INT)                msg += "VTK_INT"; \
    if (t == VTK_UNSIGNED_INT)       msg += "VTK_UNSIGNED_INT"; \
    if (t == VTK_LONG)               msg += "VTK_LONG"; \
    if (t == VTK_UNSIGNED_LONG)      msg += "VTK_UNSIGNED_LONG"; \
    if (t == VTK_FLOAT)              msg += "VTK_FLOAT"; \
    if (t == VTK_DOUBLE)             msg += "VTK_DOUBLE"; \
    if (t == VTK_ID_TYPE)            msg += "VTK_ID_TYPE"; \
    if (t == VTK_STRING)             msg += "VTK_STRING"; \
    if (t == VTK_LONG_LONG)          msg += "VTK_LONG_LONG"; \
    if (t == VTK_UNSIGNED_LONG_LONG) msg += "VTK_UNSIGNED_LONG_LONG"; \
    if (t == VTK___INT64)            msg += "VTK___INT64"; \
    if (t == VTK_UNSIGNED___INT64)   msg += "VTK_UNSIGNED___INT64"; \
    msg += ")"; \
    EG_ERR_RETURN(msg); \
  }; \
} else { \
  QString msg = QString("array '") + NAME + "' does not exist."; \
  EG_ERR_RETURN(msg); \
};

/**
 * Get an automatically casted pointer to a node data array of a grid
 * (vtkUnstructuredGrid).
 * @param TYPE the type to cast to; examples are: vtkLongArray_t, vtkDoubleArray, ...
 * @param VAR  the name for the C++ pointer variable
 * @param GRID a pointer to the vtkUnstructuredGrid
 * @param NAME the sybolic name of the data array; for example "node_index"
 */
#define EG_VTKDCN(TYPE,VAR,GRID,NAME) \
TYPE *VAR = NULL; \
if (GRID->GetPointData()->GetScalars(NAME)) { \
  VAR = dynamic_cast<TYPE*>(GRID->GetPointData()->GetScalars(NAME)); \
  if (!VAR) { \
    QString msg = "type mismatch ("; \
    int t = GRID->GetPointData()->GetScalars(NAME)->GetDataType(); \
    if (t == VTK_VOID)               msg += "VTK_VOID"; \
    if (t == VTK_BIT)                msg += "VTK_BIT"; \
    if (t == VTK_CHAR)               msg += "VTK_CHAR"; \
    if (t == VTK_SIGNED_CHAR)        msg += "VTK_SIGNED_CHAR"; \
    if (t == VTK_UNSIGNED_CHAR)      msg += "VTK_UNSIGNED_CHAR"; \
    if (t == VTK_SHORT)              msg += "VTK_SHORT"; \
    if (t == VTK_UNSIGNED_SHORT)     msg += "VTK_UNSIGNED_SHORT"; \
    if (t == VTK_INT)                msg += "VTK_INT"; \
    if (t == VTK_UNSIGNED_INT)       msg += "VTK_UNSIGNED_INT"; \
    if (t == VTK_LONG)               msg += "VTK_LONG"; \
    if (t == VTK_UNSIGNED_LONG)      msg += "VTK_UNSIGNED_LONG"; \
    if (t == VTK_FLOAT)              msg += "VTK_FLOAT"; \
    if (t == VTK_DOUBLE)             msg += "VTK_DOUBLE"; \
    if (t == VTK_ID_TYPE)            msg += "VTK_ID_TYPE"; \
    if (t == VTK_STRING)             msg += "VTK_STRING"; \
    if (t == VTK_LONG_LONG)          msg += "VTK_LONG_LONG"; \
    if (t == VTK_UNSIGNED_LONG_LONG) msg += "VTK_UNSIGNED_LONG_LONG"; \
    if (t == VTK___INT64)            msg += "VTK___INT64"; \
    if (t == VTK_UNSIGNED___INT64)   msg += "VTK_UNSIGNED___INT64"; \
    msg += ")"; \
    EG_ERR_RETURN(msg); \
  }; \
} else { \
  QString msg = QString("array '") + NAME + "' does not exist."; \
  EG_ERR_RETURN(msg); \
};

#define EG_STDINTERSLOT(OPER) \
OPER *oper = new OPER(); \
(*oper)(); \
oper->del(); \
if(m_grid->GetNumberOfPoints()) updateBoundaryCodes(false); \
updateActors(); \

#define EG_STDSLOT(OPER) \
OPER *oper = new OPER(); \
oper->setLockGui(); \
(*oper)(); \
oper->del(); \
updateActors(); \

#define EG_STDREADERSLOT(OPER) \
OPER *oper = new OPER(); \
(*oper)(); \
oper->del(); \
if(m_grid->GetNumberOfPoints()) updateBoundaryCodes(true); \
updateActors(); \
updateStatusBar(); \
zoomAll();

#define EG_STDCONNECT(OPER) \
connect(ui.action ## OPER, SIGNAL(triggered()), this, SLOT(call ## OPER ()));

#define EG_GETPTS(PTS,CELLID,GRID) \
vtkIdType *PTS; \
vtkIdType  N ## PTS; \
GRID->GetCellPoints(CELLID, N ## PTS, PTS);

inline double sqr(double x) { return x*x; }

inline int factorial_rec(int num)
{
  if (num<=1)
    return 1;
  return factorial_rec(num-1)*num; // recursive call
}

inline int factorial_it(int num)
{
  int result=1;
  for (int i=1; i<=num; ++i)
    result=result*i;
  return result;
}

inline int N_Combinations(int N,int k)
{
  return(factorial_rec(N)/(factorial_rec(N-k)*factorial_rec(k)));
}

inline int N_Permutations(int N,int k)
{
  return(factorial_rec(N)/(factorial_rec(N-k)));
}

#define LINE "========================================================================" << endl;

#define USE(X) X=X

template <class T>
struct SortedPair
{
  T v1, v2;
  SortedPair();
  SortedPair(T v1, T v2);
  bool operator==(const SortedPair<T>& P) const{ return (v1 == P.v1) && (v2 == P.v2); }
  bool operator<(const SortedPair<T>& P) const;
  bool operator>(const SortedPair<T>& P) const { return !operator<(P); }

};

template <class T>
SortedPair<T>::SortedPair(T v1, T v2)
{
  if (v1 > v2) {
    this->v1 = v2;
    this->v2 = v1;
  } else {
    this->v1 = v1;
    this->v2 = v2;
  }
}

template <class T>
SortedPair<T>::SortedPair()
{
  v1 = 0;
  v2 = 0;
}

template <class T>
bool SortedPair<T>::operator<(const SortedPair<T>& P) const
{
  if (v1 < P.v1) {
    return true;
  } else if (v1 == P.v1) {
    return v2 < P.v2;
  }
  return false;
}

template <class T>
inline uint qHash(const SortedPair<T> &P)
{
  uint h1 = qHash(P.v1);
  uint h2 = qHash(P.v2);
  return ((h1 << 16) | (h1 >> 16)) ^ h2;
}


#endif
