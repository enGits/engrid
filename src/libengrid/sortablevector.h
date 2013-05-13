// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#ifndef sortablevector_H
#define sortablevector_H

class SortableVector;

#include <QVector>

template <class T>
class SortableVector : public QVector<T>
{
  
public: // methods
  
  SortableVector() : QVector<T>() {};
  SortableVector(int size) : QVector<T>(size) {};
  SortableVector(int size, const T &value) QVector<T>(size, value) {};
  SortableVector(const SortableVector<T> &other) : QVector<T>(other) {};
  
  bool operator<(const SortableVector<T> &V) const;
  
};

template <class T>
bool SortableVector<T>::operator<(const SortableVector<T> &V) const
{
  EG_BUG;
};

#endif
