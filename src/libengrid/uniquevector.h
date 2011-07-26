// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
#ifndef uniquevector_H
#define uniquevector_H

class uniqueVector;

#include <QVector>

template <class T>
class UniqueVector : public QVector<T>
{
  
public: // methods
  
  UniqueVector() : QVector<T>() {};
  UniqueVector(int size) : QVector<T>(size) {};
  UniqueVector(int size, const T &value) : QVector<T>(size, value) {};
  UniqueVector(const UniqueVector<T> &other) : QVector<T>(other) {};
  
  bool operator==(const UniqueVector<T> &V) const;
  
};

template <class T>
bool UniqueVector<T>::operator==(const UniqueVector<T> &V) const
{
  if (QVector<T>::size() != V.size()) return false;
  QVector<bool> used(QVector<T>::size(),false);
  int N = 0;
  for (int i = 0; i < QVector<T>::size(); ++i) {
    for (int j = 0; j < QVector<T>::size(); ++j) {
      if ((V[j] == this->operator[](i)) && !used[j]) {
        ++N;
        used[j] = true;
        break;
      };
    };
  };
  return N == QVector<T>::size();
};

#endif
