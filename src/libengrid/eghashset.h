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
#ifndef SORTEDHASHSET_H
#define SORTEDHASHSET_H

template <class T> class EgHashSet;

#include <QList>

template <class T>
class EgHashSet
{

private: // data types

  struct entry_t
  {
    T   value;
    int index;
  };

private: // attributes

  QVector<QList<entry_t> > m_Buckets;
  int                      m_NumEntries;
  QVector<int>             m_Old2NewIndex;

public: // methods

  EgHashSet();
  EgHashSet(int size);

  void resize(int size);
  void clear();
  int  insert(const T &item);
  void getQVector(QVector<T> &qv);
  void updateIndexMap();
  int  newIndex(int old_index) { return m_Old2NewIndex[old_index]; }

};

template <class T>
EgHashSet<T>::EgHashSet()
{
  clear();
}

template <class T>
EgHashSet<T>::EgHashSet(int size)
{
  clear();
  m_Buckets.resize(size);
}

template <class T>
void EgHashSet<T>::resize(int size)
{
  clear();
  m_Buckets.resize(size);
}

template <class T>
void EgHashSet<T>::clear()
{
  m_Buckets.clear();
  m_NumEntries = 0;
}

template <class T>
int EgHashSet<T>::insert(const T &item)
{
  QList<entry_t> &list = m_Buckets[item.hash()];
  int idx = 0;
  while (idx < list.size()) {
    if (!(item < list[idx].value)) {
      break;
    }
    ++ idx;
  }
  if (idx < list.size()) {
    if (item == list[idx].value) {
      return list[idx].index;
    }
  }
  entry_t entry;
  entry.value = item;
  entry.index = m_NumEntries;
  list.insert(idx, entry);
  ++m_NumEntries;
  return entry.index;
}

template <class T>
void EgHashSet<T>::updateIndexMap()
{
  int idx = 0;
  m_Old2NewIndex.clear();
  m_Old2NewIndex.fill(-1, m_NumEntries);
  for (int i = 0; i < m_Buckets.size(); ++i) {
    foreach (entry_t entry, m_Buckets[i]) {
      m_Old2NewIndex[entry.index] = idx;
      ++idx;
    }
  }
}

template <class T>
void EgHashSet<T>::getQVector(QVector<T> &qv)
{
  updateIndexMap();
  qv.clear();
  qv.resize(m_NumEntries);
  QVector<bool> entry_set(m_NumEntries, false);
  for (int i = 0; i < m_Buckets.size(); ++i) {
    foreach (entry_t entry, m_Buckets[i]) {
      qv[entry.index] = entry.value;
      entry_set[entry.index] = true;
    }
  }
  foreach (bool is_set, entry_set) {
    if (!is_set) {
      EG_BUG;
    }
  }
}

#endif // SORTEDHASHSET_H
