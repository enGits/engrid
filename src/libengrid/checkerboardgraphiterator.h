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

#ifndef CHECKERBOARDGRAPHITERATOR_H
#define CHECKERBOARDGRAPHITERATOR_H

#include <vector>

using namespace std;

template <typename TGraph>
class CheckerBoardGraphIterator
{

  typedef typename TGraph::index_type index_t;

  TGraph*      m_Graph;
  index_t      m_CurrentIndex;
  index_t      m_DoneCount;
  vector<bool> m_Available;
  vector<bool> m_Done;
  bool         m_UpdateRequired;

public:

  CheckerBoardGraphIterator() { m_Graph = NULL; }
  CheckerBoardGraphIterator(TGraph* graph) { m_Graph = graph; }
  void setGraph(TGraph* graph) { m_Graph = graph; }
  bool updateRequired();
  void operator=(index_t i);
  bool operator==(index_t i) { return m_CurrentIndex == i; }
  bool operator<(index_t i) { return m_CurrentIndex < i; }
  index_t operator++();
  index_t operator*() { return m_CurrentIndex; }

};

template <typename TGraph>
void CheckerBoardGraphIterator<TGraph>::operator=(index_t i)
{
  m_CurrentIndex = i;
  m_Available.resize(m_Graph->size(), true);
  m_Done.resize(m_Graph->size(), false);
  m_DoneCount = 0;
  m_UpdateRequired = false;
}

template <typename TGraph>
bool CheckerBoardGraphIterator<TGraph>::updateRequired()
{
  bool update_required = m_UpdateRequired;
  m_UpdateRequired = false;
  return update_required;
}

template <typename TGraph>
typename CheckerBoardGraphIterator<TGraph>::index_t CheckerBoardGraphIterator<TGraph>::operator++()
{
  if (m_Graph->size() > 0) {
    m_Done[m_CurrentIndex] = true;
    ++m_DoneCount;
    for (int i = 0; i < m_Graph->getNumLinks(m_CurrentIndex); ++i) {
      m_Available[m_Graph->getLink(m_CurrentIndex,i)] = false;
    }
    do {
      ++m_CurrentIndex;
      if (m_CurrentIndex >= m_Graph->size()) {
        if (m_DoneCount == m_Graph->size()) {
          break;
        }
        for (size_t i = 0; i < m_Available.size(); ++i) {
          m_Available[i] = !m_Done[i];
        }
        m_CurrentIndex = 0;
        m_UpdateRequired = true;
      }
    } while (!m_Available[m_CurrentIndex]);
  }
  return m_CurrentIndex;
}



#endif // CHECKERBOARDGRAPHITERATOR_H
