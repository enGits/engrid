// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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

#include "timer.h"

Timer::Timer()
{
  reset(10);
  print = false;
}

Timer::Timer(int secs)
{
  reset(secs);
}

void Timer::reset(int secs)
{
  m_Secs = secs;
  m_LastTimeout = QTime::currentTime();
}

bool Timer::timeout()
{
  QTime now = QTime::currentTime();
  bool timeout = m_LastTimeout.secsTo(now) >= m_Secs;
  if (timeout) {
    m_LastTimeout = now;
  }
  return timeout;
}


