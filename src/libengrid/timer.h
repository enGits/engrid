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

#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <qdatetime.h>
#include <typeinfo>

class Timer
{

private: // attributes

  bool print;

protected: // attributes

  QTime m_LastTimeout;
  int m_Secs;


public: // attributes

  static const char endl = '\n';


public:

  Timer();
  Timer(int secs);
  void reset(int secs = 0);
  bool timeout();
  bool operator()() { return timeout(); }

  template <class T> Timer& operator<<(T t);

};

template <class T>
Timer& Timer::operator<<(T t)
{
  if (print) {
    bool endl_received = false;
    if (typeid(T) == typeid(const char)) {
      endl_received = true;
    }
    if (print) {
      if (endl_received) {
        std::cout << std::endl;
        print = false;
      } else {
        std::cout << t;
      }
    } else {
      if (endl_received) {
        if (timeout()) {
          print = true;
        }
      }
    }
  }
  return (*this);
};

#endif // TIMER_H
