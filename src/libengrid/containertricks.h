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
#ifndef containertricks_h
#define containertricks_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdarg.h>

// vlinit 
// ======

template <class C>
class vlinit_t
{
  C *c;
  
public:
  vlinit_t(C &a_c);
  vlinit_t(const vlinit_t<C> &ci);
  vlinit_t<C> add(typename C::value_type v);
  vlinit_t<C> operator=(typename C::value_type v) { return add(v); };
  vlinit_t<C> operator,(typename C::value_type v) { return add(v); };
};
  
template <class C>
inline vlinit_t<C>::vlinit_t(const vlinit_t<C> &ci) 
{ 
  c = ci.c;
};

template <class C>
inline vlinit_t<C>::vlinit_t(C &a_c) 
{ 
  c = &a_c;
};
  
template <class C>
inline vlinit_t<C> vlinit_t<C>::add(typename C::value_type v) 
{ 
  c->push_back(v);
  return *this;
};

template <class C> inline vlinit_t<C> vlinit(C &c) { return vlinit_t<C>(c); };


// clinit
// ======

template <class C>
class clinit_t
{
  C *c;
  typename C::iterator i;
  
public:
  clinit_t(C &a_c);
  clinit_t(const clinit_t<C> &ci);
  clinit_t<C> add(typename C::value_type v);
  clinit_t<C> operator=(typename C::value_type v) { return add(v); };
  clinit_t<C> operator,(typename C::value_type v) { return add(v); };
};
  
template <class C>
inline clinit_t<C>::clinit_t(const clinit_t<C> &ci) 
{ 
  c = ci.c;
  i = ci.i;
};

template <class C>
inline clinit_t<C>::clinit_t(C &a_c) 
{ 
  c = &a_c;
  i = c->begin();
};
  
template <class C>
inline clinit_t<C> clinit_t<C>::add(typename C::value_type v) 
{ 
  if (i == c->end()) {
    cerr << "array bounds exceeded" << endl;
    exit(EXIT_FAILURE);
  };
  *i = v;
  i++;
  return *this;
};

template <class C> inline clinit_t<C> clinit(C &c) { return clinit_t<C>(c); };

template <class C>
inline void clinit(C &c, typename C::value_type v, ...)
{
  if (c.size() == 0) return;
  typename C::iterator i = c.begin();
  *i = v;
  ++i;
  va_list vl;
  va_start(vl,v);
  cout << v << ' ';
  while (v = va_arg(vl,typename C::value_type)) {
    if (i == c.end()) {
      cerr << "array bounds exceeded" << endl;
      exit(EXIT_FAILURE);
    };
    *i = v;
    cout << v << ' ';
    ++i;
  };
  cout << endl;
};


// output tools
// ============

template <class C>
inline void simple_print(const C &c, ostream &s)
{
  typename C::const_iterator i = c.begin();
  s << '[';
  while (i != c.end()) {
    s << *i;
    i++;
    if (i != c.end()) s << ", ";
  };
  s << ']';
};


inline void print_table(vector<vector<double> > f, ostream &s)
{
  size_t Nj = f[0].size();
  for (size_t j = 0; j < Nj; j++) {
    for (size_t i = 0; i < f.size(); i++) {
      s.setf(iostream::scientific, iostream::floatfield);
      s.precision(4);
      s.width(11);
      s << f[i][j] << ' ';
    };
    s << endl;
  };
};

inline void print_table(vector<vector<double> > f, string file_name)
{
  file_name += ".dat";
  ofstream s(file_name.c_str());
  print_table(f, s);
};

#endif
