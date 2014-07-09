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
template <class T>
struct mv_p
{
  static T apply(const T &a, const T &b) { return a+b; };
};

template <class T>
struct mv_m
{
  static T apply(const T &a, const T &b) { return a-b; };
};

template <class T>
struct mv_ml
{
  static T apply(const T &a, const T &b) { return a*b; };
};

template <class L, class O, class R>
struct ParseNode
{
  typedef typename R::value_type value_type;
  const L &l;
  const R &r;
  ParseNode(const L &a, const R &b) : l(a), r(b) {};
  value_type operator[](const uint_t &i) const { return O::apply(l[i], r[i]); };
  uint_t size() const { return r.size(); };
  value_type abs() const;
  value_type abs2() const;
};

template <class O, class R>
struct ParseNode<double, O, R>
{
  typedef typename R::value_type value_type;
  const double l;
  const R &r;
  ParseNode(const double a, const R &b) : l(a), r(b) {};
  value_type operator[](const uint_t &i) const { return O::apply(l, r[i]); };
  uint_t size() const { return r.size(); };
  value_type abs() const;
  value_type abs2() const;
};
