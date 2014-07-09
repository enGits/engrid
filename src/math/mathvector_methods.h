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

// MathVector<V>(const MathVector<V>&)
// -----------------------------------
template <class V>
inline MathVector<V>::MathVector(const MathVector<V>& vec) : V()
{
  for (uint_t i = 0; i < this->size(); ++i) (*this)[i] = vec[i]; 
};

// MathVector<V>(const ParseNode<L,O,R>&)
// --------------------------------------
template <class V>
template <class L, class O, class R>
inline MathVector<V>::MathVector(const ParseNode<L,O,R> &expr) : V()
{
  for (uint_t i = 0; i < this->size(); ++i) (*this)[i] = expr[i]; 
};

// MathVector<V>(const value_type *v)
// ----------------------------------
template <class V>
inline MathVector<V>::MathVector(const value_type *v) : V()
{
  for (uint_t i = 0; i < this->size(); ++i) (*this)[i] = v[i]; 
};

// MathVector<V>(const value_type v1, const value_type v2, const value_type v3, const value_type v4)
// ----------------------------------------------------------------------------
template <class V>
inline MathVector<V>::MathVector(const value_type v1,
                                 const value_type v2,
                                 const value_type v3,
                                 const value_type v4) : V()
{
  (*this)[0] = v1;
  (*this)[1] = v2;
  (*this)[2] = v3;
  (*this)[3] = v4;
};

// MathVector<V>(const value_type v1, const value_type v2, const value_type v3)
// ----------------------------------------------------------------------------
template <class V>
inline MathVector<V>::MathVector(const value_type v1,
                                 const value_type v2,
                                 const value_type v3) : V()
{
  (*this)[0] = v1;
  (*this)[1] = v2;
  (*this)[2] = v3;
};

// MathVector<V>(const value_type v1, const value_type v2)
// ----------------------------------------------------------------------------
template <class V>
inline MathVector<V>::MathVector(const value_type v1,
                                 const value_type v2) : V()
{
  (*this)[0] = v1;
  (*this)[1] = v2;
};

// cross product
// -------------
template <class V>
inline MathVector<V> MathVector<V>::cross(const MathVector<V> &vec) const
{
  MathVector<V> new_vec;
  new_vec[0] = (*this)[1]*vec[2] - (*this)[2]*vec[1];
  new_vec[1] = (*this)[2]*vec[0] - (*this)[0]*vec[2];
  new_vec[2] = (*this)[0]*vec[1] - (*this)[1]*vec[0];
  return new_vec;
};

// absolute value (||v||_2 in this case)
// -------------------------------------
template <class V>
inline typename MathVector<V>::scalar_t MathVector<V>::abs() const
{
  scalar_t l = 0;
  for (uint_t i = 0; i < this->size(); ++i) l += (*this)[i]*(*this)[i];
  return sqrt(l);
};

// absolute value squared ((||v||_2)^2 in this case)
// -------------------------------------------------
template <class V>
inline typename MathVector<V>::scalar_t MathVector<V>::abs2() const
{
  scalar_t l = 0;
  for (uint_t i = 0; i < this->size(); ++i) l += (*this)[i]*(*this)[i];
  return l;
};

// absolute value (||v||_2 in this case)
// -------------------------------------
template <class L, class O, class R>
inline typename ParseNode<L,O,R>::value_type ParseNode<L,O,R>::abs() const
{
  value_type abs_value = 0;
  for (uint_t i = 0; i < this->size(); ++i) abs_value += (*this)[i]*(*this)[i];
  return sqrt(abs_value);
};

template <class O, class R>
inline typename ParseNode<double,O,R>::value_type ParseNode<double,O,R>::abs() const
{
  value_type abs_value = 0;
  for (uint_t i = 0; i < this->size(); ++i) abs_value += (*this)[i]*(*this)[i];
  return sqrt(abs_value);
};

// absolute value squared ((||v||_2)^2 in this case)
// -------------------------------------------------
template <class L, class O, class R>
inline typename ParseNode<L,O,R>::value_type ParseNode<L,O,R>::abs2() const
{
  value_type abs_value = 0;
  for (uint_t i = 0; i < this->size(); ++i) abs_value += (*this)[i]*(*this)[i];
  return abs_value;
};

template <class O, class R>
inline typename ParseNode<double,O,R>::value_type ParseNode<double,O,R>::abs2() const
{
  value_type abs_value = 0;
  for (uint_t i = 0; i < this->size(); ++i) abs_value += (*this)[i]*(*this)[i];
  return abs_value;
};

// norm the vector to a length of 1 (according to ||v||_2)
// -------------------------------------------------------
template <class V>
inline MathVector<V> MathVector<V>::normalise()
{
  scalar_t l = abs();
  for (uint_t i = 0; i < this->size(); ++i) (*this)[i] /= l;
  return(*this);
};

// create a C-style array
// ----------------------
template <class V>
typename MathVector<V>::scalar_t* MathVector<V>::c_array() const 
{
  scalar_t *t = new scalar_t[this->size()];
  for (uint_t i = 0; i < this->size(); ++i) t[i] = (*this)[i];
  return t;
};

