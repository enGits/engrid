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
#include <typeinfo>

// asssignements
// =============

// MathVector<V> = const MathVector<V>
// -----------------------------------
template <class V>
inline void MathVector<V>::operator= (const MathVector<V> &v)
{
  for (uint_t i = 0; i < this->size(); ++i) {
    (*this)[i] = v[i]; 
  };
};

// MathVector<V> = ParseNode<L,O,R>
// --------------------------------
template <class V>
template <class L, class O, class R>
inline void MathVector<V>::operator= (const ParseNode<L,O,R> &expr)
{
  for (uint_t i = 0; i < this->size(); ++i) {
    (*this)[i] = expr[i];
  };
};

// MathVector<V> = vector<V::value_type>
// -------------------------------------
template <class V>
inline void MathVector<V>::operator= (const vector<typename V::value_type> &v)
{
  for (uint_t i = 0; i < this->size(); ++i) {
    (*this)[i] = v[i];
  };
};




// sums
// ====

// MathVector<V> + MathVector<V>
// -----------------------------
template <class V>
inline ParseNode<MathVector<V>, mv_p<typename V::value_type>, MathVector<V> >
operator+ (const MathVector<V> &a,
	   const MathVector<V> &b)
{
  return ParseNode<MathVector<V>, mv_p<typename V::value_type>, MathVector<V> > (a, b);
};

// MathVector<V> + ParseNode<L,O,R>
// --------------------------------
template <class V, class L, class O, class R>
inline ParseNode<MathVector<V>, mv_p<typename ParseNode<L,O,R>::value_type>, ParseNode<L,O,R> >
operator+ (const MathVector<V>    &a,
	   const ParseNode<L,O,R> &b)
{
  return ParseNode<MathVector<V>, mv_p<typename ParseNode<L,O,R>::value_type>, ParseNode<L,O,R> > (a, b);
};

// ParseNode<L,O,R> + MathVector<V>
// --------------------------------
template <class V, class L, class O, class R>
inline ParseNode<ParseNode<L,O,R>, mv_p<typename V::value_type>, MathVector<V> >
operator+ (const ParseNode<L,O,R> &a,
	   const MathVector<V>    &b)
{
  return ParseNode<ParseNode<L,O,R>, mv_p<typename V::value_type>, MathVector<V> > (a, b);
};

// ParseNode<L1,O1,R1> + ParseNode<L2,O2,R2>
// -----------------------------------------
template <class L1, class O1, class R1, class L2, class O2, class R2>
inline ParseNode<ParseNode<L1,O1,R1>, mv_p<typename ParseNode<L2,O2,R2>::value_type>, ParseNode<L2,O2,R2> >
operator+ (const ParseNode<L1,O1,R1> &a,
	   const ParseNode<L2,O2,R2> &b)
{
  return ParseNode<ParseNode<L1,O1,R1>, mv_p<typename ParseNode<L2,O2,R2>::value_type>, ParseNode<L2,O2,R2> > (a, b);
};

// MathVector<V> -= MathVector<V>
// ------------------------------
template <class V>
inline void MathVector<V>::operator-= (const MathVector<V> &v) 
{
  for (uint_t i = 0; i < this->size(); ++i) (*this)[i] -= v[i]; 
};

// MathVector<V> += MathVector<V>
// ------------------------------
template <class V>
inline void MathVector<V>::operator+= (const MathVector<V> &v) 
{
  for (uint_t i = 0; i < this->size(); ++i) (*this)[i] += v[i]; 
};




// differences
// ===========

// MathVector<V> - MathVector<V>
// -----------------------------
template <class V>
inline ParseNode<MathVector<V>, mv_m<typename V::value_type>, MathVector<V> >
operator- (const MathVector<V> &a,
	   const MathVector<V> &b)
{
  return ParseNode<MathVector<V>, mv_m<typename V::value_type>, MathVector<V> > (a, b);
};

// MathVector<V> - ParseNode<L,O,R>
// --------------------------------
template <class V, class L, class O, class R>
inline ParseNode<MathVector<V>, mv_m<typename ParseNode<L,O,R>::value_type>, ParseNode<L,O,R> >
operator- (const MathVector<V>    &a,
	   const ParseNode<L,O,R> &b)
{
  return ParseNode<MathVector<V>, mv_m<typename ParseNode<L,O,R>::value_type>, ParseNode<L,O,R> > (a, b);
};

// ParseNode<L,O,R> - MathVector<V>
// --------------------------------
template <class V, class L, class O, class R>
inline ParseNode<ParseNode<L,O,R>, mv_m<typename V::value_type>, MathVector<V> >
operator- (const ParseNode<L,O,R> &a,
	   const MathVector<V>    &b)
{
  return ParseNode<ParseNode<L,O,R>, mv_m<typename V::value_type>, MathVector<V> > (a, b);
};

// ParseNode<L1,O1,R1> - ParseNode<L2,O2,R2>
// -----------------------------------------
template <class L1, class O1, class R1, class L2, class O2, class R2>
inline ParseNode<ParseNode<L1,O1,R1>, mv_m<typename ParseNode<L2,O2,R2>::value_type>, ParseNode<L2,O2,R2> >
operator- (const ParseNode<L1,O1,R1> &a,
	   const ParseNode<L2,O2,R2> &b)
{
  return ParseNode<ParseNode<L1,O1,R1>, mv_m<typename ParseNode<L2,O2,R2>::value_type>, ParseNode<L2,O2,R2> > (a, b);
};




// S-products
// ==========

// double * ParseNode<L,O,R>
// -------------------------
template <class L, class O, class R>
inline ParseNode<double, mv_ml<typename ParseNode<L,O,R>::value_type>, ParseNode<L,O,R> >
operator* (const double &a, 
	   const ParseNode<L,O,R> &b)
{
  return ParseNode<double, mv_ml<typename ParseNode<L,O,R>::value_type>, ParseNode<L,O,R> > (a, b);
};

// double * MathVector<V>
// ----------------------
template <class V>
inline ParseNode<double, mv_ml<typename V::value_type>, MathVector<V> >
operator* (const double        &a,
	   const MathVector<V> &b)
{
  return ParseNode<double, mv_ml<typename V::value_type>, MathVector<V> > (a, b);
};

// MathVector<V> *= value_type
// ---------------------------
template <class V>
inline void MathVector<V>::operator*= (const typename V::value_type s) 
{
  for (uint_t i = 0; i < this->size(); ++i) (*this)[i] *=s;
};




// scalar-products
// ===============

// MathVector<V> * MathVector<V>
// -----------------------------
template <class V>
inline typename V::value_type
operator* (const MathVector<V> &a,
	   const MathVector<V> &b)
{
  typename V::value_type s = 0;
  for (uint_t i = 0; i < a.size(); ++i) s += a[i]*b[i];
  return s;
};

// MathVector<V> * ParseNode<L,O,R>
// --------------------------------
template <class V, class L, class O, class R>
inline typename V::value_type
operator* (const MathVector<V>    &a,
	   const ParseNode<L,O,R> &b)
{
  typename ParseNode<L,O,R>::value_type s = 0;
  for (uint_t i = 0; i < a.size(); ++i) s += a[i]*b[i];
  return s;
};

// ParseNode<L,O,R> * MathVector<V>
// --------------------------------
template <class V, class L, class O, class R>
inline typename V::value_type
operator* (const ParseNode<L,O,R> &a,
	   const MathVector<V>    &b)
{
  typename V::value_type s = 0;
  for (uint_t i = 0; i < a.size(); ++i) s += a[i]*b[i];
  return s;
};

// ParseNode<L1,O1,R1> * ParseNode<L2,O2,R2>
// -----------------------------------------
template <class L1, class O1, class R1, class L2, class O2, class R2>
inline typename ParseNode<L2,O2,R2>::value_type
operator* (const ParseNode<L1,O1,R1> &a,
	   const ParseNode<L2,O2,R2> &b)
{
  typename ParseNode<L2,O2,R2>::value_type s = 0;
  for (uint_t i = 0; i < a.size(); ++i) s += a[i]*b[i];
  return s;
};




// stream operators
// ================

// ostream << MathVector<V>
// ------------------------
template <class V>
ostream& operator<<(ostream &s, const MathVector<V> &vec)
{
  //s << '(' << &vec << ',' << vec.size() << ')' << endl;
  s << '[';
  for (uint_t i = 0; i < vec.size(); ++i) {
    s << vec[i];
    if (i != vec.size() - 1) s << ',';
  };
  s << ']';
  return s;
};

// ostream << ParseNode<L,O,R>
// ---------------------------
template <class L, class O, class R>
ostream& operator<<(ostream &s, const ParseNode<L,O,R> &expr)
{
  s << '[';
  for (uint_t i = 0; i < expr.size(); ++i) {
    s << expr[i];
    if (i != expr.size() - 1) s << ',';
  };
  s << ']';
  return s;
};




// STL iterator operators
// ======================

// iterator - uint_t
// -----------------
template <class V>
inline typename MathVector<V>::iterator
operator- (const typename MathVector<V>::iterator &iter, 
	   uint_t                                  n)
{
  typename MathVector<V>::iterator new_iter = iter;
  new_iter -= n;
  return new_iter;
};

// iterator - uint_t
// -----------------
template <class V>
inline typename MathVector<V>::iterator
operator+ (const typename MathVector<V>::iterator &iter, 
	   uint_t                                    n)
{
  typename MathVector<V>::iterator new_iter = iter;
  new_iter += n;
  return new_iter;
};

// uint_t - iterator
// -----------------
template <class V>
inline typename MathVector<V>::iterator
operator+ (uint_t                                    n, 
	   const typename MathVector<V>::iterator   &iter)
{
  typename MathVector<V>::iterator new_iter = iter;
  new_iter += n;
  return new_iter;
};

// iterator == const_iterator
// --------------------------
template <class V>
inline bool 
MathVector<V>::iterator::operator== (const typename MathVector<V>::const_iterator &iter) const
{ 
  return (iter.vec == vec) && (iter.i == i); 
};

// const_iterator == iterator
// --------------------------
template <class V>
inline bool 
MathVector<V>::const_iterator::operator== (const typename MathVector<V>::iterator &iter) const
{ 
  return (iter.vec == vec) && (iter.i == i); 
};



