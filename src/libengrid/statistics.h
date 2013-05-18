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
#ifndef STATISTICS_H
#define STATISTICS_H

namespace Statistics {

template <typename C>
typename C::value_type meanValue(const C &container)
{
  typename C::value_type mean = typename C::value_type();
  for (typename C::const_iterator i = container.begin(); i != container.end(); ++i) {
    mean += *i;
  }
  mean /= container.size();
  return mean;
}

template <typename C>
typename C::value_type standardDeviation(const C &container, typename C::value_type mean)
{
  typename C::value_type sigma = typename C::value_type();
  for (typename C::const_iterator i = container.begin(); i != container.end(); ++i) {
    sigma += ((*i) - mean)*((*i) - mean);
  }
  sigma /= container.size();
  sigma = sqrt(sigma);
  return sigma;
}

template <typename C>
typename C::value_type standardDeviation(const C &container)
{
  typename C::value_type mean = meanValue(container);
  return standardDeviation(container, mean);
}


}

#endif // STATISTICS_H
