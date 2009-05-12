//
// C++ Interface: removepoints
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef REMOVEPOINTS_H
#define REMOVEPOINTS_H

#include <operation.h>

/**
	@author Mike Taverne <mtaverne@engits.com>
*/
class RemovePoints : public Operation
{
public:
  RemovePoints();

  ~RemovePoints();
  
  void operate();
  
};

#endif
