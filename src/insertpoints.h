//
// C++ Interface: insertpoints
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef INSERTPOINTS_H
#define INSERTPOINTS_H

#include <operation.h>

/**
	@author Mike Taverne <mtaverne@engits.com>
*/
class InsertPoints : public Operation
{
public:
  InsertPoints();

  ~InsertPoints();
  
  void operate();
  
};

#endif
