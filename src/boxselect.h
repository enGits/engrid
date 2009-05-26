//
// C++ Interface: boxselect
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BOXSELECT_H
#define BOXSELECT_H

#include <operation.h>

/**
	@author Mike Taverne <mtaverne@engits.com>
*/
class BoxSelect : public Operation
{
public:
    BoxSelect();

    ~BoxSelect();

protected: // methods
  virtual void operate();
};

#endif
