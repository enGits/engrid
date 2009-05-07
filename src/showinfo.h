//
// C++ Interface: showinfo
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SHOWINFO_H
#define SHOWINFO_H

#include <operation.h>

/**
@author Mike Taverne <mtaverne@engits.com>
*/
class ShowInfo : public Operation
{
private:
  bool CellInfo;
  vtkIdType PickedPoint;
  vtkIdType PickedCell;
  
public:
  ShowInfo(bool b, vtkIdType P, vtkIdType C);
  
//     ~ShowInfo();
  virtual void operate();
  
};

#endif
