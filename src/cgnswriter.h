#ifndef CGNSWRITER_H
#define CGNSWRITER_H

class CgnsWriter;

#include "iooperation.h"

#ifdef CGNS_SUPPORT
#include "cgnslib.h"
#endif

/**
 * Writer for CGNS files.
 */
class CgnsWriter : public IOOperation
{

protected: // attributes

  int fn;
  int B;
  int Z;
  QVector<int> eg2cgns;

protected: // methods

  void writeGrid();
  void writeBcs();
  virtual void operate();

public: // methods

  CgnsWriter();

};

#endif // CGNSWRITER_H
