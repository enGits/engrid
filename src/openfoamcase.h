#ifndef OPENFOAMCASE_H
#define OPENFOAMCASE_H

#include "iooperation.h"

class OpenFOAMcase : public IOOperation
{
public:
    OpenFOAMcase();
    virtual void operate();
};

#endif
