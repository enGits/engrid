#include "openfoamcase.h"

#include "filetemplate.h"

#include <iostream>
using namespace std;

OpenFOAMcase::OpenFOAMcase()
    : IOOperation()
{
}

void OpenFOAMcase::operate()
{
  cout << "OpenFOAMcase::operate()" << endl;
}
