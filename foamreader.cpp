#include "foamreader.h"

FoamReader::FoamReader()
{
  setFormat("Foam boundary files(boundary)");
  setExtension("");
};

void FoamReader::operate()
{
  try {
    readInputFileName();
    if (isValid()) {
    };
  } catch (Error err) {
    err.display();
  };
};
