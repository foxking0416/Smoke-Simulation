// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {9, 10, 1};
#else
const int theDim[3] = {2, 2, 1};//12, 12, 4
#endif

const double theCellSize = 0.5;

// TODO: add constants for dispersed bubble simulation