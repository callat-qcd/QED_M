#include <stdlib.h>
#include "lattice.h"

// This is code for performion the block-averaging of a field configuration to a field configuration
// half its size; there are two versions

#ifndef INCLUDED_BLOCK
#define INCLUDED_BLOCK
void BlockPAC_CS(Lattice &, Lattice &); // this version does exactly what PAC-CS does
void BlockPAC_CS(Lattice &, Lattice &, int); // this version blocks over p-cells for p < q
                                             // (q being the last argument)
#endif
