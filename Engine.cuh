#pragma once
#include "Cell.h"
template<unsigned d,unsigned q>
sogol::Cell<d, q>*  initCu(sogol::Cell<d, q> *c,int size);
template<unsigned d, unsigned q>
void runCu(sogol::Cell<d, q>* c, int size);
