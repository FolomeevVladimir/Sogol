#pragma once
#include "Cell.h"
#include "Physics.h"
#include "Stencil.h"
//#include "Vectorxd.h"
template<typename dq>
sogol::Cell<dq::d, dq::q>*  initCu(sogol::Cell<dq::d, dq::q> *c,int size);
template<typename dq,typename p>
void runCu(sogol::Cell<dq::d, dq::q>* c, int size, sogol::Vectorxd<dq::d, int> mask);
