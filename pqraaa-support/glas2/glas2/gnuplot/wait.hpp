//  (C) Copyright Karl Meerbergen 2011. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_gnuplot_wait_hpp
#define glas2_gnuplot_wait_hpp

#include <time.h>

namespace glas2 { namespace gnuplot {

void wait ( int milleseconds )
{
    clock_t endwait;
    endwait = clock () + milleseconds * CLOCKS_PER_SEC/1000 ;
    while (clock() < endwait) {}
}


} }

#endif
