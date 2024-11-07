//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_options_hpp
#define glas2_iterative_krylov_options_hpp

namespace glas2 { namespace iterative {

  struct options {
    double        absolute_tolerance_ ;
    double        relative_tolerance_ ;
    unsigned int  max_mat_vec_ ;

    options()
    : absolute_tolerance_(0.0)
    , relative_tolerance_(0)
    , max_mat_vec_(0)
    {}

  } ;
} }

#endif
