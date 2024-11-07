#! /usr/bin/python

from glas_algorithm import *

copyright_notice = "Karl Meerbergen (2012)"

x_arg = templated_argument( "X", "x" )

x_arg_const = templated_argument( "X", "x" )
x_arg_const.mutable = False

storage_size = free_concept_function( 'glas/toolbox/tensor/concept/storage_size', ['glas'] )
storage_size.set_arguments( [x_arg_const] )
write_file( storage_size, "Karl Meerbergen (2011)" )
