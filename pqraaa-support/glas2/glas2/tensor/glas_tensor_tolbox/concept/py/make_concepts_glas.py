#! /usr/bin/python

from glas_concept import *

copyright_notice = "Karl Meerbergen (2011)"
namespace = ["glas"]

collection = concept( 'glas/concept/collection', namespace, [] )
expression = concept( 'glas/concept/expresson', namespace, [] )

# DenseTensorExpressions

dense_tensor_collection = concept( 'glas/toolbox/tensor/concept/dense_tensor_collection', namespace, [] )
concepts.append( dense_tensor_collection )

dense_tensor_expression = concept( 'glas/toolbox/tensor/concept/dense_tensor_expression', namespace, [] )
concepts.append( dense_tensor_expression )

# TensorExpressions

tensor_collection = concept( 'glas/toolbox/tensor/concept/tensor_collection', namespace, [collection] )
tensor_collection.children = [dense_tensor_collection]
concepts.append( tensor_collection )

tensor_expression = concept( 'glas/toolbox/tensor/concept/tensor_expression', namespace, [expression] )
tensor_expression.children = [dense_tensor_expression]
concepts.append( tensor_expression )

handle_dir( concepts, copyright_notice )
