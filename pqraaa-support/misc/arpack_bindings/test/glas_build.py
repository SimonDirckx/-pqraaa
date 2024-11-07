build = 'exe'
variant = ['debug','release']

arpack.check_modif = True
arpack_bindings.check_modif = True

dependencies += arpack
dependencies += arpack_bindings
dependencies += glas2
dependencies += lapack
dependencies += boost_numeric_bindings
