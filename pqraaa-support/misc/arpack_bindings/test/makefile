CXXFLAGS =                 -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software/ -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0
LDFLAGS =        
LIBS =  -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//ARPACK/ -larpack -L/usr/local/gfortran/lib -lgfortran -L/usr/lib -llapack -L/usr/lib -lblas -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/lib -larpack_bindings  -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0/libs -lboost_1_59_0

NAUPD_STANDARD_DEPS = \
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/naupd.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/data.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/neupd.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/standard.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/no_op.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/copy_op.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/random_vector_op.hpp

naupd_standard_d.o : naupd_standard.cpp $(NAUPD_STANDARD_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o naupd_standard_d.o naupd_standard.cpp

naupd_standard.o : naupd_standard.cpp $(NAUPD_STANDARD_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o naupd_standard.o naupd_standard.cpp


naupd_standard_d: naupd_standard_d.o
	@c++ -o naupd_standard_d $(LDFLAGS) -g naupd_standard_d.o $(LIBS) -lstdc++ -lpthread

naupd_standard: naupd_standard.o
	@c++ -o naupd_standard $(LDFLAGS) naupd_standard.o $(LIBS) -lstdc++ -lpthread


SAUPD_DEPS = \
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/saupd.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/data.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/seupd.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/standard.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/no_op.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/copy_op.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//arpack_bindings/arpack/random_vector_op.hpp

saupd_d.o : saupd.cpp $(SAUPD_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o saupd_d.o saupd.cpp

saupd.o : saupd.cpp $(SAUPD_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o saupd.o saupd.cpp


saupd_d: saupd_d.o
	@c++ -o saupd_d $(LDFLAGS) -g saupd_d.o $(LIBS) -lstdc++ -lpthread

saupd: saupd.o
	@c++ -o saupd $(LDFLAGS) saupd.o $(LIBS) -lstdc++ -lpthread

