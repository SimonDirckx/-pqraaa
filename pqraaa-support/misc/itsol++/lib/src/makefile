CXXFLAGS =             -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0 -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software/ -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//ITSOL_2/INC -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++
LDFLAGS =      
LIBS =  -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0/libs -lboost_1_59_0  -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//ITSOL_2 -litsol -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/lib/ -litsol_pp

ARMS_DEPS = \
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/itsol/arms.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/itsol/spar_mat.hpp

arms.o : arms.cpp $(ARMS_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o arms.o arms.cpp

arms_d.o : arms.cpp $(ARMS_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o arms_d.o arms.cpp



ILUK_DEPS = \
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/itsol/iluk.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/itsol/spar_mat.hpp

iluk.o : iluk.cpp $(ILUK_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o iluk.o iluk.cpp

iluk_d.o : iluk.cpp $(ILUK_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o iluk_d.o iluk.cpp



ILUT_DEPS = \
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/itsol/ilut.hpp\
	/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//itsol++/itsol/spar_mat.hpp

ilut.o : ilut.cpp $(ILUT_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o ilut.o ilut.cpp

ilut_d.o : ilut.cpp $(ILUT_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o ilut_d.o ilut.cpp


