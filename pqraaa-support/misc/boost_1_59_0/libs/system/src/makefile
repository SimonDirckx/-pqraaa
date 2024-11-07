CXXFLAGS =       -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0
LDFLAGS =   
LIBS =  -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0/libs -lboost_1_59_0

ERROR_CODE_DEPS = 

error_code_d.o : error_code.cpp $(ERROR_CODE_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o error_code_d.o error_code.cpp

error_code.o : error_code.cpp $(ERROR_CODE_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o error_code.o error_code.cpp


