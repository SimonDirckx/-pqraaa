CXXFLAGS =       -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0
LDFLAGS =   
LIBS =  -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0/libs -lboost_1_59_0

CHRONO_DEPS = 

chrono_d.o : chrono.cpp $(CHRONO_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o chrono_d.o chrono.cpp

chrono.o : chrono.cpp $(CHRONO_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o chrono.o chrono.cpp



PROCESS_CPU_CLOCKS_DEPS = 

process_cpu_clocks_d.o : process_cpu_clocks.cpp $(PROCESS_CPU_CLOCKS_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o process_cpu_clocks_d.o process_cpu_clocks.cpp

process_cpu_clocks.o : process_cpu_clocks.cpp $(PROCESS_CPU_CLOCKS_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o process_cpu_clocks.o process_cpu_clocks.cpp



THREAD_CLOCK_DEPS = 

thread_clock_d.o : thread_clock.cpp $(THREAD_CLOCK_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o thread_clock_d.o thread_clock.cpp

thread_clock.o : thread_clock.cpp $(THREAD_CLOCK_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o thread_clock.o thread_clock.cpp


