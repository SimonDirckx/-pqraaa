CXXFLAGS =       -I/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0
LDFLAGS =   
LIBS =  -L/home/simond/hlibpro-3.0-Ubuntu22.04/hlibpro-3.0/numerics_software//boost_1_59_0/libs -lboost_1_59_0

AUTO_TIMERS_CONSTRUCTION_DEPS = 

auto_timers_construction_d.o : auto_timers_construction.cpp $(AUTO_TIMERS_CONSTRUCTION_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o auto_timers_construction_d.o auto_timers_construction.cpp

auto_timers_construction.o : auto_timers_construction.cpp $(AUTO_TIMERS_CONSTRUCTION_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o auto_timers_construction.o auto_timers_construction.cpp



CPU_TIMER_DEPS = 

cpu_timer_d.o : cpu_timer.cpp $(CPU_TIMER_DEPS)
	@c++ -c -std=c++17 -Wall  -g $(CXXFLAGS) -o cpu_timer_d.o cpu_timer.cpp

cpu_timer.o : cpu_timer.cpp $(CPU_TIMER_DEPS)
	@c++ -c -std=c++17 -Wall  -DNDEBUG -O3 -funroll-loops $(CXXFLAGS) -o cpu_timer.o cpu_timer.cpp


