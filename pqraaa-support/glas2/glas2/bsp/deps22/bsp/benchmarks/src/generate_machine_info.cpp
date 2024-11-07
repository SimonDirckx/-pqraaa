
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <mcbsp.hpp>

bool choose( const std::vector< std::string > choices, unsigned short int &choice ) {
	std::cout << "\nPlease choose:\n";
	for( size_t i = 0; i < choices.size(); ++i ) {
		std::cout << i << ": " << choices[ i ] << "\n";
	}
	std::cout << "\n";
	std::cin >> choice;
	if( std::cin.fail() ) {
		std::cerr << "Error parsing user input." << std::endl;
		return false;
	}
	if( choice >= choices.size() ) {
		std::cerr << "User input out of range." << std::endl;
		return false;
	}
	std::cout << "\n";
	return true;
}

bool write( const std::vector< std::string > content ) {
	std::ofstream machine;
	machine.open( "machine.info", std::ios::out );
	if( !machine ) {
		std::cerr << "Could not open machine.info file for writing!" << std::endl;
		return false;
	}
	std::cout << "Writing the following lines to machine.info:\n";
	for( size_t i = 0; i < content.size(); ++i ) {
		std::cout << "machine.info << " << content[ i ] << "\n";
		machine << content[ i ] << std::endl;
	}

	if( machine.bad() ) {
		std::cerr << "Error during write-out to machine.info!" << std::endl;
		machine.close();
		return false;
	}

	machine.close();
	std::cout << "Done!" << std::endl;
	return true;
}

int main( int argc, char** argv ) {

	bool fail = false;
	unsigned short int choice;
	std::vector< std::string > yes_no;
	yes_no.push_back( "No" );
	yes_no.push_back( "Yes" );

	//cache of output file
	std::vector< std::string > out;

	std::cout << "Is the target architecture an Intel Xeon Phi?" << std::endl;
	if( !choose( yes_no, choice ) ) {
		std::cerr << "Error during user input, aborting." << std::endl;
		return EXIT_FAILURE;
	}
	if( choice == 0 ) {
		std::cout << "Does the target architecture have hyper-threading enabled?" << std::endl;
		if( !choose( yes_no, choice ) ) {
			std::cerr << "Error during user input, aborting." << std::endl;
			return EXIT_FAILURE;
		}
		if( choice == 1 ) {
			out.push_back( "threads_per_core 2" );
			out.push_back( "thread_numbering wrapped" );
			std::cout << "\nArchitecture has hyperthreading enabled, but employing hyperthreading does not always lead to faster computations. Please choose whether to enable hyperthreading for MulticoreBSP:" << std::endl;
			std::vector< std::string > ht_q;
			ht_q.push_back( "Let BSP use hyperthreading" );
			ht_q.push_back( "Disallow BSP from using hyperthreading." );
			if( !choose( ht_q, choice ) ) {
				std::cerr << "Error during user input, aborting." << std::endl;
				return EXIT_FAILURE;
			}
			if( choice == 1 ) {
				out.push_back( "unused_threads_per_core 1" );
			} //other option is by default
		} else {
			assert( choice == 0 ); //no hyperthreading
			std::cout << "Is this application running on the target architecture?" << std::endl;
			if( !choose( yes_no, choice ) ) {
				std::cerr << "Error during user input, aborting." << std::endl;
				return EXIT_FAILURE;
			}
			if( choice == 0 ) {
				fail = true;
			} else {
				assert( choice == 1 ); //yes
				std::cout << "The target architecture (=this architecture) supports " << bsp_nprocs() << " hardware threads." << std::endl;
				if( !choose( yes_no, choice ) ) {
					std::cerr << "Error during user input, aborting." << std::endl;
					return EXIT_FAILURE;
				}
				if( choice == 0 ) {
					fail = true;
				}
				assert( choice == 1 );
			}
		}
	} //for the other choices, defaults are ok

	if( fail ) {
		std::cout << "Not enough information to automatically infer correct affinity settings for BSP. Please refer to the changelog for documentation on the machine.info file to write your own machine.info. Refer to your cluster documentation (or contact your system administrator) if your machine's specifics are not clear you." << std::endl;
		return EXIT_SUCCESS;
	}

	//check compact or scatter affinity
	std::cout << "When using less than the full amount of available BSP threads, do you want to maximise for bandwidth or compute speed?" << std::endl;
	std::vector< std::string > compact_q;
	compact_q.push_back( "Optimise for bandwidth-bound applications (maximise available bandwidth)" );
	compact_q.push_back( "Optimise for compute-bound applications" );
	unsigned short int compact;
	if( !choose( compact_q, compact ) ) {
		std::cerr << "Error during user input, aborting." << std::endl;
		return EXIT_FAILURE;
	}
	if( compact != 0 ) {
		if( choice < 2 ) {
			std::cout << "User request optimised strategy for compute-bound applications." << std::endl;
			out.push_back( "affinity compact" );
		}
	} //for the other choices, defaults are ok

	//questions done, write output
	if( out.size() == 0 ) {
		std::cout << "Defaults suffice, no explicit machine.info required; none will be generated." << std::endl;
		return EXIT_SUCCESS;
	} else {
		if( !write( out ) ) {
			std::cerr << "Error while writing machine.info; aborting" << std::endl;
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}

