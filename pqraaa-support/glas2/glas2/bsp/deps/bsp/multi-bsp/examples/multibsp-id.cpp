
#include "mcbsp-templates.hpp"
#include "multibsp.hpp"

#define MULTIBSP_COMPUTER 2, 1l<<33, 2, 1l<<25, 1, 1l<<15
#define N (1l<<18)

class MultiBSP_ID : public mcbsp::MultiBSP_algorithm< MultiBSP_ID > {

	private: 

		unsigned int randomID;

		virtual void leaf_spmd() {
			//get random ID
			srand( bsp_lid() * time(NULL) );
			randomID = 10 * bsp_nleafs() * (rand()/static_cast<double>(RAND_MAX));
			//go up
			bsp_up();
		}

		virtual void nested_spmd() {
			//first report subtree of root
			bsp_up();
			//await turn (does too many syncs, but simplifies code)
			for( size_t i = 0; i < bsp_nleafs(); ++i ) {
				//if it is my turn
				if( i == bsp_lid() ) {
					//report
					std::cout << "\nNested SPMD, local ID: " << bsp_pid() << "/" << bsp_nprocs() << ". ";
					std::cout << "Global ID: " << bsp_lid() << "/" << bsp_nleafs() << ".\n";
					std::cout << "Random IDs of leaf programs in my current subtree:\n";
					for( size_t i = 0; i < bsp_sleafs(); ++i ) {
						std::cout << bsp_retrieve( bsp_lid() + i ).randomID << " ";
					}
					std::cout << "\n";
				}
				//sync
				bsp_sync();
			}
			//let child nodes report
			bsp_down();
		}

};

int main() {
	//report compilation parameters
	std::cout << "Configured using the following MultiBSP computer:\n";
	std::cout << mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::getDescription() << std::endl;

	//build program and run it
	MultiBSP_ID algorithm;
	mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::bsp_begin( algorithm );

	//done
	return 0;
}

