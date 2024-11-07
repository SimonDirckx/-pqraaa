
#include "mcbsp.hpp"

#include <iostream>

class Hello_World: public mcbsp::BSP_program {

	protected:

		virtual void spmd( void );

		virtual BSP_program * newInstance();

	public:

		Hello_World();
};

void Hello_World::spmd() {
	std::cout << "Hello world from thread " << bsp_pid() << std::endl;
}

mcbsp::BSP_program * Hello_World::newInstance() {
	return new Hello_World();
}

Hello_World::Hello_World() {}

int main() {
	Hello_World p;
	p.begin( 2 );
	p.begin();
	return 0;
}

