
#ifndef _H_MULTIBSP_INTERNAL
#define _H_MULTIBSP_INTERNAL

#include <iostream>

#include "multibsp.hpp"

extern "C" {
	#define MCBSP_CPLUSPLUS
	#include <mcinternal.h>
	#include <mcutil.h>
	#undef MCBSP_CPLUSPLUS
}

#endif
