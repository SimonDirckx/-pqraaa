
#ifndef _H_MULTIBSP_INITIALISABLE
#define _H_MULTIBSP_INITIALISABLE

#include <mcbsp.hpp>

/**
 * A Multi-BSP datastructure that requires initialisation.
 * Base case: usually, data fields are only allocated at
 *            this leaf level.
 * @tparam _tail The MultiBSP computer parameters. Assumed
 *               empty for this leaf-case implementation.
 */
template< size_t ... _tail >
class Initialisable {

	public:

		/**
		 * Only to be called from within a BSP
		 * SPMD section.
		 */
		virtual void initialise() = 0;

		/**
		 * Only to be called from within a BSP
		 * SPMD section.
		 */
		virtual void destroy() = 0;

		/**
		 * Any serial initialision codes.
		 * Contents should be minimised.
		 * Empty by default.
		 */
		virtual void serial_initialise() {}

		/**
		 * Any serial destruction codes.
		 * Contents should be minimised.
		 * Empty by default.
		 */
		virtual void serial_destroy() {}
};

/**
 * A Multi-BSP datastructure that requires initialisation.
 * Recursive case.
 */

template< size_t _P, size_t _M, size_t ... _tail >
class Initialisable< _P, _M, _tail... > {

	public:

		/**
		 * Only to be called from within a BSP
		 * SPMD section.
		 */
		virtual void initialise() = 0;

		/**
		 * Only to be called from within a BSP
		 * SPMD section.
		 */
		virtual void destroy() = 0;

		/**
		 * Any serial initialision codes.
		 * Contents should be minimised.
		 * Empty by default.
		 */
		virtual void serial_initialise() {}

		/**
		 * Any serial destruction codes.
		 * Contents should be minimised.
		 * Empty by default.
		 */
		virtual void serial_destroy() {}

		/**
		 * Get a reference to the sub data structure. The data
		 * structure must be initialised. This code is usually
		 * only called from within an SPMD section.
		 *
		 * @param id Which child data structure to retrieve.
		 *
		 * @return A reference to the requested child data structure.
		 */
		virtual Initialisable< _tail... > & retrieve( const size_t id = bsp_pid() ) = 0;
};

#endif

