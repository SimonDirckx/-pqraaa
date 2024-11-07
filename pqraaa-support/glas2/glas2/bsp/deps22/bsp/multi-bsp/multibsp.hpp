
#include <limits>
#include <vector>
#include <sstream>
#include <assert.h>
#include <string.h>
#include <iostream>
#include <pthread.h>

extern "C" {
	#include <mcbsp.h>
}

#ifndef SIZE_MAX
 /**
  * Provide SIZE_MAX if this ANSI C99 standard macro 
  * does not make it into C++ space.
  */
 #define SIZE_MAX ((size_t)-1)
#endif

#ifndef _H_MULTIBSP
#define _H_MULTIBSP

/**
 * Namespace in which the C++-style extensions for
 * MulticoreBSP for C reside.
 */
namespace mcbsp {

	//for overview, a forward declaration of all classes defined in this file:

	//MultiBSP algorithm
	template< class T >
	class MultiBSP_algorithm;

	//MultiBSP template
	class MultiBSP_template;

	//MultiBSP computer
	template< size_t _P, size_t _M, size_t ... _tail >
	class MultiBSP_computer;

	//MultiBSP program root note
	template< class T, size_t V, size_t W, size_t ... U >
	class MultiBSP_program_root;

	//MultiBSP program leaf node (and internal nodes through template specialisation)
	template< class T, size_t ... U >
	class MultiBSP_program;

	//MultiBSP program code shared by all nodes, template-independent part
	class MultiBSP_program_algorithm_independent;

	//MultiBSP program code shared by all nodes
	template< class T >
	class MultiBSP_program_shared;

	//MultiBSP program code shared by internal nodes
	template< class T, size_t V, size_t W, size_t ... U >
	class MultiBSP_program_internal;

	//MultiBSP program code shared by non-root nodes
	template< class T >
	class MultiBSP_program_nonroot;

	//actual implementation follows:

	/** 
	 * Shared part of any MultiBSP program that does not depend on template
	 * parameters. Contains internal codes hidden from userland.
	 */
	class MultiBSP_program_algorithm_independent {

		template< class T, size_t ... U >
		friend class MultiBSP_program;

		template< class T, size_t V, size_t W, size_t ... U >
		friend class MultiBSP_program_root;

		template< class T, size_t V, size_t W, size_t ... U >
		friend class MultiBSP_program_internal;

		template< class T >
		friend class MultiBSP_program_nonroot;

		template< class T >
		friend class MultiBSP_algorithm;

		public:

			/**
			 * Declares the state of this MultiBSP program.
			 */
			enum STATE {

				/** The program is under construction. */
				UNDEFINED,

				/** The program has been called and is currently running. */
				RUNNING,

				/**
				 * The program has passed the control flow to its parents,
				 * and awaits the return of the control flow.
				 */
				SUSPENDED,

				/** 
				 * The program has passed the control flow to a nested SPMD run,
				 * and awaits the return of the control flow.
				 */
				WAITING,

				/** State set by parent program when this program should exit. */
				EXITING

			};

		private:

			/**
			 * Monitor construct that can be used in two modes:
			 *   1: IDLE->RUNNING: the local thread sleeps and waits for a parent program to wake it.
			 *   2: WAITING->RUNNING: the local thread sleeps and waits for its nested programs to finish.
			 */
			pthread_mutex_t mutex;

			/**
			 * The condition that, paired with its corresponding mutex,
			 * forms the monitor of this MultiBSP program.
			 *
			 * @see mutex.
			 */
			pthread_cond_t condition;

			/** Number of processes to spawn on the child level. */
			const size_t P;

			/** Local cache of bsp_lid() */
			unsigned int bsp_lid;

			/** Local cache of bsp_nleafs() */
			unsigned int bsp_nleafs;

			/** Local cache of bsp_sleafs() */
			unsigned int bsp_sleafs;

			/** Local cache of bsp_nprocs() */
			unsigned int bsp_P;

			/** Local cache of bsp_pid() */
			unsigned int bsp_s;

			/** MulticoreBSP for C thread data. */
			void * thread_data;

			/** PThread reference to child process 0. */
			pthread_t nested_process;

			/** PThread entry point. */
			static void * callback( void * const instance );

			/**
			 * Code dispatcher for this MultiBSP program instance.
			 */
			void monitor();

			/**
			 * User code dispatcher for this instance.
			 */
			virtual void spmd() {
				//sanity check
				assert( isRoot() );
				//raise error
				bsp_abort( "MultiBSP_program_algorithm_independent::spmd called on root node, but a root node has no user-defined algorithm (=SPMD section) associated with it!\n" );
			}

			/**
			 * Helper function that wakes up the parent program thread.
			 * Needed by the yield() and the monitor() functions.
			 *
			 * @see yield()
			 * @see monitor()
			 */
			virtual void wakeParent() const {
				//sanity check
				assert( isRoot() );
				//raise error
				bsp_abort( "MultiBSP_program_algorithm_independent::wakeParent called on root node, but a root node has no parent!\n" );
			}

		protected:

			/** Current program state. */
			volatile enum STATE state;

			/**
			 * When P>0, this program should execute nested SPMD sections.
			 * This function provides access to the children of this instance
			 * that are used to execute these nested sections.
			 * This access function is used to spawn child threads on the
			 * correct MultiBSP program instance.
			 *
			 * @param s Which child to retrieve. Note that 0 <= s < P, or the
			 *          call will fail.
			 * @return A pointer to the requested child instance.
			 */
			virtual MultiBSP_program_algorithm_independent * getChildPointer( const size_t s ) const {
				//sanity check
				assert( isLeaf() );
				//raise error
				bsp_abort( "getChildPointer( %lu ) called on leaf MultiBSP process!", (unsigned long int)s );
				//prevent compiler warning (about no returned value)
				exit( EXIT_FAILURE );
			}

			/**
			 * Halts current execution flow and start or continues that of all
			 * child programs. Calling thread will continue execution when all
			 * child programs have halted their execution flows.
			 */
			virtual void recurse() {
				//sanity check
				assert( isLeaf() );
				//raise error
				bsp_abort( "MultiBSP_program_algorithm_independent::recurse called on leaf node, but a leaf cannot recurse on children node since it has none.\n" );
			}

			/**
			 * Halts current execution flow to continue that of the parent
			 * MultiBSP program. Calling thread will halt execution until
			 * parent calls for recursion.
			 *
			 * @see recurse()
			 */
			virtual void yield() {
				//sanity check
				assert( isRoot() );
				//raise error
				bsp_abort( "MultiBSP_program_algorithm_independent::yield called on root node, but a root cannot yield to a parent node since it does not have one.\n" );
			}

			/**
			 * Helper function to determine if this MultiBSP program instance
			 * corresponds to the root program.
			 *
			 * @return True if and only if the current instance is the root
			 *         program.
			 */
			virtual bool isRoot() const = 0;

			/**
			 * Helper function to determine if this MultiBSP program instance
			 * corresponds to a leaf program.
			 *
			 * @return True if and only if the current instance is a leaf
			 * 	   program.
			 */
			virtual bool isLeaf() const = 0;

		public:

			/** Default constructor, should never be called. */
			MultiBSP_program_algorithm_independent();

			/**
			 * Base constructor.
			 *
			 * @param _P The number of MultiBSP threads needed for a nested run.
			 *           Setting _P=0 means this instance will be a leaf program.
			 * @param _subtree_p The number of MultiBSP leaf threads this subtree
			 *                   consists of.
			 * @param _L The total number of MultiBSP leaf programs.
			 */
			MultiBSP_program_algorithm_independent( const size_t _P, const unsigned int _subtree_p, const unsigned int _L );

			/** Virtual destructor. */
			virtual ~MultiBSP_program_algorithm_independent();

			/**
			 * Makes ready all MultiBSP program instances, as well as
			 * their corresponding threads.
			 */
			void init();

			/** Will let all child processes exit. */
			virtual void finalize() {
				//sanity check; doing nothing here is only correct for leaf nodes
				assert( isLeaf() );
			}

	};

	/**
	 * Interface shared by *all* MultiBSP program nodes.
	 *
	 * @tparam Algo The final algorithm this program will execute.
	 */
	template< class Algo >
	class MultiBSP_program_shared : public MultiBSP_program_algorithm_independent {

		template< class T >
		friend class MultiBSP_algorithm;

		private:

		protected:

			/**
			 * If this instance is not the root, this function can access the corresponding
			 * user-defined algorithm.
			 * This function is defined here, since user-defined code may require access
			 * to child program instances. From the MultiBSP_algorithm perspective, however,
			 * a MultiBSP_program has a single type, and this type must support getChild as
			 * well as getData; these two methods encompass program leafs, internal nodes,
			 * and the program root.
			 *
			 * @return The user-defined algorithm corresponding to this instance.
			 */
			virtual Algo& getData() {
				//sanity check
				assert( isRoot() );
				//signal failure
				bsp_abort( "MultiBSP_program_shared< Algo >::getData() called, but I am a root instance and have no associated data.\n" );
				//evade compiler warning (about not returning anything)
				exit( EXIT_FAILURE );
			}

			/**
			 * If this instance has children programs, this function can access them.
			 * This function is defined here, since user-defined code may require access
			 * to child program instances. From the MultiBSP_algorithm perspective, however,
			 * a MultiBSP_program has a single type, and this type must support getChild as
			 * well as getData; these two methods encompass program leafs, internal nodes,
			 * and the program root.
			 *
			 * @param s The child number to return.
			 * @return Reference to the requested child program instance.
			 */
			virtual MultiBSP_program_shared< Algo >& getChild( const size_t s ) const {
				//sanity check
				assert( isLeaf() );
				//signal failure
				bsp_abort( "Cannot get %luth child of a non-internal MultiBSP program!\n", (unsigned long int)s );
				//evade compiler warning (about not returning anything)
				exit( EXIT_FAILURE );
			}

			//TODO
			virtual const Algo& retrieve( const unsigned int ) const = 0;

		public:

			/** Default constructor, should never be called. */
			MultiBSP_program_shared() { assert( false ); }

			/**
			 * Default constructor, requires number of SPMD processes for nested section.
			 * Setting _P=0 results in a leaf MultiBSP program instance.
			 *
			 * @param _P The number of SPMD threads/processes for a nested section.
			 * @param _subP The number of leaf processes at this MultiBSP subtree.
			 * @param _T The total number of MultiBSP program leafs.
			 */
			MultiBSP_program_shared( const size_t _P, const unsigned int _subP, const unsigned int _T ) : MultiBSP_program_algorithm_independent( _P, _subP, _T ) {}

			/** Default destructor. */
			virtual ~MultiBSP_program_shared() {}
	};

	/**
	 * Interface shared by all non-root MultiBSP programs.
	 *
	 * @tparam Algo The final algorithm this program will execute.
	 */
	template< class Algo >
	class MultiBSP_program_nonroot : public virtual MultiBSP_program_shared< Algo > {

		template< class T, size_t ... _tail >
		friend class MultiBSP_program;

		private:

			/** Whether to delete the algorithm when this instance is destroyed. */
			const bool deleteAlgo;

			//TODO
			const Algo * * const data_array;

			//TODO
			virtual const Algo& retrieve( const unsigned int id ) const {
				return *(data_array[ id ]);
			}

			/**
			 * A non-root can implement the wakeParent function.
			 * @see MultiBSP_program_algorithm_independent::wakeParent()
			 */
			virtual void wakeParent() const {
				//if my PID is 0, wake up parent
				if( bsp_pid() == 0 ) {
					//lock monitor
					if( pthread_mutex_lock( &(MultiBSP_program_nonroot< Algo >::parent_program.mutex) ) != 0 ) {
						bsp_abort( "Error: could not lock parent program monitor!\n" );
					}

					//DBG
					//std::cout << &(MultiBSP_program_nonroot< Algo >::parent_program) << " signalled to start running again" << std::endl;

					//set parent to running
					MultiBSP_program_nonroot< Algo >::parent_program.state = MultiBSP_program_algorithm_independent::RUNNING;
					if( pthread_cond_signal( &(MultiBSP_program_nonroot< Algo >::parent_program.condition) ) != 0 ) {
						bsp_abort( "Error: could not signal parent program for continuation!\n" );
					}

					//unlock monitor
					if( pthread_mutex_unlock( &(MultiBSP_program_nonroot< Algo >::parent_program.mutex) ) != 0 ) {
						bsp_abort( "Error: could not release parent program monitor!\n" );
					}
				}
			}

		protected:

			/**
			 * Reference to the parent program;
			 * all non-root instances have a parent.
			*/
			MultiBSP_program_algorithm_independent &parent_program;

			/** The algorithm to execute. */
			Algo &algo;

			/**
			 * Retrieves the algorithm for internal use only.
			 * @return The local algorithm this MultiBSP program runs.
			 */
			virtual Algo& getData() {
				return algo;
			}

			/**
			 * A non-root can implement the yield function.
			 * @see MultiBSP_program_algorithm_independent::yield()
			 */
			virtual void yield() {
				//DBG
				//std::cout << "Program " << &(this->algo) << " entering SUSPENDED state" << std::endl;

				//only yield if parent is not root
				if( parent_program.isRoot() ) {
					return;
				}

				//sanity check
				assert( MultiBSP_program_algorithm_independent::state == MultiBSP_program_algorithm_independent::RUNNING );

				//set SUSPENDED state
				MultiBSP_program_algorithm_independent::state = MultiBSP_program_algorithm_independent::SUSPENDED;

				//wait until all siblings are SUSPENDED
				bsp_sync();

				//DBG
				//std::cout << "Program " << &(this->algo) << " yields control to parent process." << std::endl;

				//wake up parent
				wakeParent();

				//wait until parent program is done and passes control flow back to us
				while( MultiBSP_program_algorithm_independent::state == MultiBSP_program_algorithm_independent::SUSPENDED ) {
					if( pthread_cond_wait( &(MultiBSP_program_algorithm_independent::condition), &(this->mutex) ) != 0 ) {
						bsp_abort( "Error: MultiBSP_program::yield(): cannot put my thread to sleep!\n" );
					}
				}

				//sanity check: state should be running at this point
				assert( MultiBSP_program_algorithm_independent::state == MultiBSP_program_algorithm_independent::RUNNING );

				//set reference to our program
				MultiBSP_program_nonroot< Algo >::algo.attach( this );

				//continue user-defined execution flow
				return;
			}

			/** @see MultiBSP_program_algorithm_independent::isRoot() */
			virtual bool isRoot() const {
				return false;
			}

		public:

			/**
			 * Default constructor.
			 * @param _p Which MultiBSP_program is the parent of this instance.
			 * @param _a Which user-defined algorithm to execute.
			 * @param _del Whether to delete the local algorithm instance on destruction.
			 */
			MultiBSP_program_nonroot( MultiBSP_program_algorithm_independent &_p, const Algo * * const &_data, Algo &_a, const bool _del ) :
				parent_program( _p ),
				data_array( _data ),
				algo( _a ),
				deleteAlgo( _del )
			{
				//DBG
				//std::cout << "Node @ " << this << " received data array @ " << _data << "\n";
				//std::cout << "Node @ " << this << "   stored data array @ " << data_array << "\n";
			}

			/** Default destructor. */
			virtual ~MultiBSP_program_nonroot() {
				//delete algorithm if we were asked to
				if( deleteAlgo ) {
					delete &algo;
				}
			}
	};

	/** Base for MultiBSP_template */
	class MultiBSP_template {

		template< class T >
		friend class MultiBSP_algorithm;

		template< class T, size_t ... _tail >
		friend class MultiBSP_program;

		private:

			/** Locally cached total number of leaf processes in the entire run. */
			unsigned int bsp_totalP;

			/** Locally cached total number of leaf processes in this subtree. */
			unsigned int bsp_subP;

			/** Locally cached leaf ID unique amongst all leaf processes. */
			unsigned int bsp_leafID;

			/** Locally cached number of sibling SPMD processes. */
			unsigned int bsp_P;

			/** Locally cached SPMD process ID. */
			unsigned int bsp_s;

			/**
			 * Locally cached information on whether we are
			 * running on a leaf MultiBSP program.
			 */
			bool bsp_l;

		protected:

			/** Yield command flow to parent MultiBSP program. */
			virtual void bsp_up() const = 0;

			/** Yield command flow to the nested SPMD section. */
			virtual void bsp_down() const = 0;

			/**
			 * @return A unique MultiBSP thread ID amongst all
			 *         leaf processes.
			 */
			const unsigned int& bsp_lid() const;

			/**
			 * @return The total number of leaf processes running
			 *         this algorithm on this MultiBSP computer.
			 */
			const unsigned int& bsp_nleafs() const;

			/**
			 * @return The number of leaf processes running this
			 *         algorithm on this MultiBSP subtree.
			 */
			const unsigned int& bsp_sleafs() const;

			/**
			 * Faster implementation of the standard bsp_pid
			 * primitive.
			 *
			 * @return A unique MultiBSP thread ID on this
			 *         level.
			 */
			const unsigned int& bsp_pid() const;

			/**
			 * Faster implementation of the standard
			 * bsp_nprocs primitive.
			 *
			 * @return The number of sibling MultiBSP
			 *         processes on this level.
			 */
			const unsigned int& bsp_nprocs() const;

			/**
			 * Whether our current execution corresponds to a
			 * leaf MultiBSP program run.
			 * Note: disabled for now, could be useful within
			 *       template MultiBSP programs?
			 */
			const bool& bsp_leaf() const;

		public:

			/** Default destructor. */
			virtual ~MultiBSP_template();

	};

	/**
	 * A MultiBSP algoritm.
	 * User-supplied MultiBSP algorithms should be implemented as an
	 * extension of this class.
	 *
	 * @tparam FinalType The final algorithm this program will execute.
	 */
	template< class FinalType >
	class MultiBSP_algorithm : public virtual MultiBSP_template {

		template< class T, size_t ... _tail >
		friend class MultiBSP_program;

		template< class T >
		friend class MultiBSP_program_nonroot;

		private:

			/**
			 * Pointer to the corresponding MultiBSP program instance.
			 * This pointer will be set to an appropriate value by
			 * bsp_begin(). The user should never set this field
			 * explicitly.
			 *
			 * @see MultiBSP_computer::bsp_begin()
			 */
			MultiBSP_program_shared< FinalType > *my_program;

			/** User-supplied leaf SPMD code. */
			virtual void leaf_spmd() = 0;

			/** User-supplied nested SPMD code. */
			virtual void nested_spmd() = 0;

			/**
			 * Attaches a MultiBSP program to this algorithm.
			 *
			 * @param ref The MultiBSP program to attach to.
			 */
			void attach( MultiBSP_program_shared< FinalType > * const ref ) {
				//set reference
				my_program = ref;
				//cache bsp ID and #processes
				MultiBSP_template::bsp_totalP = ref->bsp_nleafs;
				MultiBSP_template::bsp_subP   = ref->bsp_sleafs;
				MultiBSP_template::bsp_leafID = ref->bsp_lid;
				MultiBSP_template::bsp_P      = ref->bsp_P;
				MultiBSP_template::bsp_s      = ref->bsp_s;
				MultiBSP_template::bsp_l      = ref->isLeaf();
				//sanity checks
				assert( bsp_nprocs() == ref->bsp_P );
				assert( bsp_pid() == ref->bsp_s );
			}

		protected:

			/** Yield command flow to the nested SPMD section. */
			virtual void bsp_down() const {
				//sanity check
				assert( my_program != NULL );
				//activate child program
				my_program->recurse();
			}

			/** Yield command flow to parent MultiBSP program. */
			virtual void bsp_up() const {
				//sanity check
				assert( my_program != NULL );
				//activate parent program
				my_program->yield();
			}

			/** @see MultiBSP_template::bsp_retrieve */
			virtual const FinalType& bsp_retrieve( const unsigned int id ) const {
				//sanity check
				assert( my_program != NULL );
				//range check
				if( id < MultiBSP_template::bsp_lid() || id > MultiBSP_template::bsp_lid() + MultiBSP_template::bsp_sleafs() ) {
					bsp_abort( "Requested data of leaf process ID %ud, which is not contained in the subtree of this process (%ud/%ud/%ud).\n",
						id,
						MultiBSP_template::bsp_lid(),
						MultiBSP_template::bsp_sleafs(),
						MultiBSP_template::bsp_nleafs()
					);
				}
				//ask associated program for the requested leaf ID program
				return my_program->retrieve( id );
			}

		public:

			/**
			 * Default constructor. Extending classes must keep
			 * this default constructor.
			 * Use bsp_begin to execute an algorithm instance.
			 */
			MultiBSP_algorithm() : my_program( NULL ) {}

			/** Default destructor. */
			virtual ~MultiBSP_algorithm() {}
	};

	/**
	 * Part of the MultiBSP_program implementation shared by all internal nodes.
	 */
	template< class Algo, size_t _P, size_t _M, size_t ... _tail >
	class MultiBSP_program_internal : public virtual MultiBSP_program_shared< Algo > {

		template< class T, size_t V, size_t W, size_t ... U >
		friend class MultiBSP_program_root;

		private:

		protected:

			/**
			 * Pointers to the child MultiBSP program instances that make up
			 * a nested SPMD run.
			 * Note this cannot be a raw array since MultiBSP_program nodes
			 * do not have their default constructors enabled.
			 */
			std::vector< MultiBSP_program< Algo, _tail... > * > children;

			/* @see MultiBSP_program_algorithm_shared< Algo >::getChild( const size_t s ) */
			virtual MultiBSP_program_shared< Algo >& getChild( const size_t s ) const {
				assert( s < _P );
				return *children[ s ];
			}

			/* @see MultiBSP_program_algorithm_independent::getChildPointer( const size_t s ) */
			virtual MultiBSP_program_algorithm_independent * getChildPointer( const size_t s ) const {
				assert( s < _P );
				return children[ s ];
			}

			/* @see MultiBSP_program_algorithm_independent::isLeaf() */
			virtual bool isLeaf() const {
				return false;
			}

		public:

			/**
			 * Default constructor.
			 * Note that any class implicitly inheriting from MultiBSP_program_shared,
			 * including any class extending this class, has to call the constructor
			 * of MultiBSP_program_shared first (and explicitly).
			 */
			MultiBSP_program_internal() {}

			/** Default destructor. */
			virtual ~MultiBSP_program_internal() {
				//clean up C++ objects
				for( size_t i = 0; i < children.size(); ++i ) {
					delete children[ i ];
				}
			}

			/** @see MultiBSP_program_algorithm_independent::finalize() */
			virtual void finalize() {
				//sanity check
				if( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ 0 ]->state == MultiBSP_program_algorithm_independent::SUSPENDED ) {
					bsp_abort( "Error: upper-level SPMD run exited while a nested run was waiting for a bsp_down!\n" );
				}
				assert( (MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ 0 ]->state) == MultiBSP_program_algorithm_independent::EXITING );
				//go and wait until child program is done
				if( pthread_join( MultiBSP_program_algorithm_independent::nested_process, NULL ) != 0 ) {
					bsp_abort( "Error while waiting for exit of nested SPMD run!\n" );
				}
				//done
			}

	};

	/**
	 * Actual MultiBSP_program.
	 * Automatically generated; users should not extend this class.
	 * This definition corresponds to the algorithmic leaf case.
	 */
	template< class Algo, size_t ... _tail >
	class MultiBSP_program : public MultiBSP_program_nonroot< Algo > {

		protected:

			/* @see MultiBSP_program_algorithm_independent::spmd() */
			virtual void spmd() {
				//DBG
				//std::cout << "Program " << &(this->algo) << " calls leaf SPMD." << std::endl;
				//derive our leaf ID
				MultiBSP_program_algorithm_independent::bsp_lid =
					MultiBSP_program_nonroot< Algo >::parent_program.MultiBSP_program_algorithm_independent::bsp_lid * bsp_nprocs() + bsp_pid();
				//DBG
				//std::cout << "Leaf node @ " << &(MultiBSP_program_algorithm_independent::bsp_lid) << " writes to data_array @ " << MultiBSP_program_nonroot< Algo >::data_array << " at offset " << MultiBSP_program_algorithm_independent::bsp_lid << "\n";
				//DBG
				//record this leaf data globally
				MultiBSP_program_nonroot< Algo >::data_array[ MultiBSP_program_algorithm_independent::bsp_lid ] = &(MultiBSP_program_nonroot< Algo >::algo);
				//let the algorithm drive our program
				MultiBSP_program_nonroot< Algo >::algo.attach( this );
				//sync to make sure data_arrays are initialised before proceeding
				MultiBSP_program_nonroot< Algo >::yield(); //direct call, no need for dynamic dispatching through virtual functions here
				//call user-defined code
				static_cast< MultiBSP_algorithm< Algo > * >( &(MultiBSP_program_nonroot< Algo >::algo) )->leaf_spmd();
			}	

			/* @see MultiBSP_program_algorithm_independent::isLeaf() */
			virtual bool isLeaf() const {
				return true;
			}

		public:

			/**
			 * MultiBSP program leaf node constructor. Constructs default local algorithm instance.
			 *
			 * @param _T The total number of MultiBSP leaf programs.
			 * @param _p Reference to the parent MultiBSP program node.
			 * @param _array Array of algorithms used to enable bsp_retrieve.
			 */
			MultiBSP_program( const unsigned int _T, MultiBSP_program_algorithm_independent &_p, const Algo * * const &_array ) :
				MultiBSP_program( _T, _p, _array, *(new Algo()), true ) {}

			/**
			 * MultiBSP program leaf node constructor. Uses a given algorithm as the local algorithm.
			 *
			 * @param _T The total number of MultiBSP leaf programs.
			 * @param _p Reference to the parent MultiBSP program node.
			 * @param _array Array of algorithms used to enable bsp_retrieve.
			 * @param _algo Local algorithm instance.
			 */
			MultiBSP_program( const unsigned int _T, MultiBSP_program_algorithm_independent &_p, const Algo * * const &_array, Algo &_algo ) :
				MultiBSP_program( _T, _p, _array, _algo, false ) {}

			/**
			 * MultiBSP program leaf node constructor. Uses a given algorithm as the local algorithm.
			 *
			 * @param _T The total number of MultiBSP leaf programs.
			 * @param _p Reference to the parent MultiBSP program node.
			 * @param _array Array of algorithms used to enable bsp_retrieve.
			 * @param _algo Local algorithm instance.
			 * @param _del Whether to delete the local algorithm instance on destruction.
			 */
			MultiBSP_program( const unsigned int _T, MultiBSP_program_algorithm_independent &_p, const Algo * * const &_array, Algo &_algo, const bool _del ) :
				MultiBSP_program_shared< Algo >( 0, 1, _T ),
				MultiBSP_program_nonroot< Algo >( _p, _array, _algo, _del ) {
				//DBG
				//std::cout << "Leaf node constructor was passing around data array @ " << _array << "\n";
			}

	};

	/**
	 * Actual MultiBSP_program.
	 * Automatically generated; users should not extend this class.
	 * This definition corresponds to the nested case.
	 */
	template< class Algo, size_t _P, size_t _M, size_t ... _tail >
	class MultiBSP_program< Algo, _P, _M, _tail... > : public MultiBSP_program_internal< Algo, _P, _M, _tail... >, public MultiBSP_program_nonroot< Algo > {

		private:

		protected:

			/** @see MultiBSP_program_algorithm_independent::recurse() */
			virtual void recurse() {

				//for symmetry with yield(), we sync also
				bsp_sync();

				//DBG
				//std::cout << "Recurse called!" << std::endl;

				//sanity check
				assert( MultiBSP_program_algorithm_independent::state == MultiBSP_program_algorithm_independent::RUNNING );

				//DBG
				//std::cout << "Process " << &(this->algo) << " enters WAITING state" << std::endl;

				//set my state to WAITING (for nested SPMD run)
				MultiBSP_program_algorithm_independent::state = MultiBSP_program_algorithm_independent::WAITING;

				//set all children to RUNNING state
				for( size_t k = 0; k < _P; ++k ) {

					//spinlock until child is initialised
					while( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state == MultiBSP_program_algorithm_independent::UNDEFINED ) {}

					//sanity checks
					assert( (MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state) != MultiBSP_program_algorithm_independent::WAITING );
					if( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state == MultiBSP_program_algorithm_independent::EXITING ) {
						bsp_abort( "Error: bsp_down() issued while nested SPMD program has already exited.\n" );
					}

					//get child lock
					const int rc = pthread_mutex_lock( &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->mutex) );
					if( rc != 0 ) {
						bsp_abort( "Error: MultiBSP_program::recurse(), nested case: cannot lock child monitor (%p, %s)!\n", &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->mutex), strerror( rc ) );
					}

					//sanity checks
					if( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state == MultiBSP_program_algorithm_independent::WAITING ) {
						bsp_abort( "Child MultiBSP program in illegal state while trying to wake it on recursion (WAITING)!\n" );
					}
					if( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state == MultiBSP_program_algorithm_independent::EXITING ) {
						bsp_abort( "Child MultiBSP program in illegal state while trying to wake it on recursion (EXITING)!\nAn internal SPMD section called bsp_down while the child program has already exited." );
					}

					//set child state to running
					MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state = MultiBSP_program_algorithm_independent::RUNNING;

					//wake up child thread
					if( pthread_cond_signal( &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->condition) ) != 0 ) {
						bsp_abort( "Error: can not signal child MultiBSP thread to continue execution!\n" );
					}

					//yield child lock
					if( pthread_mutex_unlock( &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->mutex) ) != 0 ) {
						bsp_abort( "Error: MultiBSP_program::recurse(), nested case: cannot unlock child monitor!\n" );
					}

				}

				//wait until signaled by child with PID 0
				while( MultiBSP_program_algorithm_independent::state != MultiBSP_program_algorithm_independent::RUNNING ) {
					//DBG
					//std::cout << "Program " << &(this->algo) << " goes to sleep" << std::endl;
					if( pthread_cond_wait( &(MultiBSP_program_algorithm_independent::condition), &(MultiBSP_program_algorithm_independent::mutex) ) != 0 ) {
						bsp_abort( "Error: waiting for child MultiBSP process failed!\n" );
					}
					//DBG
					//std::cout << "Program " << &(this->algo) << " woke up" << std::endl;
				}

				//nested run at child 0 is done, so all children should be in SUSPENDED or EXITING mode, check
				for( size_t k = 0; k < _P; ++k ) {
					assert( (MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state) != MultiBSP_program_algorithm_independent::UNDEFINED );
					assert( (MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state) != MultiBSP_program_algorithm_independent::RUNNING );
					assert( (MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state) != MultiBSP_program_algorithm_independent::WAITING );
				}
				
				//set algorithm program reference to us
				MultiBSP_program_nonroot< Algo >::algo.attach( this );

				//return to user-defined execution flow
				return;
			}

			/** @see MultiBSP_program_algorithm_independent::spmd() */
			virtual void spmd() {
				//DBG
				//std::cout << "Program " << &(this->algo) << " calls nested SPMD." << std::endl;
				//derive leaf ID (assume bsp_lid of the root program node is 0!)
				MultiBSP_program_algorithm_independent::bsp_lid =
					MultiBSP_program_nonroot< Algo >::parent_program.MultiBSP_program_algorithm_independent::bsp_lid * bsp_nprocs() + bsp_pid();
				//DBG
				//std::cout << "Internal node @ " << &(MultiBSP_program_algorithm_independent::bsp_lid);
				//std::cout << " has data_array @ " << MultiBSP_program_nonroot< Algo >::data_array << "\n";
				
				//MultiBSP programs start at the leaf level
				recurse();
				//let our algorithm drive our program
				MultiBSP_program_nonroot< Algo >::algo.attach( this );
				//sync to make sure data_arrays are initialised before proceeding
				MultiBSP_program_nonroot< Algo >::yield();
				MultiBSP_program< Algo, _P, _M, _tail... >::recurse();
				//have to cast back to base class, since nested_spmd should be private in the final algorithm type
				static_cast< MultiBSP_algorithm< Algo >* >( &(MultiBSP_program_nonroot< Algo >::algo) )->nested_spmd();
			}

		public:

			/**
			 * Constructs a MultiBSP program non-root internal node.
			 * Constructs default local algorithm instance.
			 *
			 * @param _T The total number of MultiBSP leaf programs.
			 * @param _p Reference to the parent MultiBSP program node.
			 * @param _array Array of algorithms used to enable bsp_retrieve.
			 */
			MultiBSP_program( const unsigned int _T, MultiBSP_program_algorithm_independent &_p, const Algo * * const &_array ) :
				MultiBSP_program( _T, _p, _array, *(new Algo()), true ) {}

			/**
			 * MultiBSP program non-root internal node constructor.
			 * Uses a given algorithm as the local algorithm.
			 *
			 * @param _T The total number of MultiBSP leaf programs.
			 * @param _p Reference to the parent MultiBSP program node.
			 * @param _array Array of algorithms used to enable bsp_retrieve.
			 * @param _algo Local algorithm instance.
			 */
			MultiBSP_program( const unsigned int _T, MultiBSP_program_algorithm_independent &_p, const Algo * * const &_array, Algo &_algo ) :
				MultiBSP_program( _T, _p, _array, _algo, false ) {}

			/**
			 * MultiBSP program non-root internal node constructor.
			 * Uses a given algorithm as the local algorithm.
			 *
			 * @param _T The total number of MultiBSP leaf programs.
			 * @param _p Reference to the parent MultiBSP program node.
			 * @param _array Array of algorithms used to enable bsp_retrieve.
			 * @param _algo  Local algorithm instance.
			 * @param _del Whether to delete the local algorithm instance on destruction.
			 */
			MultiBSP_program( const unsigned int _T, MultiBSP_program_algorithm_independent &_p, const Algo * * const &_array, Algo &_algo, const bool _del ) :
				MultiBSP_program_shared< Algo >( _P, MultiBSP_computer< _P, _M, _tail... >::bsp_nleafs(), _T ),
				MultiBSP_program_internal< Algo, _P, _M, _tail... >(),
				MultiBSP_program_nonroot< Algo >( _p, _array, _algo, _del )
			{
				//DBG
				//std::cout << "Internal node constructor is passing around data array @ " << MultiBSP_program_nonroot< Algo >::data_array << "\n";
				MultiBSP_program_internal< Algo, _P, _M, _tail... >::children.push_back(
					new MultiBSP_program< Algo, _tail... >( _T, *this, MultiBSP_program_nonroot< Algo >::data_array, _algo )
				);
				for( size_t i = 1; i < _P; ++i ) {
					MultiBSP_program_internal< Algo, _P, _M, _tail... >::children.push_back(
						new MultiBSP_program< Algo, _tail... >( _T, *this, MultiBSP_program_nonroot< Algo >::data_array )
					);
				}
			}

			/** Default destructor. */
			virtual ~MultiBSP_program() {}
	};

	/**
	 * Actual MultiBSP_program.
	 * Automatically generated; users should not extend this class.
	 * This definition corresponds to the root case.
	 */
	template< class Algo, size_t _P, size_t _M, size_t ... _tail >
	class MultiBSP_program_root : public MultiBSP_program_internal< Algo, _P, _M, _tail... > {

		protected:

			//TODO
			const Algo * data_array[ MultiBSP_computer< _P, _M, _tail... >::bsp_nleafs() ];

			/** @see MultiBSP_program_algorithm_independent::isRoot() */
			virtual bool isRoot() const {
				return true;
			}

			//TODO
			virtual const Algo& retrieve( const unsigned int id ) const {
				return *(data_array[ id ]);
			}

		public:

			/**
			 * Constructs a MultiBSP program root node.
			 * Constructs default local algorithm instance.
			 */
			MultiBSP_program_root() : MultiBSP_program_root( *(new Algo()) ) {}

			/**
			 * Constructs a MultiBSP program root node.
			 * Uses a given algorithm as the algorithm for the SPMD process with PID 0.
			 *
			 * @param _T The total number of MultiBSP leaf programs.
			 * @param _algo The local algorithm instance for SPMD process 0.
			 */
			MultiBSP_program_root( const unsigned int _T, Algo &_algo ) :
				MultiBSP_program_shared< Algo >( _P, MultiBSP_computer< _P, _M, _tail... >::bsp_nleafs(), _T ),
				MultiBSP_program_internal< Algo, _P, _M, _tail... >()
			{
				//DBG
				//std::cout << "MultiBSP_program_root constructor called, P was set to " << MultiBSP_program_algorithm_independent::P << std::endl;
				//std::cout << "MultiBSP root has data_array @ " << data_array << "\n";

				MultiBSP_program_internal< Algo, _P, _M, _tail... >::children.push_back( new MultiBSP_program< Algo, _tail... >( _T, *this, data_array, _algo ) );
				for( size_t i = 1; i < _P; ++i ) {
					MultiBSP_program_internal< Algo, _P, _M, _tail... >::children.push_back( new MultiBSP_program< Algo, _tail... >( _T, *this, data_array ) );
				}
				const int rc1 = pthread_mutex_init( &(MultiBSP_program_algorithm_independent::mutex), NULL );
				const int rc2 = pthread_cond_init( &(MultiBSP_program_algorithm_independent::condition), NULL );
				if( rc1 != 0 || rc2 != 0 ) {
					bsp_abort( "Error: cannot create MultiBSP root monitor!\n" );
				}
			}

			/** Default destructor. */
			virtual ~MultiBSP_program_root() {
				//destroy local monitor
				const int rc1 = pthread_mutex_destroy( &(MultiBSP_program_algorithm_independent::mutex) );
				const int rc2 = pthread_cond_destroy( &(MultiBSP_program_algorithm_independent::condition) );
				if( rc1 != 0 || rc2 != 0 ) {
					std::cerr << "Non-fatal error: could not destroy MultiBSP root monitor!\n";
				}
				//done
			}

			/** Starts the MultiBSP algorithm. */
			void start() {
				//DBG
				//std::cout << "start() called!" << std::endl;

				if( pthread_mutex_lock( &(MultiBSP_program_algorithm_independent::mutex) ) != 0 ) {
					bsp_abort( "Cannot gain lock on root MultiBSP program!\n" );
				}

				//set all children to RUNNING state
				for( size_t k = 0; k < _P; ++k ) {

					//first spinlock until child is initialised
					while( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[k]->state == MultiBSP_program_algorithm_independent::UNDEFINED ) {}

					//sanity check
					assert( (MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[k]->state) == MultiBSP_program_algorithm_independent::SUSPENDED );

					//get lock
					if( int rc = pthread_mutex_lock( &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->mutex) ) != 0 ) {
						bsp_abort( "Error: MultiBSP_program_root::bsp_begin(): cannot lock child monitor (%p, %s)!\n", &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->mutex), strerror( rc ) );
					}

					//set child to RUNNING and wake up its thread
					MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state = MultiBSP_program_algorithm_independent::RUNNING;
					if( pthread_cond_signal( &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->condition) ) != 0 ) {
						bsp_abort( "Error: can not signal child MultiBSP thread to continue execution!\n" );
					}

					//yield child lock
					if( pthread_mutex_unlock( &(MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->mutex) ) != 0 ) {
						bsp_abort( "Error: MultiBSP_program_root::bsp_begin(): cannot unlock child monitor!\n" );
					}
				}

				//set my state to WAITING (for nested SPMD run)
				MultiBSP_program_algorithm_independent::state = MultiBSP_program_algorithm_independent::WAITING;

				//enter sleep, will be woken up by child program 0 when it is EXITING
				while( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ 0 ]->state != MultiBSP_program_algorithm_independent::EXITING ) {
					//prevent yield past root
					if( MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ 0 ]->state == MultiBSP_program_algorithm_independent::SUSPENDED ) {
						bsp_abort( "Error: child program passed control flow beyond the upper-level SPMD run!\n" );
					}
					//do sleep
					if( pthread_cond_wait( &(MultiBSP_program_algorithm_independent::condition), &(MultiBSP_program_algorithm_independent::mutex) ) != 0 ) {
						bsp_abort( "Error: waiting for child MultiBSP process failed!\n" );
					}
					//DBG
					//std::cout << "Root woke up" << std::endl;
				}

				//child 0 is done, since this is the root node, all children should thus be in EXITING mode
				for( size_t k = 0; k < _P; ++k ) {
					const enum MultiBSP_program_algorithm_independent::STATE k_state = 
						MultiBSP_program_internal< Algo, _P, _M, _tail... >::children[ k ]->state;
					assert( k_state == MultiBSP_program_algorithm_independent::EXITING );
				}

				//set my state to EXITING
				MultiBSP_program_algorithm_independent::state = MultiBSP_program_algorithm_independent::EXITING;

				//call finalisation
				MultiBSP_program_internal< Algo, _P, _M, _tail... >::finalize();

				//release lock
				if( pthread_mutex_unlock( &(MultiBSP_program_algorithm_independent::mutex) ) != 0 ) {
					bsp_abort( "Cannot release MultiBSP_program_root monitor!\n" );
				}

				//done
				return;
			}

	};

	/**
	 * @param _M The number of bytes we want to print prettily.
	 *
	 * @return A pretty version of _M
	 */
	static std::string outputBytes( const size_t _M ) {
		std::ostringstream oss;
		if( _M < (1l<<10 ) ) {
			oss << _M << " B";
		} else if( _M < (1l<<20) ) {
			oss << (_M>>10) << " kB";
		} else if( _M < (1l<<30) ) {
			oss << (_M>>20) << " MB";
#if defined __x86_64__ || defined __ppc64__ || _WIN64
		} else if( _M < (1l<<40) ) {
			oss << (_M>>30) << " GB";
		} else {
			oss << (_M>>40) << " TB";
		}
#else
		} else {
			oss << (_M>>30) << " GB";
		}
#endif
		return oss.str();
	}

	/**
	 * A MultiBSP computer description; recursive definition.
	 * @tparam _P    The top-level number of processors.
	 * @tparam _M    The top-level amount of available memory.
	 * @tparam _tail The remainder MultiBSP parameters.
	 */
	template< size_t _P, size_t _M, size_t ... _tail >
	class MultiBSP_computer {

		protected:

		public:

			/** @return A pretty description of this MultiBSP computer. */
			static std::string getDescription() {
				std::ostringstream oss;
				oss << "( " << _P << ", " << outputBytes( _M ) << ") -> " << MultiBSP_computer< _tail... >::getDescription();
				return oss.str();
			}

			/**
			 * Used to create an executable MultiBSP program out of the
			 * supplied MultiBSP algorithm. The algorithm is automatically
			 * expanded to fit this MultiBSP computer.
			 *
			 * @param algo This will be the execution space for the top-level
			 *             BSP thread with PID 0.
			 */
			template< class MultiBSP_Algo >
			static void bsp_begin( MultiBSP_Algo &algo ) {
				MultiBSP_program_root< MultiBSP_Algo, _P, _M, _tail... > program = MultiBSP_program_root< MultiBSP_Algo, _P, _M, _tail... >( static_cast< unsigned int >( bsp_nleafs() ), algo );
				program.init();
				program.start();
			}

			/** @return The number of leaf nodes of this MultiBSP computer. */
			static constexpr size_t bsp_nleafs() {
				return _P * MultiBSP_computer< _tail... >::bsp_nleafs();
			}

	};

	/** Non-recursive case. */
	template< size_t _P, size_t _M >
	class MultiBSP_computer< _P, _M > {

		protected:

		public:

			/** @return A pretty description of this MultiBSP computer. */
			static std::string getDescription() {
				std::ostringstream oss;
				oss << "( " << _P << ", " << outputBytes( _M ) << ")";
				return oss.str();
			}

			/**
			 * Used to create an executable MultiBSP program out of the
			 * supplied MultiBSP algorithm. The algorithm is automatically
			 * expanded to fit this MultiBSP computer.
			 *
			 * @return A MultiBSP_program that is ready for execution.
			 */
			template< class MultiBSP_Algo >
			static MultiBSP_program< MultiBSP_Algo, _P, _M > compile( MultiBSP_Algo &algo ) {
				return MultiBSP_program< MultiBSP_Algo, _P, _M >( algo );
			}

			/**
			 * Used to create an executable MultiBSP program out of the
			 * supplied MultiBSP algorithm. The algorithm is automatically
			 * expanded to fit this MultiBSP computer.
			 *
			 * @param algo This will be the execution space for the top-level
			 *             BSP thread with PID 0.
			 */
			template< class MultiBSP_Algo >
			static void bsp_begin( MultiBSP_Algo &algo ) {
				MultiBSP_program_root< MultiBSP_Algo, _P, _M > program = MultiBSP_program_root< MultiBSP_Algo, _P, _M >( static_cast< unsigned int >( bsp_nleafs() ), algo );
				program.init();
				program.start();
			}

			/** @return The number of leaf nodes of this MultiBSP computer. */
			static constexpr size_t bsp_nleafs() {
				return _P;
			}
	};

	//------ inline function declarations: ------

	inline const unsigned int& MultiBSP_template::bsp_lid() const {
		return MultiBSP_template::bsp_leafID;
	}

	inline const unsigned int& MultiBSP_template::bsp_nleafs() const {
		return MultiBSP_template::bsp_totalP;
	}

	inline const unsigned int& MultiBSP_template::bsp_sleafs() const {
		return MultiBSP_template::bsp_subP;
	}

	inline const unsigned int& MultiBSP_template::bsp_pid() const {
		return MultiBSP_template::bsp_s;
	}

	inline const unsigned int& MultiBSP_template::bsp_nprocs() const {
		return MultiBSP_template::bsp_P;
	}

	inline const bool& MultiBSP_template::bsp_leaf() const {
		return MultiBSP_template::bsp_l;
	}

} //mcbsp

#endif

