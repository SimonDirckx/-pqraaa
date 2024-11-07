
#ifndef _H_MCBSP_POSITIONAWARE_OPS
#define _H_MCBSP_POSITIONAWARE_OPS

#include "ops.hpp"
#include "identity.hpp"

//TODO doc
namespace mcbsp {
namespace ops {

	//TODO doc
	template< typename ValueType, typename IndexType = size_t >
	class positionAware {

		protected:

			//TODO doc
			ValueType _v;

			//TODO doc
			const IndexType _i;

		public:

			//TODO doc
			positionAware( const ValueType &v, const IndexType &i ) : _v( v ), _i( i ) {}

			//TODO doc
			positionAware& operator+=( const positionAware &right ) {
				_v += right._v;
				return *this;
			}

			//TODO doc
			positionAware& operator-=( const positionAware &right ) {
				_v -= right._v;
				return *this;
			}

			//TODO doc
			positionAware& operator*=( const positionAware &right ) {
				_v *= right._v;
				return *this;
			}

			//TODO doc
			positionAware& operator/=( const positionAware &right ) {
				_v /= right._v;
				return *this;
			}

			//TODO doc
			bool operator>( const positionAware &right ) {
				//first straightforward case
				if( _v > right._v ) {
					return true;
				} else if( _v == right._v ) {
					//in case of equal values, tie-break by index
					if( _i > right._i ) {
						return true;
					}
				}
				//otherwise return false
				return false;
			}

			//TODO doc
			bool operator<( const positionAware &right ) {
				//prevent code duplication, express in terms of operator>
				return (right.operator<( *this ));
			}

	};

	template< typename ValueType, typename IndexType >
	struct IDENTITY< positionAware< ValueType, IndexType > > {

		static positionAware< ValueType, IndexType > addition() {
			return positionAware< ValueType, IndexType >( IDENTITY< ValueType >::addition(), std::numeric_limits< IndexType >::max() );
		}

		static positionAware< ValueType, IndexType > subtraction() {
			return positionAware< ValueType, IndexType >( IDENTITY< ValueType >::subtraction(), std::numeric_limits< IndexType >::max() );
		}

		static positionAware< ValueType, IndexType > multiplication() {
			return positionAware< ValueType, IndexType >( IDENTITY< ValueType >::multiplication(), std::numeric_limits< IndexType >::max() );
		}

		static positionAware< ValueType, IndexType > division() {
			return positionAware< ValueType, IndexType >( IDENTITY< ValueType >::division(), std::numeric_limits< IndexType >::max() );
		}

		static positionAware< ValueType, IndexType > max() {
			return positionAware< ValueType, IndexType >( std::numeric_limits< ValueType >::min(), std::numeric_limits< IndexType >::max() );
		}

		static positionAware< ValueType, IndexType > min() {
			return positionAware< ValueType, IndexType >( std::numeric_limits< ValueType >::max(), std::numeric_limits< IndexType >::max() );
		}

	};

} //ops
} //mcbsp

#endif

