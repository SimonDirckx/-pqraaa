//Zifan Liu, zifan.liu@cs.kuleuven.be

#ifndef glas2_bsp_distribution_equal_distributed_hpp
#define glas2_bsp_distribution_equal_distributed_hpp

#include <functional>
#include <type_traits> //This header defines a series of classes to obtain type information on compile-time.
#include <iostream>
//#include <glas2/bsp/deps/bsp/mcbsp.hpp>
//#include <glas2/bsp/deps/bsp/mcbsp-templates.hpp>
//#include<mcbsp-templates.hpp>
//#include <glas2/bsp/deps/bsp/mcbsp.h>
//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp.hpp"
//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp-templates.hpp"
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) \
	                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))

namespace glas2 { namespace bsp {

  class equal_distributed {	
  /*class row_distributed : public mcbsp::BSP_program {
   
   protected:
	 virtual mcbsp::BSP_program * newInstance() {
	 	return new row_distributed ();
         }
   	 virtual void spmd() {std::cout<<"Hello"<<std::endl;}
   */

   private:
	 unsigned int NumGlobalElements;
	 unsigned int NumLocalElements;
	 unsigned int MinLocalGid;
	 unsigned int MaxLocalGid;
	 unsigned int MinLid;
	 unsigned int MaxLid;
	 unsigned int size;
	 //unsigned int NumProcs;
	 //unsigned int n = 100;
	 unsigned int p;

   public:


	 equal_distributed (unsigned int n, unsigned int pid, unsigned int NumProcs) {
		//bsp_begin( bsp_nprocs() );
		//const unsigned int pid = bsp_pid();
		p = NumProcs;
		NumGlobalElements = n;
		//NumProcs = bsp_nprocs();
		NumLocalElements = BLOCK_SIZE(pid,NumProcs,NumGlobalElements);
		MinLocalGid = BLOCK_LOW(pid,NumProcs,NumGlobalElements);
		MaxLocalGid = BLOCK_HIGH(pid,NumProcs,NumGlobalElements);
		MinLid = 0;
		MaxLid = NumLocalElements - 1;
		size = BLOCK_SIZE(pid,NumProcs,NumGlobalElements);
		//bsp_sync();
		//bsp_end();
	 };

	 ~equal_distributed () {
//		 std::cout<<"Bye"<<std::endl;
	 };
 	 

	 unsigned int getNumGlobalElements() const {
		return NumGlobalElements;
	 };

	 unsigned int getNumLocalElements() {
		return NumLocalElements;
	 };

	 unsigned int getMinLocalGid() const {
		return MinLocalGid;
	 };

	 unsigned int getMaxLocalGid() const {
		return MaxLocalGid;
	 };

	 unsigned int getMinLid() {
		return MinLid;
	 };

	 unsigned int getMaxLid() {
		return MaxLid;
	 };

	 unsigned int getSize() const {
		return size;
	 };

	 unsigned int getRemoteStart (unsigned int i) const {
		return BLOCK_LOW(i,p, NumGlobalElements);
	 };

	 unsigned int getRemoteLength (unsigned int i) const {
		return BLOCK_SIZE(i, p, NumGlobalElements);
	 };

	 unsigned int getOwner (unsigned int i)  {
		return BLOCK_OWNER(i,p,NumGlobalElements);
	 };
  } ;


} } // namespace glas::bsp

#endif
