// 
// Project     : PQRAAA
// File        : qraaa.hh
// Description : namespace containing qr-aaa, sv-aaa (overload) and some useful helper methods
// Author      : Simon Dirckx
// Copyright   : KU Leuven Dept. CS 2024
//


#include <vector>
#include <thread>
#include <iostream>
#include <complex>
#include <cmath>
#include <glas2/vector.hpp>
#include <cork/matrix.hpp>
#include <cork/domain/ellipse.hpp>
#include <cork/domain/line_segment.hpp>
#include <cork/approximation/sv_aaa.hpp>
#include <cork/approximation/sv_aaa_max.hpp>
#include <cork/basis/functions.hpp>
#include <eigen3/Eigen/Dense>
#include <PQRAAA/linalg/truncatedQR.hh>

#include <cork/domain/line_segment.hpp>
#include <cork/eigenvalue_selector/inside_domain.hpp>
#include <cork/nonlinear_matrix.hpp>
#include <cork/eigen_pairs.hpp>
#include <cork/set_of_functions.hpp>
#include <cork/aaa_cork.hpp>
#include <cork/domain/ellipse.hpp>
#include <cork/domain/line_segment.hpp>
#include <cork/approximation/sv_aaa.hpp>
#include <cork/approximation/sv_aaa_max.hpp>
#include <cork/basis/functions.hpp>
using namespace std;
using namespace Eigen;
using namespace boost;

template< typename argument_type >
using Mat = Eigen::Matrix<argument_type,Eigen::Dynamic,Eigen::Dynamic>;

template< typename argument_type >
using Vec = Eigen::Vector<argument_type,Eigen::Dynamic>;

namespace QRAAA{
    
    //info struct
    struct infoType{
        std::chrono::nanoseconds::rep   time_count_qr       = 0;
        std::chrono::nanoseconds::rep   time_count_aaa      = 0;
        int                             qr_rank             = 0;
        int                             n_supp              = 0;
        bool qr                                             = false;
    };

    struct AAAopts{
        double tol = 1e-8;
        int max_degree = 50;
    };

    void summarize(infoType info){
        if(info.qr){
            cout<<"=================="<<endl;
            cout<<"= QR-AAA SUMMARY ="<<endl;
            cout<<"=================="<<endl;
            cout<<"Time QR  : "<<info.time_count_qr/1e9<<"s"<<endl;
            cout<<"Time AAA : "<<info.time_count_aaa/1e9<<"s"<<endl;
            cout<<"Rank QR  : "<<info.qr_rank<<endl;
            cout<<"#Z_m     : "<<info.n_supp<<endl;
            cout<<"=================="<<endl;
        }else{
            cout<<"=================="<<endl;
            cout<<"= SV-AAA SUMMARY ="<<endl;
            cout<<"=================="<<endl;
            cout<<"Time AAA : "<<info.time_count_aaa/1e9<<"s"<<endl;
            cout<<"#Z_m     : "<<info.n_supp<<endl;
            cout<<"=================="<<endl;
        }
        
    }
    
    /*
    aaa-functions
    */


    //qr-aaa
    
    template<typename argument_type,typename value_type>
    CORK::approximation::SV_AAA_approximation<argument_type,value_type> qr_aaa( Mat<value_type>& F,Eigen::Array<argument_type,-1,1>& Z ,AAAopts opts=AAAopts(),infoType& info=infoType()){
        
        int nZ  = F.rows();
        int N   = F.cols();
        info.qr = true;

        Linalg::TruncColPivHouseholderQR<Mat<value_type>> qr(F.rows(), F.cols());
        qr.setThreshold(.5*opts.tol);

        auto options = std::make_tuple( CORK::options::number_of_sample_points( nZ ),
                                        CORK::options::debug_level( 0 ),
                                        CORK::options::max_degree( opts.max_degree ),
                                        CORK::options::aaa_stop_criterion_max(opts.tol) );
        
        glas2::contiguous_vector<argument_type,int> Zglas(Z.data(),nZ);//shallow copy
        auto ticQ = std::chrono::high_resolution_clock::now();
        qr.compute(F);

        Mat<value_type> Q = Mat<value_type>::Identity(F.rows(),qr.rank());
        auto H = qr.householderQ();
        H.setLength(qr.rank());
        Q.applyOnTheLeft(H);
        
        Mat<value_type> R0 = qr.matrixR().topLeftCorner(qr.rank(),qr.rank()).template triangularView<Eigen::Upper>();
        auto tocQ = std::chrono::high_resolution_clock::now();
        
        info.time_count_qr  =   (std::chrono::duration_cast<chrono::nanoseconds>(tocQ-ticQ)).count();
        info.qr_rank        =   qr.rank();

        Mat<value_type> R   =   qr.matrixR().topLeftCorner(qr.rank(),F.cols()).template triangularView<Eigen::Upper>();

        Vec<decltype(std::abs(value_type()))> GAMM(Q.cols());
        for(int k = 0;k<Q.cols();++k){
            decltype(std::abs(value_type())) rkk = std::abs(R0(k,k));
            GAMM(k) = rkk;
        }
        Q.noalias()=Q*GAMM.asDiagonal();
        
        //check for NaNs
        if(Q.hasNaN()){
            std::__throw_runtime_error("Q has NaNs, clean the data");
        }

        //shallow copy to glas format (needed by SV-AAA)
        glas2::contiguous_matrix<value_type,int,glas2::column_major> Qglas(Q.data(),Q.rows(),Q.cols()); //shallow copy
        
        //apply SV-AAA and time
        auto tic    = std::chrono::high_resolution_clock::now();
        auto repr_q = CORK::approximation::SV_AAA( Zglas, Qglas, CORK::approximation::sv_aaa_max(), options ) ;
        auto toc    = std::chrono::high_resolution_clock::now();
        
        info.time_count_aaa =   (std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic)).count();
        auto nodes          =   repr_q.nodes();
        auto weights        =   repr_q.weights();
        int m               =   nodes.size();
        info.n_supp         =   m;
        CORK::approximation::SV_AAA_approximation<argument_type,value_type> repr_f(N,m);
        repr_f.reset(0);
        
        for(int i =0;i<m;i++){
            argument_type zm = nodes(i);
            int im = 0;
            glas2::vector<value_type> fi(N);
            //find corresponding index in Z
            for(int j = 0;j<nZ;++j){
                if(std::abs(Z(j)-zm)<1e-15){
                    im = j;
                }
            }
            //fill
            for(int j = 0;j<N;++j){
                fi(j) = F(im,j);
            }
            repr_f.add_node(weights(i),nodes(i),fi);
        }
        return repr_f;

}
        //sv-aaa overload
    
        template<typename argument_type,typename value_type>
        CORK::approximation::SV_AAA_approximation<argument_type,value_type> sv_aaa( Mat<value_type>& F,Eigen::Array<argument_type,-1,1>& Z ,AAAopts opts=AAAopts(),infoType& info=infoType()){
        
        int nZ  = F.rows();

        glas2::contiguous_vector<argument_type,int> Zglas(Z.data(),nZ);
        glas2::contiguous_matrix<value_type,int,glas2::column_major> Fglas(F.data(),F.rows(),F.cols()); //shallow copy
        auto options = std::make_tuple( CORK::options::number_of_sample_points( nZ ),
                                        CORK::options::debug_level( 0 ),
                                        CORK::options::max_degree( opts.max_degree ),
                                        CORK::options::aaa_stop_criterion_max(opts.tol) );
        //apply SV-AAA and time
        auto tic = std::chrono::high_resolution_clock::now();
        auto repr_f = CORK::approximation::SV_AAA( Zglas, Fglas, CORK::approximation::sv_aaa_max(), options ) ;
        auto toc = std::chrono::high_resolution_clock::now();
        info.time_count_aaa = (std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic)).count();
        info.n_supp = repr_f.n();
        return repr_f;

}
        template<typename argument_type,typename value_type>
        void eval(CORK::approximation::SV_AAA_approximation<argument_type,value_type> & repr,argument_type z,Vec<value_type>& v){
            int N = repr.coefficients().num_columns();
            assert(N==v.size());
            glas2::vector<value_type> vglas(N);
            repr.eval(z,vglas);
            for(int i = 0;i<N;++i){
                v(i)=vglas(i);//deep copy
            }
        }


}