#ifndef QR_AAA_HH
#define QR_AAA_HH
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
#include <set>
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
        std::chrono::nanoseconds::rep   time_total          = 0;
        int                             qr_rank             = 0;
        int                             n_supp              = 0;
        int                             n_cores             = 1;
        bool qr                                             = false;
        bool parallel                                       = false;

    };
    static infoType _default_info;
    
    //options struct
    struct AAAopts{
        double tol      = 1e-8  ;
        int max_degree  = 50    ;
        bool qr         = false ;
        int n_cores     = 8     ;//not used in qr-aaa, only for pqr-aaa
    };
    /*
    RQIZ struct
    Contains:   rational approx
                Q matrix
                indices of selected support points in Z
                support set z
    */
    template<typename argument_type,typename value_type>
    struct RQIZ{
        CORK::approximation::SV_AAA_approximation<argument_type,value_type> repr;
        Mat<value_type>     Q;
        std::vector<int>    I;
        Vec<value_type>     Z;
        
    };

    //summarize info struct
    void summarize(infoType info){
        if(info.qr){
            if(info.parallel){
                cout<<"==================="<<endl;
                cout<<"= PQR-AAA SUMMARY ="<<endl;
                cout<<"==================="<<endl;
                cout<<"# cores  : "<<info.n_cores<<endl;
                cout<<"Time QR  : "<<info.time_count_qr/1e9<<"s"<<endl;
                cout<<"Time AAA : "<<info.time_count_aaa/1e9<<"s"<<endl;
                cout<<"Time TOT.: "<<info.time_total/1e9<<"s"<<endl;
                cout<<"Max rk QR: "<<info.qr_rank<<endl;
                cout<<"#Z_m     : "<<info.n_supp<<endl;
                cout<<"==================="<<endl;
            }else{
                cout<<"=================="<<endl;
                cout<<"= QR-AAA SUMMARY ="<<endl;
                cout<<"=================="<<endl;
                cout<<"Time QR  : "<<info.time_count_qr/1e9<<"s"<<endl;
                cout<<"Time AAA : "<<info.time_count_aaa/1e9<<"s"<<endl;
                cout<<"Time TOT.: "<<info.time_total/1e9<<"s"<<endl;
                cout<<"Rk. QR   : "<<info.qr_rank<<endl;
                cout<<"#Z_m     : "<<info.n_supp<<endl;
                cout<<"=================="<<endl;
            }
            
        }else{
            cout<<"=================="<<endl;
            cout<<"= SV-AAA SUMMARY ="<<endl;
            cout<<"=================="<<endl;
            cout<<"Time AAA : "<<info.time_count_aaa/1e9<<"s"<<endl;
            cout<<"Time TOT.: "<<info.time_total/1e9<<"s"<<endl;
            cout<<"#Z_m     : "<<info.n_supp<<endl;
            cout<<"=================="<<endl;
        }
        
    }


    //qr-aaa on submatrix
    template<typename argument_type,typename value_type>
    RQIZ<argument_type,value_type> qr_aaa_RQIZ( Mat<value_type>const& F,Eigen::Array<argument_type,-1,1>& Z,int col_start,int col_stop ,AAAopts opts=AAAopts(),infoType& info=_default_info){
        auto tic_tot        = std::chrono::high_resolution_clock::now();
        int nZ  = F.rows();
        int N   = col_stop-col_start;
        info.qr = true;

        Linalg::TruncColPivHouseholderQR<Mat<value_type>> qr(nZ, N);
        qr.setThreshold(.5*opts.tol);

        auto options = std::make_tuple( CORK::options::number_of_sample_points( nZ ),
                                        CORK::options::debug_level( 0 ),
                                        CORK::options::max_degree( opts.max_degree ),
                                        CORK::options::aaa_stop_criterion_max(opts.tol) );
        
        glas2::contiguous_vector<argument_type,int> Zglas(Z.data(),nZ);//shallow copy
        auto ticQ = std::chrono::high_resolution_clock::now();
        qr.compute(F.middleCols(col_start,N));

        Mat<value_type> Q = Mat<value_type>::Identity(F.rows(),qr.rank());
        auto H = qr.householderQ();
        H.setLength(qr.rank());
        Q.applyOnTheLeft(H);
        
        Mat<value_type> R0 = qr.matrixR().topLeftCorner(qr.rank(),qr.rank()).template triangularView<Eigen::Upper>();
        auto tocQ = std::chrono::high_resolution_clock::now();
        
        info.time_count_qr  =   (std::chrono::duration_cast<chrono::nanoseconds>(tocQ-ticQ)).count();
        info.qr_rank        =   qr.rank();

        
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
        std::vector<int> I(m);
        for(int i =0;i<m;i++){
            argument_type zm = nodes(i);
            int im = 0;
            glas2::vector<value_type> fi(N);
            //find corresponding index in Z
            for(int j = 0;j<nZ;++j){
                if(std::abs(Z(j)-zm)<1e-15){
                    im = j;
                    I[i] = im;
                }
            }
            //fill
            for(int j = col_start;j<col_stop;++j){
                fi(j-col_start) = F(im,j);
            }
            repr_f.add_node(weights(i),nodes(i),fi);
        }
        RQIZ<argument_type,value_type> qraaa{repr_f,Q,I,Z};
        auto toc_tot        = std::chrono::high_resolution_clock::now();
        info.time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(toc_tot-tic_tot).count();
        return qraaa;

    }
    
    //qr-aaa on global matrix
    template<typename argument_type,typename value_type>
    RQIZ<argument_type,value_type> qr_aaa_RQIZ( Mat<value_type>const& F,Eigen::Array<argument_type,-1,1>& Z ,AAAopts opts=AAAopts(),infoType& info=_default_info){
        return  qr_aaa_RQIZ( F,Z,0,F.cols(),opts,info);   
    }
    //qr-aaa returning only final rational approx
    template<typename argument_type,typename value_type>
    CORK::approximation::SV_AAA_approximation<argument_type,value_type> qr_aaa( Mat<value_type>& F,Eigen::Array<argument_type,-1,1>& Z ,AAAopts opts=AAAopts(),infoType& info=_default_info){
        RQIZ rqiz = qr_aaa_RQIZ(F,Z,opts,info);
        return rqiz.repr;
    }
    
    //sv-aaa without qr
    template<typename argument_type,typename value_type>
        CORK::approximation::SV_AAA_approximation<argument_type,value_type> sv_aaa_no_qr( Mat<value_type>& F,Eigen::Array<argument_type,-1,1>& Z ,AAAopts opts=AAAopts(),infoType& info=infoType()){
        auto tic_tot        = std::chrono::high_resolution_clock::now();
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
        auto toc_tot        = std::chrono::high_resolution_clock::now();
        info.time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(toc_tot-tic_tot).count();
        return repr_f;
    }
    //Helper struct for Uplus
    template<typename value_type>
    struct pAAAQI{
        Mat<value_type> Q;
        std::vector<int>   I;
    };
    /*
    HELPER FUNCTIONS
    */
    std::vector<int> mockCheb(int n,int k){
        std::vector<int> pts(k);
        int m = n-1;
        for(int i = 0;i<k;++i){
            pts[i] = round(m/2 + (m/2)*cos(M_PI*(2.*i+1.)/(2.*k)));
        }
        std::sort(pts.begin(),pts.end());
        return pts;
    }


    //Uplus operation
    template<typename value_type>
    pAAAQI<value_type> Uplus(std::vector<Mat<value_type>> QVEC,std::vector<std::vector<int>>INDVEC){
        assert(QVEC.size()==INDVEC.size());
        std::set<int> dest;
        int rk = 0;
        
        std::vector<int> rkVec(INDVEC.size()+1);
        for(int i = 0;i<INDVEC.size();++i){
            rkVec[i]=rk;
            std::sort(INDVEC[i].begin(),INDVEC[i].end());
            dest.merge(std::set<int>(INDVEC[i].begin(),INDVEC[i].end()));
            rk+=QVEC[i].cols();   
        }

        std::vector<int> mock = mockCheb(QVEC[0].rows(),rk);
        std::sort(mock.begin(),mock.end());
        dest.merge( std::set<int>( mock.begin() , mock.end() ) );

        std::vector<int> I(dest.begin(),dest.end());
        rkVec[INDVEC.size()]=rk;
        Mat<value_type>Qfin = Mat<value_type>::Zero(I.size(),rk);
        for(int i = 0;i<INDVEC.size();++i){ 
            Qfin.middleCols(rkVec[i],rkVec[i+1]-rkVec[i]).noalias()=QVEC[i](I,Eigen::indexing::all);
        }
        pAAAQI<value_type> QI{Qfin,I};
        return QI;
    }


    /*
    PQR-AAA
    */

    template<typename argument_type,typename value_type>
    CORK::approximation::SV_AAA_approximation<argument_type,value_type> pqr_aaa(    Mat<value_type>& F,
                                                                                    Eigen::Array<argument_type,-1,1>& Z ,
                                                                                    AAAopts opts=AAAopts(),
                                                                                    infoType& info=_default_info){
        
        /*
        INIT
        */
        auto tic_tot        = std::chrono::high_resolution_clock::now();
        info.parallel   = true; info.qr         = true; info.n_cores    = opts.n_cores;
        
        int nZ  = F.rows();
        int N   = F.cols();

        int p = opts.n_cores;    int rem = N%p;
        int Ndivp = 0;
        
        /*
        Division of labor
        */
        
        if(rem<p/2){Ndivp = N/p;}else{Ndivp = N/p+1;}
        
        std::vector<int> indVec(p+1);
        for(int i = 0;i<p+1;++i){indVec[i]=i*Ndivp;}
        indVec[p]=N;//fix last index
        
        std::vector<int> mockIND = mockCheb(nZ,opts.max_degree);
        std::vector<Mat<value_type>> QVEC(p);   std::vector<std::vector<int>> INDVEC(p);    std::vector<QRAAA::infoType> INFOVEC(p);
        
        /*
        parallel QR-AAA
        */
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // note: not optimal->false sharing
        #pragma omp parallel num_threads(p)
        {
            #pragma omp for
            for(int INDmu = 0;INDmu<p;++INDmu){
                int ncols = indVec[INDmu+1]-indVec[INDmu];
                QRAAA::AAAopts subopts;
                QRAAA::infoType info;
                Mat<value_type> Fsub    = F.middleCols(indVec[INDmu],ncols);
                subopts.tol             = opts.tol;
                subopts.max_degree      = opts.max_degree;
                RQIZ<argument_type,value_type> rqiz = QRAAA::qr_aaa_RQIZ(Fsub,Z,subopts,info);
                QVEC[INDmu].noalias()   = rqiz.Q;
                INDVEC[INDmu]           = rqiz.I;
                INFOVEC[INDmu]          = info;
                auto repr = rqiz.repr;
            }
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /*
        Collect info
        */

        decltype(INFOVEC[0].time_count_qr)  max_qr=0;
        decltype(INFOVEC[0].time_count_aaa) max_aaa=0;
        int max_rk = 0;
        
        for(auto info:INFOVEC){
            max_qr  = std::max(max_qr,info.time_count_qr);
            max_aaa = std::max(max_aaa,info.time_count_aaa);
            max_rk=std::max(max_rk,info.qr_rank);
        }

        info.time_count_aaa=max_aaa;    info.time_count_qr=max_qr;  info.qr_rank=max_rk;
        
        /*
        Glue approximations together
        */
       
        pAAAQI QI = Uplus(QVEC,INDVEC);
        
        Mat<value_type> Qfin = QI.Q;    ArrayXd Z0 = Z(QI.I);

        QRAAA::AAAopts subopts;
        subopts.tol         =   opts.tol;
        subopts.max_degree  =   opts.max_degree;
        RQIZ<argument_type,value_type> rqiz = QRAAA::qr_aaa_RQIZ(Qfin,Z0,subopts);
        
        auto repr_q         =   rqiz.repr;
        auto nodes          =   repr_q.nodes();
        auto weights        =   repr_q.weights();
        int m               =   nodes.size();
        info.n_supp         =   m;
        
        CORK::approximation::SV_AAA_approximation<argument_type,value_type> repr_f(N,m);
        repr_f.reset(0);
        
        std::vector<int> I(m);
        for(int i =0;i<m;i++){
            argument_type zm = nodes(i);
            int im = 0;
            glas2::vector<value_type> fi(N);
            
            for(int j = 0;j<nZ;++j){
                if(std::abs(Z(j)-zm)<1e-15){
                    im = j;
                    I[i] = im;
                }
            }
            for(int j = 0;j<N;++j){
                fi(j) = F(im,j);
            }
            repr_f.add_node(weights(i),nodes(i),fi);
        }

        auto toc_tot        =   std::chrono::high_resolution_clock::now();
        info.time_total     =   std::chrono::duration_cast<std::chrono::nanoseconds>(toc_tot-tic_tot).count();
        return repr_f;
        }
    //
    template<typename argument_type,typename value_type>
    CORK::approximation::SV_AAA_approximation<argument_type,value_type> sv_aaa(    Mat<value_type>& F,
                                                                                    Eigen::Array<argument_type,-1,1>& Z ,
                                                                                    AAAopts opts=AAAopts(),
                                                                                    infoType& info=_default_info){
    if(opts.n_cores>1){
        return pqr_aaa(F,Z,opts,info);
    }else if(opts.qr)
    {
        return qr_aaa(F,Z,opts,info);
    }else
    {
        return sv_aaa_no_qr(F,Z,opts,info);
    }
    
    }
    //evaluate rational approximation
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
#endif