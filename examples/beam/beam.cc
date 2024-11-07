// 
// Project     : PQRAAA
// File        : beam.cc
// Description : test qraaa on clamped sandwich beam reference problem
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
template< typename Tval >
using Mat = Eigen::Matrix<Tval,Eigen::Dynamic,Eigen::Dynamic>;

template< typename Tval >
using Vec = Eigen::Vector<Tval,Eigen::Dynamic>;

using Tval = double;
using CTval = complex<Tval>;
int main(int argc,char** argv){
    /*
    read (vectorized, conformally-structured) sparse matrices from .dat files
    */
    char filenum[255];
    if(argc==1){
        strcpy(filenum,"168");
    }else{
        strcpy(filenum,argv[1]);
        if(!strcmp(filenum,"168")==0 & !strcmp(filenum,"840")==0 & !strcmp(filenum,"3360")==0){
            cout<<"reached1"<<endl;
            std::__throw_invalid_argument("ERROR: select 168,840 or 3360");
            return 0;
        }
    }

    bool qraaa = (strcmp(argv[2],"1")==0);

    char fmfile[255] = "include/files/beam/vm_";
    char fefile[255] = "include/files/beam/ve_";
    char fvfile[255] = "include/files/beam/vv_";

    strcat(strcat(fmfile,filenum),".dat");
    strcat(strcat(fefile,filenum),".dat");
    strcat(strcat(fvfile,filenum),".dat");
    
     
    double num; int N=0;
    ifstream fm,fe,fv;
    fm.open(fmfile);    while (fm >> num){N++;} fm.close();
    
    Vec<double> vm(N),ve(N),vv(N);
    fm.open(fmfile);    fe.open(fefile);    fv.open(fvfile);
    
    int i = 0;
    while (fm >> num){vm(i++)=num;} i=0;    fm.close();
    while (fe >> num){ve(i++)=num;} i=0;    fe.close();
    while (fv >> num){vv(i++)=num;}         fv.close(); 

    /*
    model parameters and system definition
    
    */

    double p    = .675;
    double G0   = 3.504e5;
    double Ginf = 3.062e9;
    double tau  = 8.230e-9;

    auto g = [p,G0,Ginf,tau] (Tval lam)
    {   CTval z = pow((CTval(0.,lam*tau)),p);
        return (G0 + Ginf*z)/(1. + z);
    };

    /*
    construct F

    */
    cout<<"reached"<<endl;
    int nZ = 1000;
    Eigen::ArrayXd Z = Eigen::ArrayXd::LinSpaced(nZ,200.,30000.);
    Mat<CTval> F = Mat<CTval>::Zero(nZ,N);
    for(int i=0;i<nZ;++i){
        Vec<CTval> w = (ve+g(Z(i))*vv-Z(i)*Z(i)*vm);
        F.row(i).noalias() = w.transpose();
    }
    for(int i=0;i<N;++i){
        F.col(i)/=F.col(i).norm();
    }
    /*
    set options for SV-AAA
    */
    double tol = 1e-8;
    if(qraaa){tol*=.5;}
    auto options = std::make_tuple( CORK::options::number_of_sample_points( nZ ),
                                    CORK::options::debug_level( 0 ),
                                    CORK::options::max_degree( 30 ),
                                    CORK::options::aaa_stop_criterion_max(tol)
                                    );
    glas2::contiguous_vector<Tval,int> Zglas(Z.data(),nZ);//shallow copy

    if(qraaa){
        
        /*
        compute&time QR

        */
        
        //init
        Linalg::TruncColPivHouseholderQR<Mat<CTval>> qr(F.rows(), F.cols());
        qr.setThreshold(5e-9);
        
        //compute HH factors
        auto ticQ = std::chrono::high_resolution_clock::now();
        qr.compute(F);
        
        //apply HH to Id to obtain Q
        Mat<CTval> Q = Mat<Tval>::Identity(F.rows(),qr.rank());
        auto H = qr.householderQ();
        H.setLength(qr.rank());
        Q.applyOnTheLeft(H);
        
        //extract leading subblock of R
        Mat<CTval> R0 = qr.matrixR().topLeftCorner(qr.rank(),qr.rank()).template triangularView<Eigen::Upper>();
        auto tocQ = std::chrono::high_resolution_clock::now();
        auto elapsedQ = std::chrono::duration_cast<chrono::nanoseconds>(tocQ-ticQ);
        
        // total R only used here to validate qr err., so not included in timings!
        Mat<CTval> R = qr.matrixR().topLeftCorner(qr.rank(),F.cols()).template triangularView<Eigen::Upper>();
        
        cout<<"qr done in"<<elapsedQ.count()/1e9<<endl;
        cout<<"rank : "<<qr.rank()<<endl;
        cout<<"err = "<<(F*qr.colsPermutation()-Q*R).norm()/F.norm()<<endl;

        /*
        scale and compress Q using SV-AAA
        */
        
        //scale Q
        Vec<Tval> GAMM(Q.cols());
        for(int k = 0;k<Q.cols();++k){
            Tval rkk = std::abs(R0(k,k));
            GAMM(k) = rkk;
        }
        Q.noalias()=Q*GAMM.asDiagonal();
        
        //check for NaNs
        if(Q.hasNaN()){
            std::__throw_runtime_error("Q has NaNs, clean the data");
        }

        //shallow copy to glas format (needed by SV-AAA)
        glas2::contiguous_matrix<CTval,int,glas2::column_major> Qglas(Q.data(),Q.rows(),Q.cols()); //shallow copy
        
        //apply SV-AAA and time
        auto tic = std::chrono::high_resolution_clock::now();
        auto repr_q = CORK::approximation::SV_AAA( Zglas, Qglas, CORK::approximation::sv_aaa_max(), options ) ;
        auto toc = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
        Tval timeQAAA = elapsed.count()/1e9;
        cout<<"nodes are : "<<repr_q.nodes()<<endl;
        cout<<"degr: "<<repr_q.nodes().size()<<endl;
        cout<<"computed in "<<timeQAAA<<" s"<<endl;
    }else{
        //shallow copy to glas format (needed by SV-AAA)
        glas2::contiguous_matrix<CTval,int,glas2::column_major> Fglas(F.data(),F.rows(),F.cols()); //shallow copy
        
        //apply SV-AAA and time
        auto tic = std::chrono::high_resolution_clock::now();
        auto repr_f = CORK::approximation::SV_AAA( Zglas, Fglas, CORK::approximation::sv_aaa_max(), options ) ;
        auto toc = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
        Tval timeFAAA = elapsed.count()/1e9;
        cout<<"nodes are : "<<repr_f.nodes()<<endl;
        cout<<"degr: "<<repr_f.nodes().size()<<endl;
        cout<<"computed in "<<timeFAAA<<" s"<<endl;
    }

}