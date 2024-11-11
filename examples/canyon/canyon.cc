// 
// Project     : PQRAAA
// File        : canyon.cc
// Description : test qraaa on canyon potential well Schr\"odinger reference problem
// Author      : Simon Dirckx
// Copyright   : KU Leuven Dept. CS 2024
//


#include <vector>
#include <thread>
#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include <PQRAAA/sv-aaa/qraaa.hh>

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
        strcpy(filenum,"5971");
    }else{
        strcpy(filenum,argv[1]);
        if(!strcmp(filenum,"5971")==0 & !strcmp(filenum,"7432")==0 & !strcmp(filenum,"10157")==0 & !strcmp(filenum,"15121")==0){
            cout<<"reached1"<<endl;
            std::__throw_invalid_argument("ERROR: select 5971,7432 or 15121");
            return 0;
        }
    }

    bool qraaa = (strcmp(argv[2],"1")==0);
    char matfile[255] = "include/files/canyon/M_canyon_";
    strcat(strcat(matfile,filenum),".dat");
    char bfile[255] = "include/files/canyon/bpoints_canyon_";
    strcat(strcat(bfile,filenum),".dat");
     
    double num; int n=atoi(filenum); int num_coeff = 2;
    
    ifstream fm,fb;
    fm.open(matfile);
    fb.open(bfile);
    while (fb >> num){num_coeff++;} fb.close();fb.open(bfile);
    Vec<Tval> Mvec(n*num_coeff);
    Vec<Tval> bvec(num_coeff-2);
    int i = 0;
    while (fm >> num){Mvec(i++)=num;} i=0;    fm.close();
    while (fb >> num){bvec(i++)=num;} i=0;    fb.close();
    
    Mat<Tval> M=Eigen::Map<Mat<Tval>>(Mvec.data(),n,num_coeff);
    
    int N = n;
    double m = .2;
    int nZ = 1000;
    Eigen::ArrayXd Z = Eigen::ArrayXd::LinSpaced(nZ,bvec(0)+1e-4,bvec(1)-1e-4);
    Mat<CTval> F = Mat<CTval>::Zero(nZ,N);
    
    for(int j = 0;j<nZ;++j){
        Vec<CTval> v = Vec<CTval>::Zero(num_coeff);
        v(0)=1.;
        v(1)=-Z(j);
        v(2)=exp( CTval( 0., sqrt( m*(Z(j)-bvec(0)) ) ) );
        for(int i = 3;i<num_coeff;++i){
            v(i)=exp( CTval( 0., sqrt( m*(-Z(j)+bvec(i-2)) ) ) );
        }
        F.row(j).noalias()=(M*v).transpose();
    }
    Vec<Tval> nrmVec(N);
    for(int i=0;i<N;++i){
        nrmVec(i) = F.col(i).norm(); 
        F.col(i)/=nrmVec(i);
    }
    
    double tol = 1e-8;
    
    if(qraaa){

        QRAAA::infoType info;
        QRAAA::AAAopts opts;
        opts.tol = tol;
        opts.max_degree = 30;
        auto repr_f=QRAAA::qr_aaa(F,Z,opts,info);
        QRAAA::summarize(info);

        //validate
        int nZtest = 2513;
        Eigen::ArrayXd Ztest = Eigen::ArrayXd::LinSpaced(nZtest,bvec(0)+1e-4,bvec(1)-1e-4);
        Vec<Tval> err(nZtest);
        
        Vec<CTval> ftest(N);
        Vec<CTval> rtest(N);

        for(int j = 0;j<nZtest;++j){
            Vec<CTval> v = Vec<CTval>::Zero(num_coeff);
            v(0)=1.;
            v(1)=-Ztest(j);
            v(2)=exp( CTval( 0., sqrt( m*(Ztest(j)-bvec(0)) ) ) );
            for(int i = 3;i<num_coeff;++i){
                v(i)=exp( CTval( 0., sqrt( m*(-Ztest(j)+bvec(i-2)) ) ) );
            }
            ftest=(M*v);
            //qr-aaa approximation
            QRAAA::eval(repr_f,Ztest(j),rtest);
            for(int idx = 0;idx<rtest.size();++idx){
                rtest(idx) *= nrmVec(idx);
            }
            err(j) = (ftest-rtest).array().abs().maxCoeff()/ftest.array().abs().maxCoeff();
        }
        cout<<"Maximum error = "<<err.maxCoeff()<<endl;
    }else{
        QRAAA::infoType info;
        QRAAA::AAAopts opts;
        opts.tol = tol;
        opts.max_degree = 30;
        auto repr_f=QRAAA::sv_aaa(F,Z,opts,info);
        QRAAA::summarize(info);
        
        
        //validate
        int nZtest = 2513;  Eigen::ArrayXd Ztest = Eigen::ArrayXd::LinSpaced(nZtest,bvec(0)+1e-4,bvec(1)-1e-4);
        Vec<CTval> ftest = Vec<CTval>::Zero(N);
        Vec<Tval> err(nZtest);
        Vec<CTval> rtest(N);
        
        for(int j=0;j<nZtest;++j){
            Vec<CTval> v = Vec<CTval>::Zero(num_coeff);
            v(0)=1.;
            v(1)=-Ztest(j);
            v(2)=exp( CTval( 0., sqrt( m*(Ztest(j)-bvec(0)) ) ) );
            for(int i = 3;i<num_coeff;++i){
                v(i)=exp( CTval( 0., sqrt( m*(-Ztest(j)+bvec(i-2)) ) ) );
            }
            ftest=(M*v);
            
            QRAAA::eval(repr_f,Ztest(j),rtest);
            for(int idx = 0;idx<N;++idx){rtest(idx)*=nrmVec(idx);}
            err(j) = (ftest-rtest).array().abs().maxCoeff()/ftest.array().abs().maxCoeff();
        
        }
        cout<<"Max. err.: "<<err.maxCoeff()<<endl;
    }

}