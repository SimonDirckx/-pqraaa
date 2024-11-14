// 
// Project     : PQRAAA
// File        : photonic.cc
// Description : test qraaa on photonic crystal reference problem
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
    and set parameters from user-input
    */
    ////////////////////////////////////////////////////////////////////////////////////////////
    char filenum[255];
    bool qraaa      = false;
    int  n_cores    = 1;
    double tol      = 1e-8;
    if(argc==1){
        strcpy(filenum,"168");
    }else if (argc>=2){
        strcpy(filenum,argv[1]);
        if(!strcmp(filenum,"5546")==0 & !strcmp(filenum,"29034")==0 & !strcmp(filenum,"98442")==0){
            std::__throw_invalid_argument("ERROR: select 5546, 29034 or 98442");
            return 0;
        }
        if(argc>=3){
            bool qraaa = (strcmp(argv[2],"1")==0);
        }
        if(argc>=4){
            n_cores = std::max(1,std::min(8,stoi(argv[3])));
        }
        if(argc>=5){
            tol = stod(argv[4]);
        }
    }

    char matfile[255] = "include/files/photonic/M_photonic_";
    strcat(strcat(matfile,filenum),".dat");
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /*
    construct F

    */
    //////////////////////////////////////////////////////////////
    
    double eps0=1.;
    double c=2.;
    double lamP1sq  =   1.4;
    double lamP2sq  =   1.6;
    double lam01sq  =   2.5;
    double lam02sq  =   5.;
    double gamm1    =   .001;
    double gamm2    =   .02;
    auto eps = [c,lamP1sq,lamP2sq,lam01sq,lam02sq,gamm1,gamm2] (Tval x)
    {   CTval z = CTval(0.,x);
        return c+(lamP1sq/(z*z+gamm1*z+lam01sq))+(lamP2sq/(z*z+gamm2*z+lam02sq));
    };
    double num; int n=atoi(filenum); int num_coeff = 3;
    
    ifstream fm,fb;
    fm.open(matfile);
    Vec<Tval> Mvec(n*num_coeff);
    int i = 0;
    while (fm >> num){Mvec(i++)=num;} i=0;    fm.close();
    Mat<Tval> M=Eigen::Map<Mat<Tval>>(Mvec.data(),n,num_coeff);
    
    int N = n;
    int nZ = 800;
    Eigen::ArrayXd Z = Eigen::ArrayXd::LinSpaced(nZ,0,10);
    Mat<CTval> F = Mat<CTval>::Zero(nZ,N);
    
    for(int j = 0;j<nZ;++j){
        Vec<CTval> v = Vec<CTval>::Zero(num_coeff);
        v(0)=1.;
        v(1)=-Z(j)*Z(j)*eps0;
        v(2)=-Z(j)*Z(j)*eps(Z(j));
        F.row(j).noalias()=(M*v).transpose();
    }
    Vec<Tval> nrmVec(F.cols());
    for(int i=0;i<N;++i){
        nrmVec(i) = F.col(i).norm();
        F.col(i)/=nrmVec(i);
    }
    //////////////////////////////////////////////////////////////
    /*
    approximate & output info
    */
    //////////////////////////////////////////
    QRAAA::infoType info;
    QRAAA::AAAopts  opts;
    opts.tol        = tol;
    opts.max_degree = 20;
    opts.qr         = true;
    opts.n_cores    = n_cores;
    auto repr_f     = QRAAA::sv_aaa(F,Z,opts,info);
    QRAAA::summarize(info);
    //////////////////////////////////////////

    /*
    validate
    */
    /////////////////////////////////////////////////////////////////////////////////////
    int nZtest = 2513;   
    Eigen::ArrayXd Ztest = Eigen::ArrayXd::LinSpaced(nZtest,0,10);
    Vec<CTval> ftest = Vec<CTval>::Zero(N);
    Vec<CTval> rtest = Vec<CTval>::Zero(N);
    Vec<Tval> err(nZtest);
    for(int j = 0;j<nZtest;++j){
        //exact function @ Ztest(j)
        Vec<CTval> v = Vec<CTval>::Zero(num_coeff);
        v(0)=1.;    v(1)=-Ztest(j)*Ztest(j)*eps0;   v(2)=-Ztest(j)*Ztest(j)*eps(Ztest(j));
        ftest=(M*v).transpose();
        QRAAA::eval(repr_f,Ztest(j),rtest);
        for(int idx = 0;idx<rtest.size();++idx){
            rtest(idx) *= nrmVec(idx);
        }
        err(j) = (ftest-rtest).array().abs().maxCoeff()/ftest.array().abs().maxCoeff();
    }
    cout<<"Maximum error = "<<err.maxCoeff()<<endl;
    /////////////////////////////////////////////////////////////////////////////////////
    


}