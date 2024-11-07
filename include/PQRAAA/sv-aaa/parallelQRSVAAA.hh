#include <vector>
#include <thread>
#include <iostream>
#include <complex>
#include <cmath>

#include <PQRAAA/Core>
#include <PQRAAA/linalg/Eigen-alt/Dense>

#include <glas2/vector.hpp>
#include <cork/matrix.hpp>
#include <cork/domain/ellipse.hpp>
#include <cork/domain/line_segment.hpp>
#include <cork/approximation/sv_aaa.hpp>
#include <cork/approximation/sv_aaa_max.hpp>
#include <cork/basis/functions.hpp>


// using namespace Eigen;
using namespace std;
using complex_t = std::complex<double>;
using mat = Eigen::Matrix<complex_t,Eigen::Dynamic,Eigen::Dynamic>;
using vec = Eigen::Vector<complex_t,Eigen::Dynamic>;
using vecd = Eigen::Vector<double,Eigen::Dynamic>;
namespace qr_sv_aaa{
class parallelQRSVAAA {
private:
    std::vector<mat> qrBlocks;
    std::vector<mat> aaaBlocks;
    int _p;
public:
    parallelQRSVAAA(mat& F, glas2::vector<double> Z, int p,double tol):_p(p){        
        //parallelQR(F,Z,p,tol);
    }
/*
    mat getQRBlock(int i) const {
        if (i < 0 || i >= qrBlocks.size()) {
            throw std::out_of_range("Invalid block index");
        }
        return qrBlocks[i];
    }

    mat getAAABlock(int i) const {
        if (i < 0 || i >= aaaBlocks.size()) {
            throw std::out_of_range("Invalid block index");
        }
        return aaaBlocks[i];
    }

pair<vector<mat>,vector<vec>> parallelQR(const mat F,glas2::vector<double> Z, int numberOfProcessors, double tol) {
    int m = F.cols() / numberOfProcessors;
    vector<mat> q_cell(numberOfProcessors);
    vector<vec> h_cell(numberOfProcessors);
    
    #pragma omp do parallel
    for (int ind = 0; ind < numberOfProcessors; ++ind) {
        auto Fi = F.block(0,m*ind,F.rows(),m);
        Eigen::ColPivHouseholderQR<mat>qr(F);
        qr.setThreshold(1e-8);
        auto compactQ = qr.householderQ().setLength(qr.nonzeroPivots());
        mat identity = mat::Identity(qr.nonzeroPivots(), qr.nonzeroPivots());
        mat Q = compactQ * identity;
        mat Rfull = qr.matrixQR().triangularView<Eigen::Upper>();
        mat R = Rfull.block(0, 0, qr.nonzeroPivots(), qr.nonzeroPivots());
        vec h = R.diagonal();
        for(int i = 0;i<h.size();++i){
            Q.col(i)*=h(i);
            cout<<h(i)<<endl;
        }
        auto options = std::make_tuple( CORK::options::number_of_sample_points( Z.size() )
                                , CORK::options::debug_level( 2 )
                                , CORK::options::max_degree( 100 )
                                ) ;

        glas2::matrix< std::complex<double> > Qc(Q.rows(),Q.cols());
        for(int i =0;i<Qc.num_rows();++i){
            for(int j =0;j<Qc.num_columns();++j){
                Qc(i,j) = Q(i,j);
            }
        }
        auto repr = CORK::approximation::SV_AAA( Z, Qc, CORK::approximation::sv_aaa_max(), options ) ;
        q_cell[ind] = Q;
        h_cell[ind] = h;
    }

    return make_pair(q_cell, h_cell);
}*/
};
}//namespace
