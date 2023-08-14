//
// Created by Keynesh Dongol on 6/8/23.
//

#ifndef POLYNOMIAL_COMPANIONMATRIX_H
#define POLYNOMIAL_COMPANIONMATRIX_H


#include <iostream>
#include <Eigen/Dense>

/**
 * The code is fundamentally straightforward. We need ot write a method that creates the companion matrix where needed.
 * -> In order for us to do that, we need to figure out the dimensions and then pass a polynomials along the diagonals.
 * Therfore,
 * */



class CompanionMatrix{
public:
    //Contruction and solving companion matrix
    Eigen::Matrix<double, 1, Eigen::Dynamic> findRealRoots() {

        if (coeff(4) != 0) {
            Eigen::Matrix<double,4,4> newCompanionMatrix = Eigen::Map<const Eigen::Matrix4d>(zeroVector.data());
            newCompanionMatrix.block<3,3>(1,0).setIdentity();
            newCompanionMatrix.col(3) = - coeff.head(4) / coeff(4);

            return EigenRealRoots(newCompanionMatrix);

        } else if ( coeff(3) != 0 ) {
            Eigen::Matrix<double,3,3> newCompanionMatrix = Eigen::Map<const Eigen::Matrix3d>(zeroVector.data());
            newCompanionMatrix.block<2,2>(1,0).setIdentity();
            newCompanionMatrix.col(2) = - coeff.head(3) / coeff(3);

            return EigenRealRoots(newCompanionMatrix);

        }else if ( coeff(2) != 0) {
            Eigen::Matrix<double,2,2> newCompanionMatrix = Eigen::Map<const Eigen::Matrix2d>(zeroVector.data());
            newCompanionMatrix(1,0) = 1;
            newCompanionMatrix.col(1) = - coeff.head(2) / coeff(2);

            return EigenRealRoots(newCompanionMatrix);

        } else if ( coeff(1) != 0) {
            Eigen::Matrix<double,1,Eigen::Dynamic> realRoot; //Discuss with prof to improve this
            realRoot.resize(Eigen::NoChange, 1);
            realRoot(0) = -(coeff(0))/(coeff(1));

            return realRoot;
        }
        return Eigen::Matrix<double,1,Eigen::Dynamic>();
    }

private:
    Eigen::VectorXd zeroVector = Eigen::VectorXd::Zero(16);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeff;



    //Algorithm to find real roots from companion matrix
    template <typename Derived> Eigen::Matrix<double, 1, Eigen::Dynamic> EigenRealRoots(const Eigen::MatrixBase<Derived>& companionMatrix){
        auto roots = companionMatrix.eigenvalues().eval();
        Eigen::Matrix<double, 1, Eigen::Dynamic> realRoots;
        realRoots.resize(Eigen::NoChange, (roots.imag().array() == 0).template cast<int>().sum());
        int currentRootIndex = 0;
        for(size_t i{}; i != roots.size(); ++i){
            if( roots(i).imag() == 0 ){
                realRoots(currentRootIndex++) = roots(i).real();
            }
        }
        return realRoots;
    }

};


#endif //POLYNOMIAL_COMPANIONMATRIX_H
