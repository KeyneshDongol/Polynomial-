//
// Created by Keynesh Dongol on 13/8/23.
//

#ifndef POLYNOMIAL_SPARSEFUNCTION_H
#define POLYNOMIAL_SPARSEFUNCTION_H

#include "SparseFunctionStruct.h"
#include <stdio.h>
#include <vector>
#include <iostream>
#include <algorithm>


namespace SparseFunction {

    template <int order, int nVar, template<int> typename basisType>
    struct  polynomial : public Internal::sparseFunctionStruct<order,order,nVar,basisType> {

//~~~~ Constants & type definitions
        using parentType = Internal::sparseFunctionStruct<order,order,nVar,basisType>;
        constexpr static int           numBasisFunction(){
            int n = 0;
            return parentType::numBasisFunction(n);
        }

//~~~~


//~~~~ Construction/Assignment
        std::vector<double*>    coeffVectorPtr(){
            std::vector<double*> pointerVector;
            pointerVector.resize(numBasisFunction());
            int counter = 0;
            parentType::coeffVectorPtr(pointerVector,counter);
            return pointerVector;
        }
        template <typename Derived1, typename Derived2>
        double                  fitData(const Eigen::ArrayBase<Derived1>& x, const Eigen::ArrayBase<Derived2>& y){
            assert(x.cols() == nVar);
            assert(x.rows() == y.rows());
            assert(y.cols() == 1);
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> psyx = applyAllBasisFunctions(x).matrix();
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> psyxTy = psyx.transpose()*y.matrix();
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> psyTPsy = psyx.transpose()*psyx;

            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> solution = psyTPsy.llt().solve(psyxTy);
            auto coeffPoint = coeffVectorPtr();
            for (size_t i{}; i != coeffPoint.size(); ++i) {*coeffPoint[i] = solution(i);}

            return ((operator()(x)-y)/(y.maxCoeff() - y.minCoeff())).matrix().squaredNorm()/y.size();
        }

        //~~~~ Evaluation
        template <typename Derived>
        auto                    operator()(const Eigen::ArrayBase<Derived>& x){
            assert(x.cols() == nVar);
            return parentType::operator()(x);
        }


        template <typename Derived>
        static auto             applyAllBasisFunctions(const Eigen::ArrayBase<Derived>& x){
            assert(x.cols() == nVar);
            using returnType = Eigen::Array<typename Eigen::ArrayBase<Derived>::Scalar, Eigen::ArrayBase<Derived>::RowsAtCompileTime, Eigen::Dynamic>;
            returnType result ;
            result.resize(x.rows(), numBasisFunction());
            int currentBasisFunction = 0;
            parentType::applyAllBasisFunctions(x, result, currentBasisFunction);
            return result;
        }

// **************************
// *** Substitution
// **************************

//        template <typename DerivedA, typename DerivedB>
//        auto substitute (const Eigen::ArrayBase<DerivedA>& A, const Eigen::ArrayBase<DerivedB>& B, polynomial<order, nVar, basisType>& oriFunc){
//            assert(A.rows() == B.rows() && B.rows() == nVar); assert(A.cols() == B.cols() && B.cols() == 1);
//
//            Eigen::Array<double, 1, order+1> wValues = Eigen::Array<double, 1, order+1>::LinSpaced(order+1, 0.0, 1.0);
//            auto bigXValue = (A.matrix() * wValues.matrix()).array().colwise()+ B;
//
//            polynomial<order, 1, basisType> smallPoly;
//            return smallPoly.fitData(wValues, oriFunc(bigXValue.transpose()));
//        }


        //~~~~ I/O
        void                    printCoeff(){
            std::vector<int> orders{};
            orders.resize(nVar);
            parentType::printCoeff(orders);
        }

    };



//    template <typename DerivedA, typename DerivedB>
//    static sparseFunction<order, 1, basisType> polySubstitution (const Eigen::ArrayBase<DerivedA>& A, const Eigen::ArrayBase<DerivedA>& B, sparseFunction oriFunc){
//        assert(A.rows() == 1 && B.rows() ==1); assert(A.cols() == nVar && B.cols() == nVar);
//
//        Eigen::Array<double, 1, order+1> wValues = Eigen::Array<double, 1, order+1>::LinSpaced(order+1, 0.0, 1.0);
//        auto bigXValue = (A.matrix() * wValues.matrix()).array().colwise()+ B;
//
//        sparseFunction<order, 1, basisType> smallPoly;
//
//        auto result = smallPoly.fitData(wValues, oriFunc(bigXValue));
//        return result;
//    }

} // namespace SparseFunction



#endif //POLYNOMIAL_SPARSEFUNCTION_H
