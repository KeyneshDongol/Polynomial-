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
    struct sparseFunction : public Internal::sparseFunctionStruct<order,order,nVar,basisType> {

//~~~~ Constants & type definitions
        using parentType = Internal::sparseFunctionStruct<order,order,nVar,basisType>;
        constexpr static int           numBasisFunction(){
            int n = 0;
            return parentType::numBasisFunction(n);
        }

//~~~~


//~~~~ Construction/Assignment
        std::vector<double*>    coeffVectorPtr(){
            std::vector<double*> pointerVector;//(numBasisFunction());
//            std::vector<std::unique_ptr<double>> pointerVector(numBasisFunction());

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

//        #pragma clang loop vectorize(enable)
//        #pragma clang loop interleave(enable)

            for (size_t i{}; i != coeffPoint.size(); ++i) {*coeffPoint[i] = solution(i);}

            return ((operator()(x)-y)/(y.maxCoeff() - y.minCoeff())).matrix().squaredNorm()/y.size();
        }

        //~~~~ Evaluation
        template <typename Derived>
        auto                    operator()(const Eigen::ArrayBase<Derived>& x){
            assert(x.cols() == nVar);
            return parentType::operator()(x);
        }
//        template <typename Derived>
//        static auto             applyAllBasisFunctionsVerySlow(const Eigen::ArrayBase<Derived>& x){
//            assert(x.cols() == nVar);
//            using returnType = Eigen::Array<typename Eigen::ArrayBase<Derived>::Scalar, Eigen::ArrayBase<Derived>::RowsAtCompileTime, Eigen::Dynamic>;
//            returnType result = returnType::Constant(x.rows(), numBasisFunction(), 1.);
//            sparseFunction<order,nVar,basisType> temp;
//            auto coeffPoint = temp.coeffVectorPtr();
//
//            for (auto p : coeffPoint) {*p = 0;}
//            for (int i=0; i<numBasisFunction(); i++){
//                *coeffPoint[i] = 1.;
//                result.col(i) = temp(x);
//                *coeffPoint[i] = 0.;
//            }
//            return result;
//        }
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

        //~~~~ I/O
        void                    printCoeff(){
            std::vector<int> orders{};
            orders.resize(nVar);
            parentType::printCoeff(orders);
        }

    };

} // namespace SparseFunction



#endif //POLYNOMIAL_SPARSEFUNCTION_H
