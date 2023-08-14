//
// Created by Keynesh Dongol on 13/8/23.
//

#ifndef POLYNOMIAL_BASISFUNCTIONSSETS_H
#define POLYNOMIAL_BASISFUNCTIONSSETS_H


#include <Eigen/Dense>

namespace SparseFunction {

    template <int order> struct TaylorBasis{
        template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
            return x * TaylorBasis<order-1>::basisF(x);
        }
    };
    template <> struct TaylorBasis<0>{
        template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
            return Eigen::ArrayBase<Derived>::Constant(x.rows(), x.cols(), 1);
        }
    };
    template <> struct TaylorBasis<1>{
        template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
            return x.matrix().array();
        }
    };

} // namespace SparseFunction




#endif //POLYNOMIAL_BASISFUNCTIONSSETS_H
