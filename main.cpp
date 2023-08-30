//
//  main.cpp
//  PolyNested
//
//  Created by Keynesh Dongol on 20/8/23.
//

//#include <iostream>
////#include "HornerMethod.h"
//#include <Eigen/Dense>
//#include <chrono>
//#include "Vector.h"
//#include "Math.h"
//
//
////template<typename Derived>
//
//
//
//
//int main(int argc, const char * argv[]) {
//
//    const int Order = 3;
//    const int Ndim = 1;
//
////    std::vector<double> testArray = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(5, 1);
////    std::vector<double> testArray = {2,-6,2,-1};
//
//
//    auto t1 = std::chrono::high_resolution_clock::now();
//    std::vector<double> v{1,2,3,4,5};
//    std::vector<double> v2{1,2,3,4,5};
//
////    math::expression<>
//
//    //    someVec =
//
////    Vector<double> someVec(v);
////    Vector<double> someOtherVec(v2);
////    Vector<double> vecSum = v + v2;
////    for (auto i: vecSum){
////    std:;cout << vecSum[i] << "\n";
////    }
////    std::cout << vecSum;
////    someVec
//
////    Polynomial<Order, Ndim> bla(testArray); ///Seems that the construction of this variable took no time.
////    Polynomial<std::vector<double>, double, Order, Ndim> bla;
////    std::cout << bla.BinomialCoefficient() <<"\n\n";
////    std::cout << "The evaluated polynomial is: " << bla.polyHorner(3)<<"\n\n";
//
//
////    auto temp = polyHorner(bla, 1);
////    std::cout << temp << "\n\n";
//////
////
//
//    auto t2 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double, std::nano> fp_ms = t2 - t1;
////
////    std::move(bla).BinomialCoefficent();
//////
//////
//    std::cout << "Time taken: " << fp_ms.count() << " nanoseconds \n\n" << std::endl;
////
//    return 0;
//}
# include <chrono>
# include <iomanip>
# include <iostream>
# include "Vector.h"
# include "Timer.h"
# include <boost/core/demangle.hpp>
# include <Eigen/Dense>
//# include "SparseFunction .h"
# include "SparseFunction.h"
# include "BasisFunctionsSets.h"



//using namespace std::cout;
using namespace SparseFunction;



int main(){

    const int Order = 4;
    const int Ndim = 3;
    using dtype = double;
    polynomial<Order, Ndim, TaylorBasis> testFunc;

    std::cout << "Number of basis functions: " << testFunc.numBasisFunction() << "\n\n"; ///Returns int
    constexpr int size = testFunc.numBasisFunction();


    Timer T,T2,T3;
//    T.set();

    std::vector<dtype> _vec1(size);
    std::vector<dtype> _vec2(size);
    std::vector<dtype> _vec3(size);

    // Initialize veector's contents.

    Vector<dtype> vec1(std::move(_vec1)); ///Constructor called
    Vector<dtype> vec2(std::move(_vec2));
    Vector<dtype> vec3(std::move(_vec3));

//    unsigned long start_ms_no_ets =
//        std::chrono::duration_cast<std::chrono::milliseconds>
//        (std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "\nNo-ETs evaluation starts.\n";

    T.set();
    Vector<dtype> result_no_ets = vec1 + (vec2*vec3);
    auto time = T.elapsed();

    std::cout << "No-ETs. Time eclapses: " << time << "s.\n\n";


    std::cout << "Evaluation using ETs starts.\n";

//    T2.set();


    expr::terminal<Vector<dtype>> vec4(vec1);
    expr::terminal<Vector<dtype>> vec5(vec2);
    expr::terminal<Vector<dtype>> vec6(vec3);
    T2.set();


    Vector<dtype> result_ets = (vec4 + vec5*vec6); ///This actually returns a expr::terminal<Vector<double>> type. More specifically, expr::binary_ops<expr::vec_plus_t, expr::terminal<Vector<double>>, expr::binary_ops<expr::vec_prod_t, expr::terminal<Vector<double>>, expr::terminal<Vector<double>>>>. but here becasue of the retrun type, it is giving me a Vector<double> whic implies that come conversoi nhis cgoig on tin the background.
    ///Anyways, my point is, check to see how the conversion is being carried out.
    auto time2 = T2.elapsed();
    std::cout << "With ETs. Time eclapses: " << time2
              << " s.\n" << std::endl;

    std::cout << "Evaluation using Eigen's ETs starts.\n";

    T3.set();

/*!
    Eigen::Array<dtype, size, 1> vec7(Eigen::Array<dtype, size,1>::Ones());
    Eigen::Array<dtype, size, 1> vec8(Eigen::Array<dtype, size,1>::Ones());
    Eigen::Array<dtype, size, 1> vec9(Eigen::Array<dtype, size,1>::Ones());
//    Eigen::Map<dtype, size, 1> vec8;
//    Eigen::Map<dtype, size, 1> vec9;
//    expr::terminal<Vector<dtype>> vec4(vec1);
//    expr::terminal<Vector<dtype>> vec5(vec2);
//    expr::terminal<Vector<dtype>> vec6(vec3);

    auto result_ets_eigen = (vec7 + vec8*vec9); ///Running this, we find that it takes the same time using Eigen's templated library as it does above. This is brilliant news


    ///This actually returns a expr::terminal<Vector<double>> type.
    ///More specifically, expr::binary_ops<expr::vec_plus_t, expr::terminal<Vector<double>>, expr::binary_ops<expr::vec_prod_t, expr::terminal<Vector<double>>, expr::terminal<Vector<double>>>>. but here becasue of the retrun type, it is giving me a Vector<double> whic implies that come conversoi nhis cgoig on tin the background.
    ///Anyways, my point is, check to see how the conversion is being carried out.
    auto time3 = T2.elapsed();
    std::cout << "With Eigen's ETs. Time eclapses: " << time2
                << " s.\n" << std::endl;
 */
    auto ets_ret_type = (vec4 + vec5*vec6);
    std::cout << "\nETs result's type:\n-> ";
    std::cout << boost::core::demangle( typeid(decltype(ets_ret_type)).name() ) << "\n\n";

    std::cout << boost::core::demangle(typeid(decltype(result_ets)).name() ) << "\n\n";

//    expr::binary_ops<expr::vec_plus_t, expr::termi``nal<Vector<double>>, expr::binary_ops<expr::vec_prod_t, expr::terminal<Vector<double>>, expr::terminal<Vector<double>>>>;
    ///TODO: Try to include a section detailing how to do an exterior product. This was we can confirm our undersyanding of how the above templated expressions are supposed to be written.
    ///
    ///TODO: Similarly, try and understand how these expression templates can be passed as a "structure". In other words, try to figure out how to include this into Horner's emthod.


    std::cout << "-> Now we start working on the basis function and sparseFunction in general.\n\n";


//

    auto coeffPoint = testFunc.coeffVectorPtr(); ///Returns std::__1::vector<double*, std::__1::allocator<double*>>

    testFunc.printCoeff(); std::cout << "\n\n"; ///Returns void
    std::cout << "\nPrintCoeff's result's type:\n-> ";
    std::cout << boost::core::demangle( typeid(decltype(testFunc.printCoeff())).name() ) << "\n\n";

    const int nDataPoint = 20;
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> v = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Random(nDataPoint, Ndim);
    ///I think we can rewrite the above in a manner that we have our own random matrix function thus keeeping the mathematical operatios "local".

    Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> y = 3.1*v.col(0) + 2.*v.col(1) + 3.*v.col(1)*v.col(1) + 1.2;

//    std::cout << "The fitting data: \n" << testFunc.fitData(v, y) << "\n\n";
//    testFunc.printCoeff();

    ///Now we evaluate the polynomials
    std::cout << "testFunc(v) value: \n" << testFunc(v);
    std::cout << "\nAnd the fitting error is: " << testFunc.fitData(v, y)<< "\n"; //This is giving me some error.
//    std::cout << "\noperator overlaoded type is:\n->";
//    std::cout << "\nEvaluated sparseFunc is: \n->"<< boost::core::demangle(typeid(decltype(testFunc(v))).name()) << "\n\n";
    std::cout << boost::core::demangle(typeid(decltype(testFunc.fitData(v, y))).name()) << "\n\n";

    std::cout << "\ntestFunc.printCoeff() type ->" << boost::core::demangle(typeid(decltype(testFunc.printCoeff())).name()) << "\n\n"; //void
    std::cout << "\ntestFunc type ->"<< boost::core::demangle(typeid(decltype(testFunc)).name()) << "\n\n"; //SparseFunction::sparseFunction<4,2,SparseFunction::TaylorBasis>





// **************************
// *** Substitution
// **************************
//    std::cout << "-\n\n\n>Now we will see how the substitution is done:\n\n";
//
//    Eigen::Array<double, order+1, 1> wValues = Eigen::Array<double, order+1, 1>::LinSpaced(order+1, 0.0, 1.0);
//    Eigen::Array<double, numVar, 1> a = Eigen::Array<double, numVar, 1>::Random();
//    Eigen::Array<double, numVar, 1> b = Eigen::Array<double, numVar, 1>::Random();
//
//    std::cout << "\n -> wValues are: \n" << wValues.transpose() << "\n\n";
//    std::cout << "\n -> Value of vector A: \n" << a.transpose() << "\n\n";
//    std::cout << "\n -> Value of vector B: \n" << b.transpose() << "\n\n";
//
//    std::cout << (a.matrix() * wValues.matrix().transpose()).array().colwise() + b << "\n\n";
//
//    auto bigXValue = (a.matrix() * wValues.matrix().transpose()).array().colwise()+ b;
//    std::cout << "\n Polynomial evaluated on these \"points\": \n" << bla(bigXValue.transpose()) << "\n\n\n";
//
//    sparseFunction<order, 1, TaylorBasis> smallPoly;
//    std::cout << "This is the new 1D vector's fitting error: \n" << smallPoly.fitData(wValues, bla(bigXValue.transpose())) << "\n\n\n";
//    smallPoly.printCoeff();
//
//    std::cout << "\n\n\n";

//    Eigen::Array<double, 1, numVar> a = Eigen::Array<double, 1, numVar>::Random();
//    Eigen::Array<double, 1, numVar> b = Eigen::Array<double, 1, numVar>::Random();
//    auto temp = bla.polySubstitution(a, b, bla);
//    std::cout << temp << "\n\n";

//    SparseFunction::sparseFunction<order, 1, TaylorBasis> test;
//    auto _new = bla.polySubstitution(a, b);
//    std::cout << _new << "\n\n";
    //What we need to do is, figute out a way to call a funnction




    return 0;
}
