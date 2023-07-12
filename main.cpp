//
//  main.cpp
//  Polynoimial Operator
//
//  Created by Keynesh on 15/2/23.
//

//#define EIGEN_USE_BLAS

#include <iostream>
#include "Polynomial.hpp"
//#include <Eigen/Dense>
//#include <typeinfo>
#include <chrono>
#include <vector>
#include <numeric>
#include <utility>

using namespace std::chrono_literals;


//auto leftRotated(Eigen::Matrix<double, 3, 1> row) {
//    Eigen::Matrix<double,3,1> newMatrix ({{row(1)},{row(2)},{row(0)}});
//    return newMatrix;
//}
void print(std::vector<int> const &input)
{
    std::copy(input.begin(),
            input.end(),
            std::ostream_iterator<int>(std::cout, " "));
}

int main(int argc, const char * argv[]) {
    // insert code here...

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    


//    while (true)
//    {
//        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> a;
//        a.resize(2000, 2000);
//        a.setRandom();
//        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> b;
//        b.resize(2000, 2000);
//        b.setRandom();
//        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> c = a * b;
//        std::cout << c(0,0) << std::endl;
//    }
//    return 0;
    int n = 5;
    std::vector<double> myList;
     
    for (int i=0; i < n ; i++) {
        //Start time
        auto t1 = high_resolution_clock::now();
        Eigen::Matrix<double,35,1> coeff = Eigen::Matrix<double,35,1>::Random();
//            Eigen::Matrix<double, 35, 1> coeff  = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,-19,20,21,22,23,-24,25,26,0,28,29,30,31,32,33,34,35};
//            Eigen::Matrix<double, 15, 1> coeff2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        //    Eigen::Matrix<double, 15, 1> coeff2 = Eigen::Matrix<double, 15, 1>::Ones();
        //    Eigen::Matrix<double, 2, 1> points = Eigen::Matrix<double, 2, 1>::Ones();
        
        //
        //    Eigen::Matrix<double, 3, 50> random = Eigen::Matrix<double, 3, 50>::Random();
        const int numVar = 4;
        const int nDataPoints = 5;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> random = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Random(numVar, nDataPoints);
        std::cout << random << "\n";
        Polynomial<3> poly33(coeff);
        //    std::cout << poly33 <<"\n";
        
        //        std::cout << poly33.evalAll(random) << "\n";
        Eigen::Matrix<double, 1, Eigen::Dynamic> bla = poly33.evalAll(random);
        Eigen::Matrix<double, 3, 1> A_,a_;
        A_ = Eigen::Matrix<double, 3, 1>::Ones();
        a_ = Eigen::Matrix<double, 3, 1>::Ones();
        std::cout << poly33.reductionOperator(A_, a_)<< "\n";
        //
        ////
        //    Eigen::Matrix<double, 2, 1> random2 = Eigen::Matrix<double, 2, 1>::Ones();
        //    Polynomial<2> poly22(coeff2);
        //    std::cout << poly22.evalAll(random2) << "\n";
        //    Eigen::Matrix<double, 2, 1> A_2,a_2;
        //    A_2 = Eigen::Matrix<double, 2, 1>::Zero();
        //    a_2 = Eigen::Matrix<double, 2, 1>::Ones();
        //    std::cout << poly22.reductionOperator(A_2, a_2)<< "\n";
        //
        //    std::cout << (poly22.reductionOperator(A_2, a_2)).findRealRoots()<< "\n";
        
        
        
            Eigen::Matrix<double,5,1> findroots = {120, -154, 71, -14, 1};
            Eigen::Matrix<double,5,1> findroots1 = {4, -5, 0, -2,1};
            Eigen::Matrix<double,5,1> findroots2 = {4, -5, 0, -2, 0};
            Eigen::Matrix<double,5,1> findroots3 = {-1, 0, 0, 0, 0};
            Eigen::Matrix<double,5,1> findroots4 = {-6, 11, -6, 1,0};
            Eigen::Matrix<double,5,1> findroots5 = {2, -3, 1, 0,0};
            Eigen::Matrix<double,5,1> findroots6 = {-24, 1, 0, 0, 0};
        //
            Polynomial<1> poly1(findroots);
            Polynomial<1> poly2(findroots1);
            Polynomial<1> poly3(findroots2);
            Polynomial<1> poly4(findroots3);
            Polynomial<1> poly5(findroots4);
            Polynomial<1> poly6(findroots5);
            Polynomial<1> poly7(findroots6);
            std::cout << poly1.findRealRoots() << "\n";
            std::cout << poly2.findRealRoots() << "\n";
            std::cout << poly3.findRealRoots() << "\n";
            std::cout << poly4.findRealRoots() << "\n";
            std::cout << poly5.findRealRoots() << "\n";
            std::cout << poly6.findRealRoots() << "\n";
            std::cout << poly7.findRealRoots() << "\n";
        
        
        //    auto stop = high_resolution_clock::now();
        //
        //    std::cout << "Hello, World!\n";
        //    auto duration = duration_cast<microseconds>(stop - start);
        //
        //    std::cout << "Time taken by function: "
        //             << duration.count() << " microseconds plus the build time is somewhere else(Check in build log)" << "\n";
        //
        //
        
        //End time
        auto t2 = high_resolution_clock::now();
        
        //     Getting number of milliseconds as an integer.
        //    auto ms_int = duration_cast<milliseconds>(t2 - t1);
        
        //     Getting number of milliseconds as a double.
        duration<double, std::milli> ms_double = t2 - t1;
        
        
        myList.push_back(ms_double.count());
        
    }

    for (std::vector<int>::size_type i = 0; i < myList.size(); i++) {
        std::cout << myList.at(i) << ' '<< "\n";
    }
    auto result = std::reduce(myList.begin(), myList.end())/n;
    std::cout << "The average time to run the script is :" << result << "\n";
//    std::cout << ms_int.count() << "ms\n";
//    std::cout << ms_double.count() << "ms\n";
    
    return 0;
}



