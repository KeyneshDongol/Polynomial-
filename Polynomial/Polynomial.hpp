//
//  Polynomial.hpp
//  Polynoimial Operator
//
//  Created by Keynesh on 15/2/23.
//

#ifndef Polynomial_h
#define Polynomial_h


#include <iostream>
#include <Eigen/Dense>
#include <cassert>

template<int ndim> class Polynomial;


//############################
//Polynomial of 1 dimension
//############################


template<> class Polynomial<1> {
// Polynomial of degree 4 of order 1
public:
    template <typename Derived>  explicit Polynomial(const Eigen::DenseBase<Derived>& coeffpass): coeff(coeffpass) { // CoeffPass should be a 5x1 Eigen Vector
        assert(coeffpass.rows() == 5);  assert( coeffpass.cols() == 1);
    }
    Polynomial(int degree){
        assert(degree>=0 & degree <5);
        coeff = Eigen::Matrix <double, 5, 1>::Zero();
        coeff(degree) = 1.;
    }
    
    //********************************
    //* Arithmetic
    //********************************
    Polynomial<1>& operator+=(const Polynomial<1>& other){ coeff += other.coeff; return *this; }
    Polynomial<1>& operator-=(const Polynomial<1>& other){ coeff -= other.coeff; return *this; }
    Polynomial<1>& operator+=(const double scalar){ coeff(0) += scalar; return *this; }
    Polynomial<1>& operator-=(const double scalar){ coeff(0) -= scalar; return *this; }
    Polynomial<1>& operator*=(const double scalar){ coeff *= scalar; return *this; }
    Polynomial<1>& operator/=(const double scalar){ coeff /= scalar; return *this; }

    friend std::ostream &operator<<(std::ostream &t_os, Polynomial<1> const& t_poly);
    
// =============================
// IMPLEMENTATION
// =============================
private:
    // The polynomial is expressed as
    // 0th order a0, 1th order a1, etc
    // a0, a1, etc are scalars
//Data
    Eigen::Matrix<double, 5, 1> coeff; // Polyn coefficients of increasing degree
//    double b0, b1, b2, b3, b4;

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    friend Polynomial<1> operator+(Polynomial<1> lhs, const Polynomial<1>& rhs){ lhs += rhs; return lhs;}
    friend Polynomial<1> operator-(Polynomial<1> lhs, const Polynomial<1>& rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<1> operator+(Polynomial<1> lhs, const double rhs){ lhs += rhs; return lhs;}
    friend Polynomial<1> operator+(const double lhs, Polynomial<1> rhs){ rhs += lhs; return rhs;}
    friend Polynomial<1> operator-(Polynomial<1> lhs, const double rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<1> operator-(const double lhs, Polynomial<1> rhs){ rhs -= lhs; rhs *=  -1.; return rhs;}
    friend Polynomial<1> operator*(Polynomial<1> lhs, const double rhs){ lhs *= rhs; return lhs;}
    friend Polynomial<1> operator*(const double lhs, Polynomial<1> rhs){ rhs *= lhs; return rhs;}
    friend Polynomial<1> operator/(Polynomial<1> lhs, const double rhs){ lhs /= rhs; return lhs;}
    

    template <typename Derived> auto evalAll(const Eigen::DenseBase<Derived>& points) const {
        assert(points.rows() == 1);

        return (coeff(0) +
                coeff(1) * points.array() +
                coeff(2) * points.array()* points.array() +
                coeff(3) * points.array()* points.array()* points.array() +
                coeff(4) * points.array()* points.array()* points.array()* points.array()).matrix();
    }

private:
    //Algorithm to find real roots from companion matrix
    template <typename Derived> Eigen::Matrix<double, 1, Eigen::Dynamic> EigenRealRoots(const Eigen::MatrixBase<Derived>& companionMatrix){
        auto roots = companionMatrix.eigenvalues().eval();
        Eigen::Matrix<double, 1, Eigen::Dynamic> realRoots;
        realRoots.resize(Eigen::NoChange, (roots.imag().array() == 0).template cast<int>().sum());
        int currentRootIndex = 0;
        for(int i=0; i< roots.size(); ++i){
            if( roots(i).imag() == 0 ){
                realRoots(currentRootIndex++) = roots(i).real();
            }
        }
        return realRoots;
    }
    
    const Eigen::VectorXd zeroVector = Eigen::VectorXd::Zero(16);


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
};

const char PM(double number){
    if(number >= 0){
        char p = '+';
        return p;
    }else {
        char m = '-';
        
        return m;
    }
}

 std::ostream &operator<<(std::ostream &t_os, Polynomial<1> const& t_poly){
     t_os << t_poly.coeff(0) << " " << PM(t_poly.coeff(1)) << t_poly.coeff(1) << "x " << PM(t_poly.coeff(2)) << t_poly.coeff(2) << "x^2 " << PM(t_poly.coeff(3)) << t_poly.coeff(3) << "x^3 " << PM(t_poly.coeff(4)) << t_poly.coeff(4) << "x^4\n";
        
     return t_os;
 }


//############################
// Polynomial of 2 dimension
//############################

template<> class
Polynomial<2> {//TODO: Operator overload for 2D case plus write out a better way to cout the monomials
// Polynomial of degree 4 of order 2
    
    // The polynomial is linearly expressed as
    // 0th order     a0 +
    // 1st order     a1 * P +
    // 2nd order     P.transpose * a2 * P [a2 is an Upper Triangle]
    // 3rd order     P^2.transpose * a3 * P
    // 4th order     P^3.transpose * a4 * P + a4m * (x * y)^2

public:
    
   template <typename Derived>
    explicit Polynomial(const Eigen::DenseBase<Derived>& coeffpass): coeff(coeffpass) {// CoeffPass should be a 15x1 Eigen Vector.
       assert(coeffpass.rows() == 15);  assert( coeffpass.cols() == 1);
    }// CoeffPass should be a 15x1 Eigen Vector

    Polynomial(int degreeX, int degreeY);
    
    
private:
    //********************************
    //* Arithmetic
    //********************************
    Polynomial<2>& operator+=(const Polynomial<2>& other){ coeff += other.coeff; return *this; }
    Polynomial<2>& operator-=(const Polynomial<2>& other){ coeff -= other.coeff; return *this; }
    Polynomial<2>& operator+=(const double scalar){ coeff(0) += scalar; return *this; }
    Polynomial<2>& operator-=(const double scalar){ coeff(0) -= scalar; return *this; }
    Polynomial<2>& operator*=(const double scalar){ coeff *= scalar; return *this; }
    Polynomial<2>& operator/=(const double scalar){ coeff /= scalar; return *this; }

    friend std::ostream &operator<<(std::ostream &t_os, Polynomial<2> const& t_poly);

// =============================
// IMPLEMENTATION
// =============================
private:
    // The polynomial is expressed as
    // 0th order a0, 1th order a1, etc
    // a0, a1, etc are scalars
//Data
    Eigen::Matrix<double, 15, 1> coeff;// Polyn coefficients of increasing degree
    
    //****************************************
    // -> coeff vector will be rehaped to smaller sub matrixes to carry out linear algebriac expression
    // a0 is scaler
    // a1 is [1,Ndim]
    // a2 is [Ndim, Ndim] Upper Triangular
    // a3 is [Ndim, Ndim]
    // a4 is [Ndim,Ndim]
    // a4m is scalar
    //****************************************
    
public:

   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    friend Polynomial<2> operator+(Polynomial<2> lhs, const Polynomial<2>& rhs){ lhs += rhs; return lhs;}
    friend Polynomial<2> operator-(Polynomial<2> lhs, const Polynomial<2>& rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<2> operator+(Polynomial<2> lhs, const double rhs){ lhs += rhs; return lhs;}
    friend Polynomial<2> operator+(const double lhs, Polynomial<2> rhs){ rhs += lhs; return rhs;}
    friend Polynomial<2> operator-(Polynomial<2> lhs, const double rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<2> operator-(const double lhs, Polynomial<2> rhs){ rhs -= lhs; rhs *=  -1.; return rhs;}
    friend Polynomial<2> operator*(Polynomial<2> lhs, const double rhs){ lhs *= rhs; return lhs;}
    friend Polynomial<2> operator*(const double lhs, Polynomial<2> rhs){ rhs *= lhs; return rhs;}
    friend Polynomial<2> operator/(Polynomial<2> lhs, const double rhs){ lhs /= rhs; return lhs;}
    

   //values are respective points
    
    template <typename Derived>Eigen::Matrix<double, 1, Eigen::Dynamic> evalAll(const Eigen::MatrixBase<Derived>& points) { //points is 2xD Eigen Matrix
        assert(points.rows() == 2);
        
        const double a0 = coeff(0);
        const Eigen::Map<Eigen::RowVector2d> a1 (coeff.data()+1);
        const Eigen::Map<Eigen::RowVector2d> a2_p1 (coeff.data()+3);
        const double a2_p2 = coeff(5);
        const Eigen::Map<Eigen::Matrix2d> a3 (coeff.data()+6);
        const Eigen::Map<Eigen::Matrix2d> a4 (coeff.data()+10);
        const double a4m = coeff(14);
        
        auto xyProd = (points.template topRows(1).array() * points.template topRows(1).array() * points.template bottomRows(1).array() * points.template bottomRows(1).array()).matrix();

        auto value2 = a1 * points;
        auto value3 = (points.transpose().template leftCols<1>() * (a2_p1 * points)).diagonal().transpose();
        auto value4 = (points.transpose().template rightCols<1>() * (a2_p2* points.template bottomRows<1>())).diagonal().transpose();
        auto value5 = (((points.array() * points.array()).matrix()).transpose() * a3 * points).diagonal().transpose();
        auto value6 = (((points.array() *points.array()*points.array()).matrix()).transpose() * a4 * points).diagonal().transpose();
        auto value7 = a4m * xyProd;

        return (a0 + (value2 + value3 + value4 + value5 + value6 + value7).array()).matrix();
    }

    template <typename Derived, typename DerivedOther> Polynomial<1> reductionOperator(const Eigen::MatrixBase<Derived>& A_, const Eigen::MatrixBase<DerivedOther>& a_){
        //A_,a_ are both 2x1 Eigen Column Vector
        assert(A_.rows() && a_.rows() == 2); assert(A_.cols() && a_.cols() == 1);

        const double a0 = coeff(0);
        const Eigen::Map<Eigen::RowVector2d> a1 (coeff.data()+1);
        const Eigen::Map<Eigen::RowVector2d> a2_p1 (coeff.data()+3);
        const double a2_p2 = coeff(5);
        const Eigen::Map<Eigen::Matrix2d> a3 (coeff.data()+6);
        const Eigen::Map<Eigen::Matrix2d> a4 (coeff.data()+10);
        const double a4m = coeff(14);


        auto                            A = A_(0);
        auto                            B = A_(1);
        auto                            a = a_(0);
        auto                            b = a_(1);
        auto                            Atwo = (A_.array()*A_.array()).matrix();
        auto                            atwo = (a_.array()*a_.array()).matrix();
        auto                            Athree = (A_.array()*A_.array()*A_.array()).matrix();
        auto                            athree = (a_.array()*a_.array()*a_.array()).matrix();



// Eigen expressions -> double
        auto w_0 = a0  + (a1*a_) + a4m*a*a*b*b +
                   (a_.transpose().template head<1>() * (a2_p1 * a_)).value()
                        + (atwo.transpose()*a3*a_).value()
                        + (athree.transpose()*a4*a_).value();

// Eigen expressions -> double
        auto w_1 = (a1*A_).value()+
                    ((a2_p1 * a_) * A_.transpose().template head<1>()).value() +
                   (A_.transpose().template tail<1>() * (a2_p2 * a_.template tail<1>())).value() +
                   (a_.transpose().template tail<1>() * (a2_p2 * A_.template tail<1>())).value()
                        + (atwo.transpose()*a3*A_).value()
                        + ((2*A_.array()*a_.array()).matrix().transpose()*a3*a_).value()
                        + ((3*A_.array()*atwo.array()).matrix().transpose()*a4*a_).value()
                        + ((athree).transpose()*a4*A_).value()
                            + (a4m * (2*A*B*B*a + 2*B*b*a*a + 2*A*A*B*b));
// Eigen expressions -> double
        auto w_2 = ((a2_p1 * A_)*A_.transpose().template head<1>()).value() +
                    (A_.transpose().template tail<1>() * (a2_p2 * A_.template tail<1>())).value()
                        + ((2*A_.array()*a_.array()).matrix().transpose()*a3*A_).value()
                        + (Atwo.transpose()*a3*a_).value()
                        + (3*(Atwo.array()*a_.array()).matrix().transpose()*a4*a_).value()
                        + (3*(A_.array()*atwo.array()).matrix().transpose()*a4*A_).value()
                            + a4m*(A*A*b*b + 4*A*B*a*b + B*B*a*a);

// Eigen expressions -> double
        auto w_3 = (Atwo.transpose()*a3*A_).value()
                        + (Athree.transpose()*a4*a_).value()
                        + (3*(Atwo.array()*a_.array()).matrix().transpose()*a4*A_).value()
                                + a4m*2*A*B*B*a;

// Eigen expressions -> double
        auto w_4 = (Athree.transpose()*a4*A_).value()
                          + a4m*A*A*B*B;

        //returns a 5x1 Eigen Row vector which will be passed as coefficients for Poly<1> Class.
        return Polynomial<1>(Eigen::Matrix <double, 5, 1>{w_0 , w_1 , w_2 , w_3 , w_4});
    }

    
    
};

std::ostream &operator<<(std::ostream &t_os, Polynomial<2> const& t_poly){
    t_os << t_poly.coeff(0) << " " << PM(t_poly.coeff(1)) << t_poly.coeff(1) << "x " << PM(t_poly.coeff(2)) << t_poly.coeff(2) << "y "<< PM(t_poly.coeff(3)) << t_poly.coeff(3) << "x^2 "<< PM(t_poly.coeff(4)) << t_poly.coeff(4) << "xy "<< PM(t_poly.coeff(5)) << t_poly.coeff(5) << "y^2 "<< PM(t_poly.coeff(6)) << t_poly.coeff(6) << "x^3 "<< PM(t_poly.coeff(7)) << t_poly.coeff(7) << "y^2x "<< PM(t_poly.coeff(8)) << t_poly.coeff(8) << "x^2y "<< PM(t_poly.coeff(9)) << t_poly.coeff(9) << "y^3 "<< PM(t_poly.coeff(10)) << t_poly.coeff(10) << "x^4 "<< PM(t_poly.coeff(11)) << t_poly.coeff(11) << "y^3x "<< PM(t_poly.coeff(12)) << t_poly.coeff(12) << "x^3y "<< PM(t_poly.coeff(13)) << t_poly.coeff(13) << "y^4 "<< PM(t_poly.coeff(14)) << t_poly.coeff(14) << "x^2y^2 \n";
       
    return t_os;
}



//############################
//Polynomial of 3 dimension
//############################


template<> class Polynomial<3> {
public:
    // Polynomial of degree 4 of order 3
    // The polynomial is expressed as
    // 0th degree   c0 +
    // 1st degree   c1 * P +
    // 2nd degree   P.transpose * c2 * P +   //c2 is an upper triangular matrix
    // 3rd degree   P^2.transpose * c3 * P + c3m * x * y * z
    // 4th degree   P^3.transpose * c4 * P + c4m1 * P * x * y * z + c4m2 * (P * RotateLeft(P))^2

    template <typename Derived>
    explicit Polynomial(const Eigen::DenseBase<Derived>& coeffpass): coeff(coeffpass) {// CoeffPass should be a 35x1 Eigen Vector.
        assert(coeffpass.rows() == 35);  assert( coeffpass.cols() == 1);
     }// CoeffPass should be a 35x1 Eigen Vector.

    Polynomial(int degreeX, int degreeY, int degreeZ);

//********************************
//* Arithmetic
//********************************
    Polynomial<3>& operator+=(const Polynomial<3>& other){ coeff += other.coeff; return *this; }
    Polynomial<3>& operator-=(const Polynomial<3>& other){ coeff -= other.coeff; return *this; }
    Polynomial<3>& operator+=(const double scalar){ coeff(0) += scalar; return *this; }
    Polynomial<3>& operator-=(const double scalar){ coeff(0) -= scalar; return *this; }
    Polynomial<3>& operator*=(const double scalar){ coeff *= scalar; return *this; }
    Polynomial<3>& operator/=(const double scalar){ coeff /= scalar; return *this; }

    friend std::ostream &operator<<(std::ostream &t_os, Polynomial<3> const& t_poly);

// =============================
// IMPLEMENTATION
// =============================

private:

// The polynomial is expressed as
// 0th order a0, 1th order a1, etc
// a0, a1, etc are scalars

//Data
    Eigen::Matrix<double, 35, 1> coeff;// Polyn coefficients of increasing degree
    
    //****************************************
    // coeff vector will be rehaped to smaller sub matrixes to carry out linear algebriac expression
    // c0 is scalar
    // c1 is [1,Ndim]
    // c2 is [Ndim,Ndim]  Upper Triangular
    // c3 is [Ndim,Ndim]
    // c3m is scalar
    // c4 is [Ndim,Ndim]
    // c4m1 is [1,Ndim]
    // c4m2 is [1,Ndim]
    //****************************************

    
public:

//Data
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    friend Polynomial<3> operator+(Polynomial<3> lhs, const Polynomial<3>& rhs){ lhs += rhs; return lhs;}
    friend Polynomial<3> operator-(Polynomial<3> lhs, const Polynomial<3>& rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<3> operator+(Polynomial<3> lhs, const double rhs){ lhs += rhs; return lhs;}
    friend Polynomial<3> operator+(const double lhs, Polynomial<3> rhs){ rhs += lhs; return rhs;}
    friend Polynomial<3> operator-(Polynomial<3> lhs, const double rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<3> operator-(const double lhs, Polynomial<3> rhs){ rhs -= lhs; rhs *=  -1.; return rhs;}
    friend Polynomial<3> operator*(Polynomial<3> lhs, const double rhs){ lhs *= rhs; return lhs;}
    friend Polynomial<3> operator*(const double lhs, Polynomial<3> rhs){ rhs *= lhs; return rhs;}
    friend Polynomial<3> operator/(Polynomial<3> lhs, const double rhs){ lhs /= rhs; return lhs;}
//
//
//   //values at respective points
//
//
//
    template <typename Derived>Eigen::Matrix<double, 1, Eigen::Dynamic> evalAll(const Eigen::MatrixBase<Derived>& points) { //points is 3xD Eigen Matrix
        assert(points.rows() == 3);
       
        const double                                    c0      = coeff(0);
        const Eigen::Map<const Eigen::RowVector3d>      c1      (coeff.data()+1);
        const Eigen::Map<const Eigen::RowVector3d>      c2_p1   (coeff.data()+4);
        const Eigen::Map<const Eigen::RowVector2d>      c2_p2   (coeff.data()+7);
        const double                                    c2_p3   = coeff(9);
        const Eigen::Map<const Eigen::Matrix3d>         c3      (coeff.data()+10);
        const double                                    c3m     = coeff(19);
        const Eigen::Map<const Eigen::Matrix3d>         c4      (coeff.data()+20);
        const Eigen::Map<const Eigen::RowVector3d>      c4m1    (coeff.data()+29);
        const Eigen::Map<const Eigen::RowVector3d>      c4m2    (coeff.tail(3).data());



    auto xyzProd = ((points.topRows(1).array()) * (points.bottomRows(1).array()) * (points.template middleRows<1>(1)).array()).matrix();

    auto xyPoint  = (points.topRows(1).array() * points.topRows(1).array()
                         * points.template middleRows<1>(1).array()*points.template middleRows<1>(1).array()).matrix();

    auto yzPoint  = (points.template middleRows<1>(1).array() * points.template middleRows<1>(1).array()
                         * points.bottomRows(1).array() * points.bottomRows(1).array()).matrix();

    auto zxPoint  = (points.bottomRows(1).array() * points.bottomRows(1).array()
                         * points.topRows(1).array() * points.topRows(1).array()).matrix();


//     auto value1 = c0;
     auto value2 = c1 * points;
     auto value3 = (points.transpose().template leftCols<1>() * (c2_p1 * points)).diagonal().transpose();
     auto value4 = (points.transpose().template middleCols<1>(1)* (c2_p2 * points.template bottomRows<2>())).diagonal().transpose();
     auto value5 = (points.transpose().template rightCols<1>() * (c2_p3 * points.derived().template bottomRows<1>())).diagonal().transpose() ;
     auto value6 = ((points.array() * points.array()).matrix().transpose() * c3 * points).diagonal().transpose();
     auto value7 = c3m * xyzProd;
     auto value8 = ((points.array() *points.array() *points.array()).matrix().transpose() * c4 * points).diagonal().transpose();
     auto value9 = ((c4m1 * points).array() * xyzProd.array()).matrix();
     auto value10 = (c4m2(0)*xyPoint);
     auto value11 = (c4m2(1)*yzPoint);
     auto value12 = (c4m2(2)*zxPoint);
     return (c0 + (value2 + value3 + value4 + value5 + value6+ value7+ value8 + value9 + value10 + value11+ value12).array()).matrix();
   }
//
    
    template <typename Derived> auto leftRotated(const Eigen::MatrixBase<Derived>& row) { //row is a 3xD Eigen matrix
        assert(row.rows() == 3);
        Eigen::Matrix<double,3,1> newMatrix ({{row(1)},{row(2)},{row(0)}});
        return newMatrix;
    }
    
    template <typename Derived, typename DerivedOther> Polynomial<1> reductionOperator(const Eigen::MatrixBase<Derived>& A_, const Eigen::MatrixBase<DerivedOther>& a_){
        //A_,a_ are both 3x1 Eigen Column Vector
        assert(A_.rows() && a_.rows() == 3); assert(A_.cols() && a_.cols() == 1);
        
        const double                              c0      = coeff(0);
        Eigen::Map<const Eigen::RowVector3d>      c1      (coeff.data()+1);
        Eigen::Map<const Eigen::RowVector3d>      c2_p1   (coeff.data()+4);
        Eigen::Map<const Eigen::RowVector2d>      c2_p2   (coeff.data()+7);
        const  double                             c2_p3   = coeff(9);
        Eigen::Map<const Eigen::Matrix3d>         c3      (coeff.data()+10);
        const double                              c3m     = coeff(19);
        Eigen::Map<const Eigen::Matrix3d>         c4      (coeff.data()+20);
        Eigen::Map<const Eigen::RowVector3d>      c4m1    (coeff.data()+29);
        Eigen::Map<const Eigen::RowVector3d>      c4m2    (coeff.tail(3).data());
        
        
        auto                            A = A_(0);
        auto                            B = A_(1);
        auto                            C = A_(2);
        auto                            a = a_(0);
        auto                            b = a_(1);
        auto                            c = a_(2);
        auto                            ALR = leftRotated(A_);
        auto                            aLR = leftRotated(a_);
        auto                            Atwo = (A_.array()*A_.array()).matrix();
        auto                            atwo = (a_.array()*a_.array()).matrix();
        auto                            Athree = (A_.array()*A_.array()*A_.array()).matrix();
        auto                            athree = (a_.array()*a_.array()*a_.array()).matrix();
        auto                            aLRtwo = (aLR.array()*aLR.array()).matrix();
        auto                            ALRtwo = (ALR.array()*ALR.array()).matrix();


// Eigen expressions -> double
        auto w_0 = c0  + (c1*a_) + c3m*a*b*c +
                   ((a_.transpose().template leftCols<1>() * (c2_p1 * a_)).eval() +
                   a_.transpose().template middleCols<1>(1) * (c2_p2 * a_.template bottomRows<2>()) +
                   a_.transpose().template rightCols<1>() * (c2_p3 * a_.template bottomRows<1>())
                        + atwo.transpose()*c3*a_
                        + athree.transpose()*c4*a_
                        + c4m1*(a_*(a*b*c))
                        + c4m2*(atwo.array()*aLRtwo.array()).matrix()).value();

// Eigen expressions -> double
        auto w_1 = (c1*A_).value()+
                    ((c2_p1 * a_) * A_.transpose().template head<1>()).value() +
                   (A_.transpose().template segment<1>(1) * (c2_p2 * a_.template tail<2>())).value() +
                   (A_.transpose().template tail<1>() * (c2_p3 * a_.template tail<1>())).value() +
                    ((c2_p1 * A_) * a_.transpose().template head<1>() ).value() +
                   (a_.transpose().template segment<1>(1) * (c2_p2 * A_.template tail<2>())).value() +
                   (a_.transpose().template tail<1>() * (c2_p3 * A_.template tail<1>())).value()
                        + (aLRtwo.transpose()*c3*A_).value()
                        + ((2*A_.array()*a_.array()).matrix().transpose()*c3*a_).value()
                        + ((3*A_.array()*aLRtwo.array()).matrix().transpose()*c4*a_).value()
                        + ((3*athree).transpose()*c4*A_).value()
                        + (c4m1*(A_*(a*b*c))).value() + (c4m1*(a_*(a*b*C))).value() + (c4m1*(a_*(a*c*B))).value() + (c4m1*(A_*(b*c))).value()
                        + (c4m2*(a_.array()*aLRtwo.array()*A_.array()).matrix()).value()
                            + (c3m * (a*b*C + a*c*B + b*c*A));
// Eigen expressions -> double
        auto w_2 = ((c2_p1 * A_)*A_.transpose().template head<1>()).value() +
                    (A_.transpose().template segment<1>(1) * (c2_p2 * A_.template tail<2>())).value() +
                    (A_.transpose().template tail<1>() * (c2_p3 * A_.template tail<1>())).value()
                        + ((2*A_.array()*a_.array()).matrix().transpose()*c3*A_).value()
                        + (Atwo.transpose()*c3*a_).value()
                        + (3*(Atwo.array()*a_.array()).matrix().transpose()*c4*a_).value()
                        + (3*(A_.array()*atwo.array()).matrix().transpose()*c4*A_).value()
                        + (c4m1*(A_*(a*b*c))).value() + (c4m1*(A_*(a*c*B))).value() + (c4m1*(A_*(b*c*A))).value() + (c4m1*(a_*(a*B*C))).value()
                              + (c4m1*(a_*(b*A*C))).value() + (c4m1*(a_*(c*A*B))).value() + (c4m1*(a_*(A*B*C))).value()
                        + (c4m2*(aLRtwo.array()*ALRtwo.array()).matrix()).value()
                        + (c4m2*(atwo.array()*ALRtwo.array()).matrix()).value()
                        + c4m2*(4*a_.array()*aLR.array()*A_.array()*ALR.array()).matrix()
                            + (c3m*(a*B*C + b*A*C + c*A*B));

// Eigen expressions -> double
        auto w_3 = (Atwo.transpose()*c3*A_).value()
                        + c3m*A*B*C
                        + (Athree.transpose()*c4*a_).value()
                        + (3*(Atwo.array()*a_.array()).matrix().transpose()*c4*A_).value()
                        + (c4m1*(A_*(a*B*C)) + c4m1*(A_*(b*A*C)) + c4m1*(A_*(c*A*B))).value()
                        + (c4m2*(2*aLR.array()*Atwo.array()*ALR.array()).matrix()).value()
                        + (c4m2*(2*A_.array()*Atwo.array()).matrix()).value();

// Eigen expressions -> double
        auto w_4 = (Athree.transpose()*c4*A_).value()
                        + (c4m1*(A_*(A*B*C))).value()
                        + (c4m2*(Atwo.array()*ALRtwo.array()).matrix()).value();

        //returns a 5x1 Eigen Row vector which will be passed as coefficients for Poly<1> Class.
        return Polynomial<1>(Eigen::Matrix <double, 5, 1>{w_0 , w_1 , w_2 , w_3 , w_4});
    }

    

};
std::ostream &operator<<(std::ostream &t_os, Polynomial<3> const& t_poly){
    t_os << t_poly.coeff(0) << " " << PM(t_poly.coeff(1)) << t_poly.coeff(1) << "x " << PM(t_poly.coeff(2)) << t_poly.coeff(2) << "y "<< PM(t_poly.coeff(3)) << t_poly.coeff(3) << "z "<< PM(t_poly.coeff(4)) << t_poly.coeff(4) << "x^2 "<< PM(t_poly.coeff(5)) << t_poly.coeff(5) << "xy "<< PM(t_poly.coeff(6)) << t_poly.coeff(6) << "xz "<< PM(t_poly.coeff(7)) << t_poly.coeff(7) << "y^2 "<< PM(t_poly.coeff(8)) << t_poly.coeff(8) << "yz "<< PM(t_poly.coeff(9)) << t_poly.coeff(9) << "z^2 "<< PM(t_poly.coeff(10)) << t_poly.coeff(10) << "x^3 "<< PM(t_poly.coeff(11)) << t_poly.coeff(11) << "y^2x "<< PM(t_poly.coeff(12)) << t_poly.coeff(12) << "z^2x "<< PM(t_poly.coeff(13)) << t_poly.coeff(13) << "x^2y "<< PM(t_poly.coeff(14)) << t_poly.coeff(14) << "y^3 "<< PM(t_poly.coeff(15)) << t_poly.coeff(15) << "z^2y " << PM(t_poly.coeff(16)) << t_poly.coeff(16) << "x^2z "<< PM(t_poly.coeff(17)) << t_poly.coeff(17) << "y^2z "<< PM(t_poly.coeff(18)) << t_poly.coeff(18) << "z^3 "<< PM(t_poly.coeff(19)) << t_poly.coeff(19) << "xyz "<< PM(t_poly.coeff(20)) << t_poly.coeff(20) << "x^4 "<< PM(t_poly.coeff(21)) << t_poly.coeff(21) << "y^3x "<< PM(t_poly.coeff(22)) << t_poly.coeff(22) << "z^3x "<< PM(t_poly.coeff(23)) << t_poly.coeff(23) << "x^3y "<< PM(t_poly.coeff(24)) << t_poly.coeff(24) << "y^4 "<< PM(t_poly.coeff(25)) << t_poly.coeff(25) << "z^3y " << PM(t_poly.coeff(26)) << t_poly.coeff(26) << "x^3z "<< PM(t_poly.coeff(27)) << t_poly.coeff(27) << "y^z " << PM(t_poly.coeff(28)) << t_poly.coeff(28) << "z^4 "<< PM(t_poly.coeff(29)) << t_poly.coeff(29) << "x^2yz "<< PM(t_poly.coeff(30)) << t_poly.coeff(30) << "xy^2z "<< PM(t_poly.coeff(31)) << t_poly.coeff(31) << "xyz^2 "<< PM(t_poly.coeff(32)) << t_poly.coeff(32) << "x^2y^2 "<< PM(t_poly.coeff(33)) << t_poly.coeff(33) << "y^2z^2 "<< PM(t_poly.coeff(34)) << t_poly.coeff(34) << "z^2x^2 \n";
    
    return t_os;
}

#endif /* Polynomial_h */
