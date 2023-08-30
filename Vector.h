//
//  Vector.h
//  PolyNested
//
//  Created by Keynesh Dongol on 22/8/23.
//

#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <numeric>
#include <algorithm>
#include <functional>
// For random number generator only (would normally be a separate <Random.h>)
#include <random>
#include <vector>
#include <cassert>
#include <utility>

template<typename ScalarType>
class Vector
{
public:
    explicit Vector() { std::cout << "ctor called.\n"; };
    explicit Vector(int size): _vec(size) { std::cout << "ctor called.\n"; };
    explicit Vector(const std::vector<ScalarType> &vec): _vec(vec)
    { std::cout << "ctor called.\n"; };

    Vector(const Vector<ScalarType> & vec): _vec{vec()}
    { std::cout << "copy ctor called.\n"; };
    Vector(Vector<ScalarType> && vec): _vec(std::move(vec()))
    { std::cout << "move ctor called.\n"; };

    Vector<ScalarType> & operator=(const Vector<ScalarType> &) = default;
    Vector<ScalarType> & operator=(Vector<ScalarType> &&) = default;

    decltype(auto) operator[](int indx) { return _vec[indx]; }
    decltype(auto) operator[](int indx) const { return _vec[indx]; }

    decltype(auto) operator()() & { return (_vec); };
    decltype(auto) operator()() const & { return (_vec); };
    Vector<ScalarType> && operator()() && { return std::move(*this); }

    int size() const { return _vec.size(); }

private:
    std::vector<ScalarType> _vec;
};

///
/// These are conventional overloads of operator + (the vector plus operation)
/// and operator * (the vector inner product operation) without using the expression
/// templates. They are later used for bench-marking purpose.
///

// + (vector plus) operator.
template<typename ScalarType>
auto operator+(const Vector<ScalarType> &lhs, const Vector<ScalarType> &rhs)
{
    assert(lhs().size() == rhs().size() &&
           "error: ops plus -> lhs and rhs size mismatch.");

    std::vector<ScalarType> _vec;
    _vec.resize(lhs().size());
    std::transform(std::cbegin(lhs()), std::cend(lhs()),
                   std::cbegin(rhs()), std::begin(_vec),
                   std::plus<>());
    return Vector<ScalarType>(std::move(_vec));
}

// * (vector inner product) operator.
template<typename ScalarType>
auto operator*(const Vector<ScalarType> &lhs, const Vector<ScalarType> &rhs)
{
    assert(lhs().size() == rhs().size() &&
           "error: ops multiplies -> lhs and rhs size mismatch.");

    std::vector<ScalarType> _vec;
    _vec.resize(lhs().size());
    std::transform(std::cbegin(lhs()), std::cend(lhs()),
                   std::cbegin(rhs()), std::begin(_vec),
                   std::multiplies<>());
    return Vector<ScalarType>(std::move(_vec));
}








/// Fwd declaration.
template<typename> class Vector;

namespace expr
{


/// -----------------------------------------
///
/// The first section is a base class template for all kinds of expression. It
/// employs the Curiously Recurring Template Pattern, which enables its instantiation
/// to any kind of expression structure inheriting from it.
///
/// -----------------------------------------


    /// Base class for all expressions.
    template<typename Expr> class expr_base
    {
    public:
        const Expr& self() const { return static_cast<const Expr&>(*this); }
        Expr& self() { return static_cast<Expr&>(*this); }

    protected:
        explicit expr_base() {};
        int size() const { return self().size_impl(); }
        auto operator[](int indx) const { return self().at_impl(indx); }
        auto operator()() const { return self()(); };
    };
/// -----------------------------------------
///
/// Now we are going to implement implementation of an abstraction onf pure algebraic
/// expressions. The idea being thast any PAE can be converted to a rael object
/// instance using operator(). It is here that the real computation is being carreid out.
///
/// -----------------------------------------

    /// Generic wrapper fro underlying data structure
    template<typename DataType> class terminal: expr_base<terminal<DataType>>
    {
    public:
        using base_type = expr_base<terminal<DataType>>;
        using base_type::size;
        using base_type::operator[];
        friend base_type;

        explicit terminal(const DataType &val): _val(val) {}
        int size_impl() const { return _val.size(); };
        auto at_impl(int indx) const { return _val[indx]; };
        decltype(auto) operator()() const { return (_val); }

    private:
        const DataType &_val;
    };

/// -----------------------------------------
///
/// Binary Operation Expression
///
/// This is a PAE abstraction of any binary exporession. Simiarly,  it inherits from expr_base.
///
/// It provides the size() method, indexed acess through at_impl() _ and a conversion
/// to referenced objects through the () operator. Each call to the at_impl() _method is
/// a element wise computation.
///
/// -----------------------------------------

    ///Genereic wrapper for binary operations (that are element-wise).
    template<typename Ops, typename lExpr, typename rExpr>
    class binary_ops: public expr_base<binary_ops<Ops,lExpr,rExpr>>
    {
    public:
        using base_type = expr_base<binary_ops<Ops,lExpr,rExpr>>;
        using base_type::size;
        using base_type::operator[];
        friend base_type;

        explicit binary_ops(const Ops &ops, const lExpr &lxpr, const rExpr &rxpr)
                : _ops(ops), _lxpr(lxpr), _rxpr(rxpr) {};
        int size_impl() const { return _lxpr.size(); };

        /// This does the element-wise computation for index indx.
        auto at_impl(int indx) const { return _ops(_lxpr[indx], _rxpr[indx]); };

        /// Conversion from arbitrary expr to concrete data type. It evaluates
        /// element-wise computations for all indices.
        template<typename DataType> operator DataType()
        {
            DataType _vec(size());
            for(int _ind = 0; _ind < _vec.size(); ++_ind)
                _vec[_ind] = (*this)[_ind];
            return _vec;
        }

    private: /// Ops and expr are assumed cheap to copy.
        Ops   _ops;
        lExpr _lxpr;
        rExpr _rxpr;
    };



/// -----------------------------------------
/// Now, we will define rwo structure that defines the algrebraic operations on the PAE.
/// Methods implemented are:
///     - vector plus
///     - vector inner product
///
/// -----------------------------------------

    ///Element wise operation
    struct vec_plus_t
    {
        constexpr explicit vec_plus_t() = default;
        template<typename LType, typename RType>
        auto operator()(const LType &lhs, const RType &rhs) const
        { return lhs+rhs; }
    };

    /// Element-wise inner product operation.
    struct vec_prod_t
    {
        constexpr explicit vec_prod_t() = default;
        template<typename LType, typename RType>
        auto operator()(const LType &lhs, const RType &rhs) const
        { return lhs*rhs; }
    };

    /// Constant plus and inner product operator objects.
    constexpr vec_plus_t vec_plus{};
    constexpr vec_prod_t vec_prod{};

    /// Plus operator overload on expressions: return binary expression.
    template<typename lExpr, typename rExpr>
    auto operator+(const lExpr &lhs, const rExpr &rhs)
    { return binary_ops<vec_plus_t,lExpr,rExpr>(vec_plus,lhs,rhs); }

    /// Inner prod operator overload on expressions: return binary expression.
    template<typename lExpr, typename rExpr>
    auto operator*(const lExpr &lhs, const rExpr &rhs)
    { return binary_ops<vec_prod_t,lExpr,rExpr>(vec_prod,lhs,rhs); }




} //!expr












//namespace Container {
//
//    /*!
//     *Our goal here is to simply created a templated class called containter, that is builds off
//     *std::array's fixed-size array of memory allocation for storing the "component" of some
//     *structure one plans on passing. Here we will implememt methods to enable arthmetic
//     *operations.
//     *
//     *We can think of this as a 1-Dimensional Vector of N elements templated class of scalar types.
//     *Thus, we can expect the default call to be:
//     *
//     *-> Vector<double, 5> someVec;
//     *
//     *and we can construct our vector using curly brackets, as such:
//     *
//     *-> Vector<double, 5 > someVec = {5,4,3,2,1};
//     */
//
//template <typename Scalar, size_t Ndim>
//struct Vector;
//
////    This is template friendly manner to overload the stream operator. https://stackoverflow.com/questions/4660123
//template <typename Scalar, size_t Ndim>
//std::ostream& operator<< (std::ostream&, const Vector<Scalar, Ndim>&);
//
////Declaration of random number genration for testing of code later
//template <typename T> class RandUniformReal;
//template <typename T> class RandUniformInt;

//template <typename ScalarType>
//class Vector {
//public:
//    explicit Vector() { std::cout << "ctor called/ \n";}
//    explicit Vector(int size): _vec(size) { std::cout << "ctor call. \n"; }
//    explicit Vector(const std::vector<ScalarType>& vec): _vec(vec){ std::cout << "ctor called. \n"; }
//
//    Vector(const Vector<ScalarType>& vec): _vec(vec()) { std::cout << "copy ctor called. \n"; }
//    Vector(const Vector<ScalarType>&& vec): _vec(std::move(vec())) { std::cout << "move ctor called. \n"; }
//
//    Vector<ScalarType>& operator= (const Vector<ScalarType> &) = default;
//    Vector<ScalarType>& operator= (Vector<ScalarType> &&) = default;
//
//    decltype(auto) operator[](int indx) { return _vec[indx]; }
//    decltype(auto) operator[](int indx) const { return _vec[indx]; }
//
//    decltype(auto) operator()() & { return (_vec); }
//    decltype(auto) operator()() const & { return (_vec); }
//    Vector<ScalarType> && operator()() && { return std::move(*this); }
//
//    int size() const { return _vec.size(); }
//
//
//private:
//    std::vector<ScalarType> _vec;
//
//};
//
//// + (vector plus) operator.
//template <typename ScalarType>
//auto operator+ (const Vector<ScalarType>& lhs, const Vector<ScalarType>& rhs) {
//    assert(lhs.size() == rhs.size() &&
//           "error: binaryOps plus(+) -> lhs and rhs size are not the same.");
//
//    std::vector<ScalarType> _vec;
//    _vec.resize(lhs.size());
//    std::transform(std::cbegin(lhs()), std::cend(lhs()),
//                   std::cbegin(rhs()), std::begin(_vec),
//                   std::plus<>());
//    return Vector<ScalarType>(std::move(_vec));
//}
//
//// * (vector inner product) operator
//template <typename ScalarType>
//auto operator* (const Vector<ScalarType>& lhs, const Vector<ScalarType>& rhs) {
//
//}

//} // namespace Container


