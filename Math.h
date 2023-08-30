//
//  Math.h
//  Polynomial_Methods
//
//  Created by Keynesh Dongol on 22/7/23.
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

namespace math
{
    template<class E>
    class expression
    {
    public:
        auto size() const {
            return static_cast<E const&>(*this).size();
        }

        auto operator[](std::size_t i) const
        {
            if (i >= size())
                throw std::length_error("");
            return static_cast<E const&>(*this)[i];
        }

        operator E&() { return static_cast<E&>(*this); }
        operator E const&() const { return static_cast<E const&>(*this); }
    }; // class expression

    template<typename T, class Allocator = std::allocator<T>>
    class vector
            : public expression<vector<T>>
    {
    private:
        using data_type = std::vector<T, Allocator>;
        data_type m_data;

    public:
        using value_type = T;
        using allocator_type = Allocator;
        using size_type = typename data_type::size_type;
        using difference_type = typename data_type::difference_type;
        using reference = typename data_type::reference;
        using const_reference = typename data_type::const_reference;
        using pointer = typename data_type::pointer ;
        using const_pointer = typename data_type::const_pointer;

        vector(size_type d)
                : m_data(d)
        { }
        vector(std::initializer_list<value_type> init)
                : m_data(init)
        { }
        template<class E>
        vector(expression<E> const& expression)
                : m_data(expression.size())
        {
            for (size_type i = 0; i < expression.size(); ++i)
                m_data[i] = expression[i];
        }

        size_type size() const {
            return m_data.size();
        }

        value_type  operator[](size_type i) const { return m_data[i]; }
        value_type& operator[](size_type i)       { return m_data[i]; };
    }; // class vector

    namespace detail
    {
        template<typename T>
        class scalar
                : public expression<scalar<T>>
        {
        public:
            using value_type = T;
            using allocator_type = std::allocator<void>;
            using size_type = typename std::allocator<T>::size_type;
            using difference_type = typename std::allocator<T>::difference_type;
            using reference = typename std::allocator<T>::reference;
            using const_reference = typename std::allocator<T>::const_reference;
            using pointer = typename std::allocator<T>::pointer;
            using const_pointer = typename std::allocator<T>::const_pointer;

            scalar(value_type value)
                    : m_value(value)
            { }

            size_type size() const {
                return 0;
            }

            operator value_type&() { return static_cast<value_type&>(*this); }
            operator value_type const&() const { return static_cast<value_type const&>(*this); }

            value_type  operator[](size_type i) const { return m_value; }
            value_type& operator[](size_type i)       { return m_value; }

        private:
            value_type m_value;
        }; // class scalar

        template<class>
        struct is_scalar : std::false_type { };

        template<class T>
        struct is_scalar<scalar<T>> : std::true_type { };
    } // namespace detail

    template<class E1, class E2, class BinaryOperation>
    class vector_binary_operation
            : public expression<vector_binary_operation<E1, E2, BinaryOperation>>
    {
    public:
        using value_type = decltype(BinaryOperation()(typename E1::value_type(), typename E2::value_type()));
        using allocator_type = std::conditional_t<
                detail::is_scalar<E1>::value,
                typename E2::allocator_type::template rebind<value_type>::other,
                typename E1::allocator_type::template rebind<value_type>::other>;

    private:
        using vector_type = vector<value_type, allocator_type>;

    public:
        using size_type = typename vector_type::size_type;
        using difference_type = typename vector_type::difference_type;
        using reference = typename vector_type::reference;
        using const_reference = typename vector_type::const_reference;
        using pointer = typename vector_type::pointer;
        using const_pointer = typename vector_type::const_pointer;

        vector_binary_operation(expression<E1> const& e1, expression<E2> const& e2, BinaryOperation op)
                : m_e1(e1), m_e2(e2),
                  m_op(op)
        {
            if (e1.size() > 0 && e2.size() > 0 && !(e1.size() == e2.size()))
                throw std::logic_error("");
        }

        size_type size() const {
            return m_e1.size(); // == m_e2.size()
        }

        value_type operator[](size_type i) const {
            return m_op(m_e1[i], m_e2[i]);
        }

    private:
        E1 m_e1;
        E2 m_e2;
        //E1 const& m_e1;
        //E2 const& m_e2;
        BinaryOperation m_op;
    }; // class vector_binary_operation

    template<class E1, class E2>
    vector_binary_operation<E1, E2, std::plus<>>
    operator+(expression<E1> const& e1, expression<E2> const& e2) {
        return{ e1, e2, std::plus<>() };
    }
    template<class E1, class E2>
    vector_binary_operation<E1, E2, std::minus<>>
    operator-(expression<E1> const& e1, expression<E2> const& e2) {
        return{ e1, e2, std::minus<>() };
    }
    template<class E1, class E2>
    vector_binary_operation<E1, E2, std::multiplies<>>
    operator*(expression<E1> const& e1, expression<E2> const& e2) {
        return{ e1, e2, std::multiplies<>() };
    }
    template<class E1, class E2>
    vector_binary_operation<E1, E2, std::divides<>>
    operator/(expression<E1> const& e1, expression<E2> const& e2) {
        return{ e1, e2, std::divides<>() };
    }

    template<class E, typename T>
    vector_binary_operation<E, detail::scalar<T>, std::divides<>>
    operator/(expression<E> const& expr, T val) {
        return{ expr, detail::scalar<T>(val), std::divides<>() };
    }
    template<class E, typename T>
    vector_binary_operation<E, detail::scalar<T>, std::multiplies<>>
    operator*(T val, expression<E> const& expr) {
        return{ expr, detail::scalar<T>(val), std::multiplies<>() };
    }
    template<class E, typename T>
    vector_binary_operation<E, detail::scalar<T>, std::multiplies<>>
    operator*(expression<E> const& expr, T val) {
        return{ expr, detail::scalar<T>(val), std::multiplies<>() };
    }
} // namespace math



