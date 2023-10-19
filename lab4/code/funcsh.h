#pragma once

#include<iostream>
#include <array>
#include <type_traits>
#include <tuple>

 

template <typename RealType, unsigned int N> struct nodes;

template <typename RealType> struct nodes<RealType, 1> {
    static constexpr std::array<RealType, 1> p{ 0.0 };
    static constexpr std::array<RealType, 1> w{ 2.0 };
};

template <typename RealType> struct nodes<RealType, 2> {
    static constexpr std::array<RealType, 2> p{ -0.5773502692, 0.5773502692 };
    static constexpr std::array<RealType, 2> w{ 1.0, 1.0 };
};

template <typename RealType> struct nodes<RealType, 3> {
    static constexpr std::array<RealType, 3> p{ -0.7745966692, 0, 0.7745966692 };
    static constexpr std::array<RealType, 3> w{ 0.5555555556, 0.8888888889, 0.5555555556 };
};

template <typename RealType> struct nodes<RealType, 4> {
    static constexpr std::array<RealType, 4> p{ -0.8611363116, -0.3399810436,
                                               0.3399810436, 0.8611363116 };
    static constexpr std::array<RealType, 4> w{ 0.3478548451,0.6521451549,
                                                0.6521451549, 0.3478548451 };
};

template <typename RealType> struct nodes<RealType, 5> {
    static constexpr std::array<RealType, 5> p{ -0.9061798459, -0.5384693101, 0.0,
                                               0.5384693101, 0.9061798459 };
    static constexpr std::array<RealType, 5> w{ 0.2369268851, 0.4786286705, 0.5688888889,
                                                0.4786286705, 0.2369268851 };
};


template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

 
template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate(
    const Callable& func,   
    const typename ArgumentGetter<Callable>::Argument& start,   
    const typename ArgumentGetter<Callable>::Argument& end,   
    const std::array<RealType, N>& points,   
    const std::array<RealType, N>& weights  
)
{
    RealType Int = 0;

    RealType semi_dif = (end - start) / 2;
    RealType semi_sum = (end + start) / 2;

    for (int i = 0; i < N; i++)
    {
        Int += weights[i] * func(semi_sum + semi_dif * points[i]);
    }

    return Int * semi_dif;
};

 
template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate(
    const Callable& func,   
    const typename ArgumentGetter<Callable>::Argument& start,  
    const typename ArgumentGetter<Callable>::Argument& end,   
    const std::array<RealType, N>& points,  
    const std::array<RealType, N>& weights,  
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx   
)
{
    unsigned int L = (end - start) / dx + 1;
    double step = (end - start) / L;

    RealType Int = 0;

    for (int i = 0; i < L; i++) {
        Int += integrate<Callable, RealType, N>(func, start + i * step,
            start + (i + 1) * step, points, weights);
    }

    return Int;
};

template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrateRR(
    const Callable& func,   
    const typename ArgumentGetter<Callable>::Argument& start,   
    const typename ArgumentGetter<Callable>::Argument& end,   
    const std::array<RealType, N>& points,   
    const std::array<RealType, N>& weights, 
    const RealType err  
)
{
    double h = end - start;

    RealType Ih = integrate<Callable, RealType, N>(func, start, end, points, weights, h);
    RealType Ih_2 = integrate<Callable, RealType, N>(func, start, end, points, weights, h / 2);
    RealType delta = h;
    RealType delta_2 = std::abs(Ih - Ih_2);

    for (double step = h / 4; delta > err; step = step / 2) {
        delta = delta_2;
        Ih = Ih_2;
        Ih_2 = integrate<Callable, RealType, N>(func, start, end, points, weights, step);
        delta_2 = std::abs(Ih - Ih_2);
    }

    double p = std::log(delta / delta_2) / std::log(2);

    return Ih_2 + (Ih - Ih_2) / (pow(2, p) - 1);
};
