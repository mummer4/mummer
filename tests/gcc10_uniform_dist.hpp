// Modified from GCC 10 uniform_int_dist.h Simplified and reduced functionality.
// Used as replacement for std::uniform_int_distribution using the same
// algorithm as GCC 10 for reproducibility of testing results. Original license
// below.

// Class template uniform_int_distribution -*- C++ -*-

// Copyright (C) 2009-2020 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/**
 * @file bits/uniform_int_dist.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{random}
 */

#ifndef _GCC10_UNIFORM_DIST_H
#define _GCC10_UNIFORM_DIST_H

#ifndef __glibcxx_assert
#define __glibcxx_assert(x)
#endif

#include <type_traits>
#include <limits>

/**
 * @brief Uniform discrete distribution for random numbers.
 * A discrete random distribution on the range @f$[min, max]@f$ with equal
 * probability throughout the range.
 */
template<typename _IntType = int>
class uniform_int_distribution
{
	static_assert(std::is_integral<_IntType>::value,
								"template argument must be an integral type");

public:
	/** The type of the range of the distribution. */
	typedef _IntType result_type;
	/** Parameter type. */
	struct param_type
	{
		typedef uniform_int_distribution<_IntType> distribution_type;

		param_type() : param_type(0) { }

		explicit
		param_type(_IntType __a,
							 _IntType __b = std::numeric_limits<_IntType>::max())
			: _M_a(__a), _M_b(__b)
			{
				__glibcxx_assert(_M_a <= _M_b);
			}

		result_type
		a() const
			{ return _M_a; }

		result_type
		b() const
			{ return _M_b; }

		friend bool
		operator==(const param_type& __p1, const param_type& __p2)
			{ return __p1._M_a == __p2._M_a && __p1._M_b == __p2._M_b; }

		friend bool
		operator!=(const param_type& __p1, const param_type& __p2)
			{ return !(__p1 == __p2); }

private:
		_IntType _M_a;
		_IntType _M_b;
	};

public:
	/**
	 * @brief Constructs a uniform distribution object.
	 */
	uniform_int_distribution() : uniform_int_distribution(0) { }

	/**
	 * @brief Constructs a uniform distribution object.
	 */
	explicit
	uniform_int_distribution(_IntType __a,
													 _IntType __b = std::numeric_limits<_IntType>::max())
		: _M_param(__a, __b)
		{ }

	explicit
	uniform_int_distribution(const param_type& __p)
		: _M_param(__p)
		{ }


	inline result_type a() const { return _M_param.a(); }
	inline result_type b() const { return _M_param.b(); }

	/**
	 * @brief Returns the inclusive lower bound of the distribution range.
	 */
	inline result_type min() const { return this->a(); }

	/**
	 * @brief Returns the inclusive upper bound of the distribution range.
	 */
	result_type max() const { return this->b(); }

	/**
	 * @brief Generating functions.
	 */
	template<typename _UniformRandomNumberGenerator>
	result_type
	operator()(_UniformRandomNumberGenerator& __urng)
		{ return this->operator()(__urng, _M_param); }

	template<typename _UniformRandomNumberGenerator>
	result_type
	operator()(_UniformRandomNumberGenerator& __urng,
						 const param_type& __p);

	/**
	 * @brief Return true if two uniform integer distributions have
	 *        the same parameters.
	 */
	friend bool
	operator==(const uniform_int_distribution& __d1,
						 const uniform_int_distribution& __d2)
		{ return __d1._M_param == __d2._M_param; }

private:
	param_type _M_param;
};

template<typename _IntType>
template<typename _UniformRandomNumberGenerator>
typename uniform_int_distribution<_IntType>::result_type
uniform_int_distribution<_IntType>::
operator()(_UniformRandomNumberGenerator& __urng,
					 const param_type& __param)
{
	typedef typename _UniformRandomNumberGenerator::result_type
	  _Gresult_type;
	typedef typename std::make_unsigned<result_type>::type __utype;
	typedef typename std::common_type<_Gresult_type, __utype>::type
	  __uctype;

	const __uctype __urngmin = __urng.min();
	const __uctype __urngmax = __urng.max();
	const __uctype __urngrange = __urngmax - __urngmin;
	const __uctype __urange
	  = __uctype(__param.b()) - __uctype(__param.a());

	__uctype __ret;

	if (__urngrange > __urange)
	{
		// downscaling
		const __uctype __uerange = __urange + 1; // __urange can be zero
		const __uctype __scaling = __urngrange / __uerange;
		const __uctype __past = __uerange * __scaling;
		do
			__ret = __uctype(__urng()) - __urngmin;
		while (__ret >= __past);
		__ret /= __scaling;
	}
	else if (__urngrange < __urange)
	{
		// upscaling
		/*
			Note that every value in [0, urange]
			can be written uniquely as

			(urngrange + 1) * high + low

			where

			high in [0, urange / (urngrange + 1)]

			and

			low in [0, urngrange].
		*/
		__uctype __tmp; // wraparound control
		do
		{
			const __uctype __uerngrange = __urngrange + 1;
			__tmp = (__uerngrange * operator()
							 (__urng, param_type(0, __urange / __uerngrange)));
			__ret = __tmp + (__uctype(__urng()) - __urngmin);
		}
		while (__ret > __urange || __ret < __tmp);
	}
	else
	  __ret = __uctype(__urng()) - __urngmin;

	return __ret + __param.a();
}

#endif
