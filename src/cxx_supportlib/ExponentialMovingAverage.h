/*
 *  Phusion Passenger - https://www.phusionpassenger.com/
 *  Copyright (c) 2016 Phusion Holding B.V.
 *
 *  "Passenger", "Phusion Passenger" and "Union Station" are registered
 *  trademarks of Phusion Holding B.V.
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */
#ifndef _PASSENGER_EXPONENTIAL_MOVING_AVERAGE_H_
#define _PASSENGER_EXPONENTIAL_MOVING_AVERAGE_H_

#include <oxt/macros.hpp>
#include <boost/config.hpp>
#include <algorithm>
#include <cmath>

namespace Passenger {

using namespace std;


/**
 * Implements discontiguous exponential averaging, as described by John C. Gunther 1998.
 * Can be used to compute moving exponentially decaying averages and standard deviations.
 * Unlike normal exponential moving average, this algorithm also works when the data has
 * gaps, and it also avoids initial value bias and postgap bias. See
 * http://www.drdobbs.com/tools/discontiguous-exponential-averaging/184410671
 *
 * ## Template parameters
 *
 * ### alpha
 *
 * This specifies by what factor data should decay. Its range is [0, 1000]. Lower values
 * cause data to decay more quickly, higher values cause data to decay more slowly.
 *
 * ### maxAge
 *
 * This represents an educational guess as to how long (in microseconds) it takes for
 * the sampled data sequence to change significantly. If you don't expect large random
 * variations then you should set this to a large value. For a data sequence dominated by
 * large random variations, setting this to 1000000 (1 second) might be appropriate.
 *
 * If the time interval between updates is `dt`, using a `maxAge` of `N * dt` will cause
 * each update to fill in up to `N - 1` of any preceeding skipped updates with the current
 * data value.
 *
 * ### alphaTimeUnit
 *
 * This specifies the time, in microseconds, after which the data should decay
 * by a factor of exactly `alpha`. For example, if `alpha = 0.5` and `alphaTimeUnit = 2000000`,
 * then data decays by 0.5 per 2 seconds.
 *
 * The default value is 1 second.
 */
template<
	unsigned int alpha,
	unsigned long long maxAge,
	unsigned long long alphaTimeUnit = 1000000
>
class DiscExponentialAverage {
private:
	double sumOfWeights, sumOfData, sumOfSquaredData;
	unsigned long long prevTime;

	static BOOST_CONSTEXPR double floatingAlpha() {
		return alpha / 1000.0;
	}

	static BOOST_CONSTEXPR double newDataWeightUpperBound() {
		return 1 - pow(floatingAlpha(), maxAge / (double) alphaTimeUnit);
	}

public:
	DiscExponentialAverage()
		: sumOfWeights(0),
		  sumOfData(0),
		  sumOfSquaredData(0),
		  prevTime(0)
		{ }

	void update(double value, unsigned long long now) {
		if (OXT_UNLIKELY(now <= prevTime)) {
			return;
		}

		double weightReductionFactor = pow(floatingAlpha(),
			(now - prevTime) / (double) alphaTimeUnit);
		double newDataWeight = std::min(1 - weightReductionFactor,
			newDataWeightUpperBound());
		sumOfWeights = weightReductionFactor * sumOfWeights + newDataWeight;
	  	sumOfData = weightReductionFactor * sumOfData + newDataWeight * value;
  		sumOfSquaredData = weightReductionFactor * sumOfSquaredData
  			+ newDataWeight * pow(value, 2.0);
  		prevTime = now;
	}

	bool available() const {
		return sumOfWeights > 0;
	}

	double completeness(unsigned long long now) const {
		return pow(floatingAlpha(), now - prevTime) * sumOfWeights;
	}

	double average() const {
		return sumOfData / sumOfWeights;
	}

	double stddev() const {
		return sqrt(sumOfSquaredData / sumOfWeights - pow(average(), 2));
	}
};


/**
 * Calculates a (normal) exponential moving average. This algorithm is not timing sensitive:
 * it doesn't take into account gaps in the data over time, and treats all values
 * equally regardless of when the value was collected. See also DiscExponentialAverage.
 *
 * You should initialize the the average value with a value equal to `nullValue`.
 * If `prevAverage` equals `nullValue` then this function simply returns `currentValue`.
 */
inline double
exponentialMovingAverage(double prevAverage, double currentValue, double alpha, double nullValue = -1) {
	if (OXT_UNLIKELY(prevAverage == nullValue)) {
		return currentValue;
	} else {
		return alpha * currentValue + (1 - alpha) * prevAverage;
	}
}


} // namespace Passenger

#endif /* _PASSENGER_EXPONENTIAL_MOVING_AVERAGE_H_ */
