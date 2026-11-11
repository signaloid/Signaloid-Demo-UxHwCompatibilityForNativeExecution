/*
 *	Copyright (c) 2023–2026, Signaloid.
 *
 *	Permission is hereby granted, free of charge, to any person obtaining a copy
 *	of this software and associated documentation files (the "Software"), to deal
 *	in the Software without restriction, including without limitation the rights
 *	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the Software is
 *	furnished to do so, subject to the following conditions:
 *
 *	The above copyright notice and this permission notice shall be included in all
 *	copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *	SOFTWARE.
 */

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include "uxhw.h"

enum
{
	kMicrosecondsPerSecond = 1000000,
};

static gsl_rng * gGSLr;

/**
 *	@brief	Initializes the libGSL random number generator. The generator type and seed
 *		can be set with the GSL_RNG_TYPE and GSL_RNG_SEED environment variables.
 *		If GSL_RNG_SEED is not set, the seed comes from the current time.
 *
 */
static void
initializeGenerators(void)
{
	/*
	 *	gsl_rng_env_setup() reads GSL_RNG_TYPE and GSL_RNG_SEED into
	 *	gsl_rng_default and gsl_rng_default_seed, so gsl_rng_default must be
	 *	read after it. gsl_rng_alloc() seeds with gsl_rng_default_seed.
	 */
	gsl_rng_env_setup();
	gGSLr = gsl_rng_alloc(gsl_rng_default);

	if (getenv("GSL_RNG_SEED") == NULL)
	{
		struct timeval t = { .tv_sec = 0, .tv_usec = 0 };

		gettimeofday(&t, NULL);
		gsl_rng_set(gGSLr, (unsigned long) t.tv_sec * kMicrosecondsPerSecond + (unsigned long) t.tv_usec);
	}

	return;
}

/**
 *	@brief	Returns a random index in the interval [0,`n`).
 *
 *	@param	n	: The upper bound of the interval of the possible indices to be returned.
 *	@return	ssize_t	: Returns an index `i`, where 0 <= `i` < `n`. Returns -1 if `n` is 0
 *			  or exceeds the range of the random number generator.
 */
static ssize_t
randomFromRange(size_t n)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	const unsigned long generatorRange = gsl_rng_max(gGSLr) - gsl_rng_min(gGSLr);

	if ((n == 0) || (n > generatorRange))
	{
		return -1;
	}

	return (ssize_t) gsl_rng_uniform_int(gGSLr, n);
}

double
UxHwDoubleSample(double value)
{
	return value;
}

float
UxHwFloatSample(float value)
{
	return value;
}

void
UxHwDoubleSampleBatch(double value, double *  destSampleArray, size_t numberOfRandomSamples)
{
	if (destSampleArray == NULL)
	{
		fprintf(stderr, "UxHwDoubleSampleBatch: destSampleArray is NULL.\n");

		return;
	}

	for (size_t ii = 0; ii < numberOfRandomSamples; ii++)
	{
		destSampleArray[ii] = value;
	}

	return;
}

void
UxHwFloatSampleBatch(float value, float *  destSampleArray, size_t numberOfRandomSamples)
{
	if (destSampleArray == NULL)
	{
		fprintf(stderr, "UxHwFloatSampleBatch: destSampleArray is NULL.\n");

		return;
	}

	for (size_t ii = 0; ii < numberOfRandomSamples; ii++)
	{
		destSampleArray[ii] = value;
	}

	return;
}

double
UxHwDoubleDistFromSamples(double * samples, size_t sampleCount)
{
	if (samples == NULL)
	{
		fprintf(stderr, "UxHwDoubleDistFromSamples: samples is NULL. Returning NAN...\n");

		return NAN;
	}

	const ssize_t randomIndex = randomFromRange(sampleCount);

	if (randomIndex < 0)
	{
		fprintf(stderr, "UxHwDoubleDistFromSamples: sampleCount is 0 or too large. Returning NAN...\n");

		return NAN;
	}

	return samples[randomIndex];
}

float
UxHwFloatDistFromSamples(float * samples, size_t sampleCount)
{
	if (samples == NULL)
	{
		fprintf(stderr, "UxHwFloatDistFromSamples: samples is NULL. Returning NAN...\n");

		return NAN;
	}

	const ssize_t randomIndex = randomFromRange(sampleCount);

	if (randomIndex < 0)
	{
		fprintf(stderr, "UxHwFloatDistFromSamples: sampleCount is 0 or too large. Returning NAN...\n");

		return NAN;
	}

	return samples[randomIndex];
}

double
UxHwDoubleDistFromWeightedSamples(WeightedDoubleSample * samples, size_t weightedSampleCount)
{
	double *    normalizedCumulativeWeights;
	double      weightSum = 0.0;
	double      probability;
	size_t      index = 0;

	if (samples == NULL)
	{
		fprintf(stderr, "UxHwDoubleDistFromWeightedSamples: samples is NULL. Returning NAN...\n");

		return NAN;
	}

	if (weightedSampleCount == 0)
	{
		fprintf(stderr, "UxHwDoubleDistFromWeightedSamples: weightedSampleCount is 0. Returning NAN...\n");

		return NAN;
	}

	normalizedCumulativeWeights = (double *) malloc(weightedSampleCount * sizeof(double));

	if (normalizedCumulativeWeights == NULL)
	{
		fprintf(stderr, "UxHwDoubleDistFromWeightedSamples: malloc failed. Returning NAN...\n");

		return NAN;
	}

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		if ((samples[ii].sampleWeight < 0.0) || !isfinite(samples[ii].sampleWeight))
		{
			fprintf(stderr, "UxHwDoubleDistFromWeightedSamples: weights must be finite and non-negative. Returning NAN...\n");
			free(normalizedCumulativeWeights);

			return NAN;
		}

		weightSum += samples[ii].sampleWeight;
		normalizedCumulativeWeights[ii] = weightSum;
	}

	/*
	 *	A positive, finite sum makes the last cumulative weight exactly 1.0
	 *	after normalization. probability is in [0, 1), so the search below
	 *	always sets index. The strict comparison never picks a zero-weight
	 *	sample.
	 */
	if ((weightSum <= 0.0) || !isfinite(weightSum))
	{
		fprintf(stderr, "UxHwDoubleDistFromWeightedSamples: sum of weights must be positive and finite. Returning NAN...\n");
		free(normalizedCumulativeWeights);

		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	probability = gsl_ran_flat(gGSLr, 0.0, 1.0);

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		normalizedCumulativeWeights[ii] /= weightSum;

		if (probability < normalizedCumulativeWeights[ii])
		{
			index = ii;
			break;
		}
	}

	free(normalizedCumulativeWeights);

	return samples[index].sample;
}

float
UxHwFloatDistFromWeightedSamples(WeightedFloatSample * samples, size_t weightedSampleCount)
{
	float * normalizedCumulativeWeights;
	float   weightSum = 0.0;
	double  probability;
	size_t  index = 0;

	if (samples == NULL)
	{
		fprintf(stderr, "UxHwFloatDistFromWeightedSamples: samples is NULL. Returning NAN...\n");

		return NAN;
	}

	if (weightedSampleCount == 0)
	{
		fprintf(stderr, "UxHwFloatDistFromWeightedSamples: weightedSampleCount is 0. Returning NAN...\n");

		return NAN;
	}

	normalizedCumulativeWeights = (float *) malloc(weightedSampleCount * sizeof(float));

	if (normalizedCumulativeWeights == NULL)
	{
		fprintf(stderr, "UxHwFloatDistFromWeightedSamples: malloc failed. Returning NAN...\n");

		return NAN;
	}

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		if ((samples[ii].sampleWeight < 0.0f) || !isfinite(samples[ii].sampleWeight))
		{
			fprintf(stderr, "UxHwFloatDistFromWeightedSamples: weights must be finite and non-negative. Returning NAN...\n");
			free(normalizedCumulativeWeights);

			return NAN;
		}

		weightSum += samples[ii].sampleWeight;
		normalizedCumulativeWeights[ii] = weightSum;
	}

	/*
	 *	A positive, finite sum makes the last cumulative weight exactly 1.0
	 *	after normalization. probability is in [0, 1), so the search below
	 *	always sets index. The strict comparison never picks a zero-weight
	 *	sample. probability is a double, because rounding it to float could
	 *	give exactly 1.0.
	 */
	if ((weightSum <= 0.0f) || !isfinite(weightSum))
	{
		fprintf(stderr, "UxHwFloatDistFromWeightedSamples: sum of weights must be positive and finite. Returning NAN...\n");
		free(normalizedCumulativeWeights);

		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	probability = gsl_ran_flat(gGSLr, 0.0, 1.0);

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		normalizedCumulativeWeights[ii] /= weightSum;

		if (probability < normalizedCumulativeWeights[ii])
		{
			index = ii;
			break;
		}
	}

	free(normalizedCumulativeWeights);

	return samples[index].sample;
}

void
UxHwDoubleDistFromMultidimensionalSamples(double * destinationArray, void *  samples, size_t sampleCount, size_t sampleCardinality)
{
	if ((destinationArray == NULL) || (samples == NULL))
	{
		fprintf(stderr, "UxHwDoubleDistFromMultidimensionalSamples: destinationArray or samples is NULL.\n");

		return;
	}

	double * *      castSamples = (double * *) samples;
	const ssize_t   randomIndex = randomFromRange(sampleCount);

	if (randomIndex < 0)
	{
		fprintf(stderr, "UxHwDoubleDistFromMultidimensionalSamples: sampleCount is 0 or too large.\n");

		return;
	}

	for (size_t ii = 0; ii < sampleCardinality; ii++)
	{
		destinationArray[ii] = castSamples[randomIndex][ii];
	}

	return;
}

void
UxHwFloatDistFromMultidimensionalSamples(float * destinationArray, void *  samples, size_t sampleCount, size_t sampleCardinality)
{
	if ((destinationArray == NULL) || (samples == NULL))
	{
		fprintf(stderr, "UxHwFloatDistFromMultidimensionalSamples: destinationArray or samples is NULL.\n");

		return;
	}

	float * *       castSamples = (float * *) samples;
	const ssize_t   randomIndex = randomFromRange(sampleCount);

	if (randomIndex < 0)
	{
		fprintf(stderr, "UxHwFloatDistFromMultidimensionalSamples: sampleCount is 0 or too large.\n");

		return;
	}

	for (size_t ii = 0; ii < sampleCardinality; ii++)
	{
		destinationArray[ii] = castSamples[randomIndex][ii];
	}

	return;
}

double
UxHwDoubleExponentialDist(double mu)
{
	if (mu < 0.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_exponential(gGSLr, mu);
}

float
UxHwFloatExponentialDist(float mu)
{
	return (float) UxHwDoubleExponentialDist((double) mu);
}

double
UxHwDoubleGumbel1Dist(double mu, double beta)
{
	if (beta < 0.0)
	{
		return NAN;
	}

	if (beta == 0.0)
	{
		return mu;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_gumbel1(gGSLr, 1.0 / beta, 1.0) + mu;
}

float
UxHwFloatGumbel1Dist(float mu, float beta)
{
	return (float) UxHwDoubleGumbel1Dist((double) mu, (double) beta);
}

double
UxHwDoubleUniformDist(double a, double b)
{
	if (a > b)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_flat(gGSLr, a, b);
}

float
UxHwFloatUniformDist(float a, float b)
{
	return (float) UxHwDoubleUniformDist((double) a, (double) b);
}

double
UxHwDoubleGaussDist(double mu, double sigma)
{
	double sample;

	if (sigma < 0.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	sample = gsl_ran_gaussian_ziggurat(gGSLr, sigma);

	return mu + sample;
}

float
UxHwFloatGaussDist(float mu, float sigma)
{
	return (float) UxHwDoubleGaussDist((double) mu, (double) sigma);
}

double
UxHwDoubleLogisticDist(double location, double scale)
{
	double sample;

	if (scale < 0.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	sample = gsl_ran_logistic(gGSLr, scale);

	return location + sample;
}

float
UxHwFloatLogisticDist(float location, float scale)
{
	return (float) UxHwDoubleLogisticDist((double) location, (double) scale);
}

double
UxHwDoubleLaplaceDist(double mu, double b)
{
	double sample;

	if (b < 0.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	sample = gsl_ran_laplace(gGSLr, b);

	return mu + sample;
}

float
UxHwFloatLaplaceDist(float mu, float b)
{
	return (float) UxHwDoubleLaplaceDist((double) mu, (double) b);
}

double
UxHwDoubleWeibullDist(double k, double lambda)
{
	if ((k <= 0.0) || (lambda < 0.0))
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_weibull(gGSLr, lambda, k);
}

float
UxHwFloatWeibullDist(float k, float lambda)
{
	return (float) UxHwDoubleWeibullDist((double) k, (double) lambda);
}

double
UxHwDoubleLognormalDist(double mu, double sigma)
{
	if (sigma < 0.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_lognormal(gGSLr, mu, sigma);
}

float
UxHwFloatLognormalDist(float mu, float sigma)
{
	return (float) UxHwDoubleLognormalDist((double) mu, (double) sigma);
}

double
UxHwDoubleBoundedparetoDist(double alpha, double xMin, double xMax)
{
	double sample = xMax;

	/*
	 *	Reject parameters outside the support of the distribution.
	 */
	if ((alpha <= 0.0) || (xMin <= 0.0) || (xMax <= xMin))
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	/*
	 *	Inverse-CDF sampling, so one draw is enough even when xMax is close to
	 *	xMin or alpha is small. supportMass is 1 - (xMin / xMax)^alpha. expm1()
	 *	and log1p() keep precision when it is close to 0. Rounding can still give
	 *	a value at or above xMax, so draw again in that case.
	 *
	 *	log(xMin / xMax) is the most accurate form when xMax is close to xMin.
	 *	If the ratio is subnormal or underflows to 0, use the difference of the
	 *	logs instead.
	 */
	const double    supportRatio    = xMin / xMax;
	const double    logSupportRatio = (supportRatio >= DBL_MIN) ? log(supportRatio) : (log(xMin) - log(xMax));
	const double    supportMass     = -expm1(alpha * logSupportRatio);

	while (sample >= xMax)
	{
		const double    uniformSample   = gsl_rng_uniform(gGSLr);
		const double    logGrowth       = -log1p(-uniformSample * supportMass) / alpha;
		const double    growth          = exp(logGrowth);

		/*
		 *	The sample is xMin * exp(logGrowth). exp() can overflow even when
		 *	the sample is finite, so in that case add the logs instead. The
		 *	product is kept otherwise, as it is the more accurate form near xMin.
		 */
		if (isinf(growth))
		{
			sample = exp(log(xMin) + logGrowth);
		}
		else
		{
			sample = xMin * growth;
		}
	}

	return sample;
}

float
UxHwFloatBoundedparetoDist(float alpha, float xMin, float xMax)
{
	return (float) UxHwDoubleBoundedparetoDist((double) alpha, (double) xMin, (double) xMax);
}

double
UxHwDoubleBetaDist(double a, double b)
{
	if ((a <= 0.0) || (b <= 0.0))
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_beta(gGSLr, a, b);
}

float
UxHwFloatBetaDist(float a, float b)
{
	return (float) UxHwDoubleBetaDist((double) a, (double) b);
}

double
UxHwDoubleGammaDist(double k, double theta)
{
	if ((k <= 0.0) || (theta < 0.0))
	{
		return NAN;
	}

	if (theta == 0.0)
	{
		return 0.0;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_gamma(gGSLr, k, theta);
}

float
UxHwFloatGammaDist(float k, float theta)
{
	return (float) UxHwDoubleGammaDist((double) k, (double) theta);
}

double
UxHwDoubleInverseGammaDist(double alpha, double beta)
{
	if ((alpha <= 1.0) || (beta < 0.0))
	{
		return NAN;
	}

	if (beta == 0.0)
	{
		return 0.0;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return 1.0 / gsl_ran_gamma(gGSLr, alpha, 1.0 / beta);
}

float
UxHwFloatInverseGammaDist(float alpha, float beta)
{
	return (float) UxHwDoubleInverseGammaDist((double) alpha, (double) beta);
}

double
UxHwDoubleChiSquaredDist(double k)
{
	if (k <= 0.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_chisq(gGSLr, k);
}

float
UxHwFloatChiSquaredDist(float k)
{
	return (float) UxHwDoubleChiSquaredDist((double) k);
}

double
UxHwDoubleFDist(double d1, double d2)
{
	if ((d1 <= 0.0) || (d2 <= 2.0))
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_fdist(gGSLr, d1, d2);
}

float
UxHwFloatFDist(float d1, float d2)
{
	return (float) UxHwDoubleFDist((double) d1, (double) d2);
}

double
UxHwDoubleStudentsTDist(double nu)
{
	if (nu <= 1.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_tdist(gGSLr, nu);
}

float
UxHwFloatStudentsTDist(float nu)
{
	return (float) UxHwDoubleStudentsTDist((double) nu);
}

double
UxHwDoubleGEVDist(double xi)
{
	double uniformSample;

	if (xi >= 1.0)
	{
		return NAN;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	/*
	 *	Sample the standard GEV distribution (zero location, unit scale) via
	 *	inverse-transform sampling. We reject the open-interval endpoints so
	 *	that the nested logarithms below remain finite.
	 */
	do
	{
		uniformSample = gsl_ran_flat(gGSLr, 0.0, 1.0);
	}
	while ((uniformSample <= 0.0) || (uniformSample >= 1.0));

	const double    logLog          = log(-log(uniformSample));
	const double    shapeEpsilon    = 1e-12;

	/*
	 *	Use expm1() rather than pow(..., -xi) - 1.0 to avoid catastrophic
	 *	cancellation for small |xi|. The result tends to -logLog as xi -> 0,
	 *	which is the standard Gumbel limit handled by the epsilon branch.
	 */
	if (fabs(xi) < shapeEpsilon)
	{
		return -logLog;
	}

	return expm1(-xi * logLog) / xi;
}

float
UxHwFloatGEVDist(float xi)
{
	return (float) UxHwDoubleGEVDist((double) xi);
}

double
UxHwDoubleMixture(double a, double b, double aScale)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	if (gsl_rng_uniform(gGSLr) < aScale)
	{
		return a;
	}

	return b;
}

float
UxHwFloatMixture(float a, float b, float aScale)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	if (gsl_rng_uniform(gGSLr) < aScale)
	{
		return a;
	}

	return b;
}

double
UxHwDoubleNthMoment(double value, size_t n)
{
	if (n == 0)
	{
		return 1.0;
	}

	if (n == 1)
	{
		return value;
	}

	return 0.0;
}

float
UxHwFloatNthMoment(float value, size_t n)
{
	if (n == 0)
	{
		return 1.0;
	}

	if (n == 1)
	{
		return value;
	}

	return 0.0;
}

double
UxHwDoubleNthMode(double value, size_t n)
{
	(void) n;

	return value;
}

float
UxHwFloatNthMode(float value, size_t n)
{
	(void) n;

	return value;
}

double
UxHwDoubleSupportMin(double value)
{
	return value;
}

float
UxHwFloatSupportMin(float value)
{
	return value;
}

double
UxHwDoubleSupportMax(double value)
{
	return value;
}

float
UxHwFloatSupportMax(float value)
{
	return value;
}

double
UxHwDoubleProbabilityGT(double value, double cutoff)
{
	return value > cutoff ? 1 : 0;
}

float
UxHwFloatProbabilityGT(float value, float cutoff)
{
	return value > cutoff ? 1 : 0;
}

double
UxHwDoubleGetIndependentCopy(double value)
{
	return value;
}

float
UxHwFloatGetIndependentCopy(float value)
{
	return value;
}

void
UxHwDoubleGetIndependentJointMultidimensionalCopy(double * srcDistArray, double *  destDistArray, size_t numberOfDistributions)
{
	if ((srcDistArray == NULL) || (destDistArray == NULL))
	{
		fprintf(stderr, "UxHwDoubleGetIndependentJointMultidimensionalCopy: srcDistArray or destDistArray is NULL.\n");

		return;
	}

	for (size_t ii = 0; ii < numberOfDistributions; ii++)
	{
		destDistArray[ii] = srcDistArray[ii];
	}

	return;
}

void
UxHwFloatGetIndependentJointMultidimensionalCopy(float * srcDistArray, float *  destDistArray, size_t numberOfDistributions)
{
	if ((srcDistArray == NULL) || (destDistArray == NULL))
	{
		fprintf(stderr, "UxHwFloatGetIndependentJointMultidimensionalCopy: srcDistArray or destDistArray is NULL.\n");

		return;
	}

	for (size_t ii = 0; ii < numberOfDistributions; ii++)
	{
		destDistArray[ii] = srcDistArray[ii];
	}

	return;
}

double
UxHwDoubleLimitDistributionSupport(double value, double supportMin, double supportMax)
{
	if ((value >= supportMin) && (value <= supportMax))
	{
		return value;
	}

	return NAN;
}

float
UxHwFloatLimitDistributionSupport(float value, float supportMin, float supportMax)
{
	if ((value >= supportMin) && (value <= supportMax))
	{
		return value;
	}

	return NAN;
}

double
UxHwDoubleQuantile(double value, double probability)
{
	(void) probability;

	return value;
}

float
UxHwFloatQuantile(float value, float probability)
{
	(void) probability;

	return value;
}

double
UxHwDoubleEvaluatePDF(double value, double domainValue)
{
	if (value == domainValue)
	{
		return INFINITY;
	}
	else
	{
		return 0.0;
	}
}

float
UxHwFloatEvaluatePDF(float value, float domainValue)
{
	if (value == domainValue)
	{
		return INFINITY;
	}
	else
	{
		return 0.0;
	}
}

double
UxHwDoubleBayesLaplace(double ( * statisticalModel )(void *, double), void * modelParams, double prior, double observedData, size_t numberOfObservations)
{
	(void) statisticalModel;
	(void) modelParams;
	(void) observedData;
	(void) numberOfObservations;

	fprintf(stderr, "Warning: UxHwDoubleBayesLaplace is not supported in native execution mode! Returning prior...\n");

	return prior;
}

float
UxHwFloatBayesLaplace(float ( * statisticalModel )(void *, float), void * modelParams, float prior, float observedData, size_t numberOfObservations)
{
	(void) statisticalModel;
	(void) modelParams;
	(void) observedData;
	(void) numberOfObservations;

	fprintf(stderr, "Warning: UxHwFloatBayesLaplace is not supported in native execution mode! Returning prior...\n");

	return prior;
}

/*
 *	Native execution carries no distributions to differentiate, so the
 *	derivative is unavailable and the stub returns NaN, as the platform
 *	does when a derivative is unavailable.
 */
double
UxHwDoubleDerivative(double value, double *  variables, const uint8_t *  orders, size_t numberOfVariables)
{
	(void) value;
	(void) variables;
	(void) orders;
	(void) numberOfVariables;

	return NAN;
}

float
UxHwFloatDerivative(float value, float *  variables, const uint8_t *  orders, size_t numberOfVariables)
{
	(void) value;
	(void) variables;
	(void) orders;
	(void) numberOfVariables;

	return NAN;
}

double
UxHwDoubleArgmin(
	double      functionOutput,
	double *    argumentArray,
	size_t      numberOfArguments,
	double *    minimizingArgumentInstanceArray)
{
	if ((minimizingArgumentInstanceArray != NULL) && (argumentArray != NULL))
	{
		for (size_t ii = 0; ii < numberOfArguments; ii++)
		{
			minimizingArgumentInstanceArray[ii] = argumentArray[ii];
		}
	}

	return functionOutput;
}

float
UxHwFloatArgmin(
	float   functionOutput,
	float * argumentArray,
	size_t  numberOfArguments,
	float * minimizingArgumentInstanceArray)
{
	if ((minimizingArgumentInstanceArray != NULL) && (argumentArray != NULL))
	{
		for (size_t ii = 0; ii < numberOfArguments; ii++)
		{
			minimizingArgumentInstanceArray[ii] = argumentArray[ii];
		}
	}

	return functionOutput;
}

double
UxHwDoubleArgminTotallyCorrelatedInputs(
	double      functionOutput,
	double *    argumentArray,
	size_t      numberOfArguments,
	double *    minimizingArgumentInstanceArray)
{
	return UxHwDoubleArgmin(
		functionOutput,
		argumentArray,
		numberOfArguments,
		minimizingArgumentInstanceArray
	);
}

float
UxHwFloatArgminTotallyCorrelatedInputs(
	float   functionOutput,
	float * argumentArray,
	size_t  numberOfArguments,
	float * minimizingArgumentInstanceArray)
{
	return UxHwFloatArgmin(
		functionOutput,
		argumentArray,
		numberOfArguments,
		minimizingArgumentInstanceArray
	);
}

void
UxHwFloatGeneratePath(
	void * parameterStructArray[],
	size_t numberOfParameterStructs,
	size_t pathLength,
	size_t numberOfStateVariables,
	void (* stateGeneratorFuncPtr)(void * parameterStruct, float * *  paths, size_t iterationStepCount, size_t numberOfStateVariables),
	float * * *   resultArray)
{
	/*
	 *	Sanity checking of input arguments
	 */
	if (pathLength == 0)
	{
		fprintf(stderr, "Error: UxHwFloatGeneratePath - Function must be called with a pathLength greater than 0.\n");

		return;
	}

	if (numberOfStateVariables == 0)
	{
		fprintf(stderr, "Error: UxHwFloatGeneratePath - Function must be called with numberOfStateVariables greater than 0.\n");

		return;
	}

	if (numberOfParameterStructs == 0)
	{
		fprintf(stderr, "Error: UxHwFloatGeneratePath - Function must be called with numberOfParameterStructs greater than 0.\n");

		return;
	}

	if (stateGeneratorFuncPtr == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatGeneratePath - stateGeneratorFuncPtr cannot be NULL.\n");

		return;
	}

	if (parameterStructArray == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatGeneratePath - parameterStructArray cannot be NULL.\n");

		return;
	}

	if (resultArray == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatGeneratePath - resultArray cannot be NULL.\n");

		return;
	}

	for (size_t opt = 0; opt < numberOfParameterStructs; opt++)
	{
		if (resultArray[opt] == NULL)
		{
			fprintf(stderr, "Error: UxHwFloatGeneratePath - resultArray cannot contain NULL.\n");

			return;
		}

		for (size_t step = 0; step < pathLength; step++)
		{
			if (resultArray[opt][step] == NULL)
			{
				fprintf(stderr, "Error: UxHwFloatGeneratePath - resultArray cannot contain NULL.\n");

				return;
			}
		}
	}

	/*
	 *	Loop through each set of parameters
	 */
	for (size_t ii = 0; ii < numberOfParameterStructs; ii++)
	{
		/*
		 *	Call the stepping function the required number of times
		 */
		for (size_t step = 0; step < pathLength; step++)
		{
			stateGeneratorFuncPtr(parameterStructArray[ii], resultArray[ii], step, numberOfStateVariables);
		}
	}

	return;
}

void
UxHwDoubleGeneratePath(
	void * parameterStructArray[],
	size_t numberOfParameterStructs,
	size_t pathLength,
	size_t numberOfStateVariables,
	void (* stateGeneratorFuncPtr)(void * parameterStruct, double * *  paths, size_t iterationStepCount, size_t numberOfStateVariables),
	double * * *  resultArray)
{
	/*
	 *	Sanity checking of input arguments
	 */
	if (pathLength == 0)
	{
		fprintf(stderr, "Error: UxHwDoubleGeneratePath - Function must be called with a pathLength greater than 0.\n");

		return;
	}

	if (numberOfStateVariables == 0)
	{
		fprintf(stderr, "Error: UxHwDoubleGeneratePath - Function must be called with numberOfStateVariables greater than 0.\n");

		return;
	}

	if (numberOfParameterStructs == 0)
	{
		fprintf(stderr, "Error: UxHwDoubleGeneratePath - Function must be called with numberOfParameterStructs greater than 0.\n");

		return;
	}

	if (stateGeneratorFuncPtr == NULL)
	{
		fprintf(stderr, "Error: UxHwDoubleGeneratePath - stateGeneratorFuncPtr cannot be NULL.\n");

		return;
	}

	if (parameterStructArray == NULL)
	{
		fprintf(stderr, "Error: UxHwDoubleGeneratePath - parameterStructArray cannot be NULL.\n");

		return;
	}

	if (resultArray == NULL)
	{
		fprintf(stderr, "Error: UxHwDoubleGeneratePath - resultArray cannot be NULL.\n");

		return;
	}

	for (size_t opt = 0; opt < numberOfParameterStructs; opt++)
	{
		if (resultArray[opt] == NULL)
		{
			fprintf(stderr, "Error: UxHwDoubleGeneratePath - resultArray cannot contain NULL.\n");

			return;
		}

		for (size_t step = 0; step < pathLength; step++)
		{
			if (resultArray[opt][step] == NULL)
			{
				fprintf(stderr, "Error: UxHwDoubleGeneratePath - resultArray cannot contain NULL.\n");

				return;
			}
		}
	}

	/*
	 *	Loop through each set of parameters
	 */
	for (size_t ii = 0; ii < numberOfParameterStructs; ii++)
	{
		/*
		 *	Call the stepping function the required number of times
		 */
		for (size_t step = 0; step < pathLength; step++)
		{
			stateGeneratorFuncPtr(parameterStructArray[ii], resultArray[ii], step, numberOfStateVariables);
		}
	}

	return;
}

void
UxHwFloatGenerateConvergingPath(
	void * parameterStructArray[],
	size_t numberOfParameterStructs,
	size_t maxPathLength,
	size_t numberOfStateVariables,
	bool (* stateGeneratorFuncPtr)(void * parameterStruct, float * *  paths, size_t iterationStepCount, size_t numberOfStateVariables),
	float * *    resultArray)
{
	float * * paths;

	/*
	 *	Sanity checking of input arguments
	 */
	if (maxPathLength == 0)
	{
		fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - Function must be called with a maxPathLength greater than 0.\n");

		return;
	}

	if (numberOfStateVariables == 0)
	{
		fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - Function must be called with numberOfStateVariables greater than 0.\n");

		return;
	}

	if (numberOfParameterStructs == 0)
	{
		fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - Function must be called with numberOfParameterStructs greater than 0.\n");

		return;
	}

	if (stateGeneratorFuncPtr == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - stateGeneratorFuncPtr cannot be NULL.\n");

		return;
	}

	if (parameterStructArray == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - parameterStructArray cannot be NULL.\n");

		return;
	}

	if (resultArray == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - resultArray cannot be NULL.\n");

		return;
	}

	for (size_t opt = 0; opt < numberOfParameterStructs; opt++)
	{
		if (resultArray[opt] == NULL)
		{
			fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - resultArray cannot contain NULL.\n");

			return;
		}
	}

	/*
	 *	Allocate the "working" array
	 */
	paths = (float * *) calloc(maxPathLength, sizeof(float *));

	if (paths == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - Could not allocate memory for buffer.\n");

		return;
	}

	for (size_t step = 0; step < maxPathLength; step++)
	{
		paths[step] = (float *) calloc(numberOfStateVariables, sizeof(float));

		if (paths[step] == NULL)
		{
			fprintf(stderr, "Error: UxHwFloatGenerateConvergingPath - Could not allocate memory for buffer.\n");

			for (size_t jj = 0; jj < step; jj++)
			{
				free(paths[jj]);
			}

			free(paths);

			return;
		}
	}

	/*
	 *	Loop through each set of parameters
	 */
	for (size_t ii = 0; ii < numberOfParameterStructs; ii++)
	{
		/*
		 *	If the function never converges, return the final step
		 */
		size_t lastStep = maxPathLength - 1;

		/*
		 *	Clear the buffer so no values carry over from the previous parameter
		 *	struct. calloc() has already cleared it for the first one.
		 */
		if (ii > 0)
		{
			for (size_t step = 0; step < maxPathLength; step++)
			{
				memset(paths[step], 0, numberOfStateVariables * sizeof(float));
			}
		}

		/*
		 *	Call the stepping function the required number of times
		 */
		for (size_t step = 0; step < maxPathLength; step++)
		{
			if (stateGeneratorFuncPtr(parameterStructArray[ii], paths, step, numberOfStateVariables))
			{
				/*
				 *	Function is indicating it has converged
				 */
				lastStep = step;
				break;
			}
		}

		memcpy(resultArray[ii], paths[lastStep], numberOfStateVariables * sizeof(float));
	}

	for (size_t step = 0; step < maxPathLength; step++)
	{
		free(paths[step]);
	}

	free(paths);

	return;
}

void
UxHwDoubleGenerateConvergingPath(
	void * parameterStructArray[],
	size_t numberOfParameterStructs,
	size_t maxPathLength,
	size_t numberOfStateVariables,
	bool (* stateGeneratorFuncPtr)(void * parameterStruct, double * *  paths, size_t iterationStepCount, size_t numberOfStateVariables),
	double * *   resultArray)
{
	double * * paths;

	/*
	 *	Sanity checking of input arguments
	 */
	if (maxPathLength == 0)
	{
		fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - Function must be called with a maxPathLength greater than 0.\n");

		return;
	}

	if (numberOfStateVariables == 0)
	{
		fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - Function must be called with numberOfStateVariables greater than 0.\n");

		return;
	}

	if (numberOfParameterStructs == 0)
	{
		fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - Function must be called with numberOfParameterStructs greater than 0.\n");

		return;
	}

	if (stateGeneratorFuncPtr == NULL)
	{
		fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - stateGeneratorFuncPtr cannot be NULL.\n");

		return;
	}

	if (parameterStructArray == NULL)
	{
		fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - parameterStructArray cannot be NULL.\n");

		return;
	}

	if (resultArray == NULL)
	{
		fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - resultArray cannot be NULL.\n");

		return;
	}

	for (size_t opt = 0; opt < numberOfParameterStructs; opt++)
	{
		if (resultArray[opt] == NULL)
		{
			fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - resultArray cannot contain NULL.\n");

			return;
		}
	}

	/*
	 *	Allocate the "working" array
	 */
	paths = (double * *) calloc(maxPathLength, sizeof(double *));

	if (paths == NULL)
	{
		fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - Could not allocate memory for buffer.\n");

		return;
	}

	for (size_t step = 0; step < maxPathLength; step++)
	{
		paths[step] = (double *) calloc(numberOfStateVariables, sizeof(double));

		if (paths[step] == NULL)
		{
			fprintf(stderr, "Error: UxHwDoubleGenerateConvergingPath - Could not allocate memory for buffer.\n");

			for (size_t jj = 0; jj < step; jj++)
			{
				free(paths[jj]);
			}

			free(paths);

			return;
		}
	}

	/*
	 *	Loop through each set of parameters
	 */
	for (size_t ii = 0; ii < numberOfParameterStructs; ii++)
	{
		/*
		 *	If the function never converges, return the final step
		 */
		size_t lastStep = maxPathLength - 1;

		/*
		 *	Clear the buffer so no values carry over from the previous parameter
		 *	struct. calloc() has already cleared it for the first one.
		 */
		if (ii > 0)
		{
			for (size_t step = 0; step < maxPathLength; step++)
			{
				memset(paths[step], 0, numberOfStateVariables * sizeof(double));
			}
		}

		/*
		 *	Call the stepping function the required number of times
		 */
		for (size_t step = 0; step < maxPathLength; step++)
		{
			if (stateGeneratorFuncPtr(parameterStructArray[ii], paths, step, numberOfStateVariables))
			{
				/*
				 *	Function is indicating it has converged
				 */
				lastStep = step;
				break;
			}
		}

		memcpy(resultArray[ii], paths[lastStep], numberOfStateVariables * sizeof(double));
	}

	for (size_t step = 0; step < maxPathLength; step++)
	{
		free(paths[step]);
	}

	free(paths);

	return;
}

void
UxHwFloatPropagateFunction(
	void ( * functionPtr )(void * args, float * input, size_t sizeOfInput, float * output, size_t sizeOfOutput),
	void *      args,
	float *     input,
	size_t sizeOfInput,
	float *     output,
	size_t sizeOfOutput)
{
	if (functionPtr == NULL)
	{
		fprintf(stderr, "Error: UxHwFloatPropagateFunction - functionPtr cannot be NULL.\n");

		return;
	}

	functionPtr(args, input, sizeOfInput, output, sizeOfOutput);

	return;
}

void
UxHwDoublePropagateFunction(
	void ( * functionPtr )(void * args, double * input, size_t sizeOfInput, double * output, size_t sizeOfOutput),
	void *      args,
	double *    input,
	size_t sizeOfInput,
	double *    output,
	size_t sizeOfOutput)
{
	if (functionPtr == NULL)
	{
		fprintf(stderr, "Error: UxHwDoublePropagateFunction - functionPtr cannot be NULL.\n");

		return;
	}

	functionPtr(args, input, sizeOfInput, output, sizeOfOutput);

	return;
}
