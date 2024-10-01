/*
 *	Copyright (c) 2023–2024, Signaloid.
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
#include <time.h>
#include <stddef.h>
#include <sys/time.h>
#include "uxhw.h"

gsl_rng *	gGSLr;

/**
 *	@brief	Initializes the libGSL random number generator. Uses current time to initialize the seed.
 *
 */
static void
initializeGenerators(void)
{
	struct timeval		t;
	const gsl_rng_type *	gslRngType = gsl_rng_default;

	gettimeofday(&t, NULL);
	srandom(t.tv_usec);

	gsl_rng_env_setup();
	gGSLr = gsl_rng_alloc(gslRngType);
	gettimeofday(&t, NULL);
	gsl_rng_set(gGSLr, t.tv_usec);

	return;
}

/**
 *	@brief	Returns a random integer in the interval [0,`n`).
 *
 * 	@param	n	: The upper bound of the interval of the possible integers to be returned.
 * 	@return	int	: Returns an integer `i`, where 0 <= `i` < `n`.
 */
static int
randomFromRange(int n)
{
	int	limit;
	int	r;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	limit = RAND_MAX - (RAND_MAX % n);
	while ((r = random()) >= limit) ;

	return r % n;
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

	for (size_t i = 0; i < numberOfRandomSamples; i++)
	{
		destSampleArray[i] = value;
	}

	return;
}

void
UxHwFloatSampleBatch(float value, float *  destSampleArray, size_t numberOfRandomSamples)
{
	if (destSampleArray == NULL)
	{
		fprintf(stderr, "UxHwDoubleSampleBatch: destSampleArray is NULL.\n");

		return;
	}

	for (size_t i = 0; i < numberOfRandomSamples; i++)
	{
		destSampleArray[i] = value;
	}

	return;
}

double
UxHwDoubleDistFromSamples(double* samples, size_t samplesCount)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return samples[randomFromRange(samplesCount)];
}

float
UxHwFloatDistFromSamples(float *  samples, size_t samplesCount)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return samples[randomFromRange(samplesCount)];
}

double
UxHwDoubleDistFromWeightedSamples(WeightedDoubleSample *  samples, size_t sampleCount, size_t unweightedSampleCount)
{
	double *	normalizedCumulativeWeights;
	double		weightSum = 0.0;
	double		probability;
	size_t		index;

	normalizedCumulativeWeights = (double *)malloc(sampleCount * sizeof(double));
	if (normalizedCumulativeWeights == NULL)
	{
		fprintf(stderr, "UxHwDoubleDistFromWeightedSamples: malloc failed. Returning NAN...");

		return NAN;
	}

	for (size_t i = 0; i < sampleCount; i++)
	{
		weightSum += samples[i].sampleWeight;
		normalizedCumulativeWeights[i] = weightSum;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	probability = gsl_ran_flat(gGSLr, 0.0, 1.0);

	for (size_t i = 0; i < sampleCount; i++)
	{
		normalizedCumulativeWeights[i] /= weightSum;

		if (probability <= normalizedCumulativeWeights[i])
		{
			index = i;
			break;
		}
		
	}

	free(normalizedCumulativeWeights);

	return samples[index].sample;
}

float
UxHwFloatDistFromWeightedSamples(WeightedFloatSample *  samples, size_t sampleCount, size_t unweightedSampleCount)
{
	float *		normalizedCumulativeWeights;
	float		weightSum = 0.0;
	float		probability;
	size_t		index;

	normalizedCumulativeWeights = (float *)malloc(sampleCount * sizeof(float));
	if (normalizedCumulativeWeights == NULL)
	{
		fprintf(stderr, "UxHwFloatDistFromWeightedSamples: malloc failed. Returning NAN...");

		return NAN;
	}

	for (size_t i = 0; i < sampleCount; i++)
	{
		weightSum += samples[i].sampleWeight;
		normalizedCumulativeWeights[i] = weightSum;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	probability = gsl_ran_flat(gGSLr, 0.0, 1.0);

	for (size_t i = 0; i < sampleCount; i++)
	{
		normalizedCumulativeWeights[i] /= weightSum;

		if (probability <= normalizedCumulativeWeights[i])
		{
			index = i;
			break;
		}
		
	}

	free(normalizedCumulativeWeights);

	return samples[index].sample;
}

void
UxHwDoubleDistFromMultidimensionalSamples(double *  destinationArray, void *  samples, size_t sampleCount, size_t sampleCardinality)
{
	double **	castSamples = (double **)samples;
	size_t		randomIndex;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	randomIndex = randomFromRange(sampleCount);

	for (size_t i = 0; i < sampleCardinality; i++)
	{
		destinationArray[i] = castSamples[randomIndex][i];
	}

	return;
}

void
UxHwFloatDistFromMultidimensionalSamples(float *  destinationArray, void *  samples, size_t sampleCount, size_t sampleCardinality)
{
	float **	castSamples = (float **)samples;
	size_t		randomIndex;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	randomIndex = randomFromRange(sampleCount);

	for (size_t i = 0; i < sampleCardinality; i++)
	{
		destinationArray[i] = castSamples[randomIndex][i];
	}

	return;
}

double
UxHwDoubleExponentialDist(double mu)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_exponential(gGSLr, mu);
}

float
UxHwFloatExponentialDist(float mu)
{
	return (float)UxHwDoubleExponentialDist((double)mu);
}

double
UxHwDoubleGumbel1Dist(double mu, double beta)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_gumbel1(gGSLr, 1.0 / beta, 1.0) + mu;
}

float
UxHwFloatGumbel1Dist(float mu, float beta)
{
	return (float)UxHwDoubleGumbel1Dist((double)mu, (double)beta);
}

double
UxHwDoubleUniformDist(double min, double max)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return min + ((max - min) * ((double)random() / (double)RAND_MAX));
}

float
UxHwFloatUniformDist(float min, float max)
{
	return (float)UxHwDoubleUniformDist((double)min, (double)max);
}

double
UxHwDoubleGaussDist(double mean, double stddev)
{
	double	sample;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	sample = gsl_ran_gaussian_ziggurat(gGSLr, stddev);

	return mean + sample;
}

float
UxHwFloatGaussDist(float mean, float stddev)
{
	return (float)UxHwDoubleGaussDist((double)mean, (double)stddev);
}

double
UxHwDoubleLogisticDist(double location, double scale)
{
	double sample;

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
	return (float)UxHwDoubleLogisticDist((double)location, (double)scale);
}

double
UxHwDoubleLaplaceDist(double mu, double b)
{
	double sample;

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
	return (float)UxHwDoubleLaplaceDist((double)mu, (double)b);
}

double
UxHwDoubleWeibullDist(double k, double lambda)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_weibull(gGSLr, lambda, k);
}

float
UxHwFloatWeibullDist(float k, float lambda)
{
	return (float)UxHwDoubleWeibullDist((double)k, (double)lambda);
}

double
UxHwDoubleLognormalDist(double mu, double sigma)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return gsl_ran_lognormal(gGSLr, mu, sigma);
}

float
UxHwFloatLognormalDist(float mu, float sigma)
{
	return (float)UxHwDoubleLognormalDist((double)mu, (double)sigma);
}

double
UxHwDoubleBoundedparetoDist(double alpha, double min, double max)
{
	double	sample = max;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	while (sample >= max)
	{
		sample = gsl_ran_pareto(gGSLr, alpha, min);
	}

	return sample;
}

float
UxHwFloatBoundedparetoDist(float alpha, float min, float max)
{
	return (float)UxHwDoubleBoundedparetoDist((double)alpha, (double)min, (double)max);
}

double
UxHwDoubleMixture(double dist1, double dist2, double p)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	if (random() < p * RAND_MAX)
	{
		return dist1;
	}

	return dist2;
}

float
UxHwFloatMixture(float dist1, float dist2, float p)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	if (random() < p * RAND_MAX)
	{
		return dist1;
	}

	return dist2;
}

double
UxHwDoubleNthMoment(double value, size_t n)
{
	if (n == 1)
	{
		return value;
	}

	return 0.0;
}

float
UxHwFloatNthMoment(float value, size_t n)
{
	if (n == 1)
	{
		return value;
	}

	return 0.0;
}

double
UxHwDoubleNthMode(double value, size_t n)
{
	return value;
}

float
UxHwFloatNthMode(float value, size_t n)
{
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
	return (value > cutoff ? 1 : 0);
}

float
UxHwFloatProbabilityGT(float value, float cutoff)
{
	return (value > cutoff ? 1 : 0);
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
UxHwDoubleQuantile(double value, double quantileProbability)
{
	return value;
}

float
UxHwFloatQuantile(float value, float quantileProbability)
{
	return value;
}

double
UxHwDoubleBayesLaplace(double (*likelihood)(double), double prior, double evidence)
{
	fprintf(stderr, "Warning: UxHwDoubleBayesLaplace is not supported in native execution mode! Returning prior...");
	
	return prior;
}

float
UxHwFloatBayesLaplace(float (*likelihood)(float), float prior, float evidence)
{
	fprintf(stderr, "Warning: UxHwFloatBayesLaplace is not supported in native execution mode! Returning prior...");
	
	return prior;
}
