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
#include <stdio.h>
#include <time.h>
#include <stddef.h>
#include <sys/time.h>
#include "uxhw.h"

gsl_rng * gGSLr;

/**
 *	@brief	Initializes the libGSL random number generator. Uses current time to initialize the seed.
 *
 */
static void
initializeGenerators(void)
{
	struct timeval          t;
	const gsl_rng_type *    gslRngType = gsl_rng_default;

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
	int limit;
	int r;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	limit = RAND_MAX - (RAND_MAX % n);

	while ((r = random()) >= limit)
	{
		/* Rejection sampling: discard values that would cause modulo bias. */
	}

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
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return samples[randomFromRange(sampleCount)];
}

float
UxHwFloatDistFromSamples(float * samples, size_t sampleCount)
{
	if (!gGSLr)
	{
		initializeGenerators();
	}

	return samples[randomFromRange(sampleCount)];
}

double
UxHwDoubleDistFromWeightedSamples(WeightedDoubleSample * samples, size_t weightedSampleCount)
{
	double *    normalizedCumulativeWeights;
	double      weightSum = 0.0;
	double      probability;
	size_t      index;

	normalizedCumulativeWeights = (double *) malloc(weightedSampleCount * sizeof(double));

	if (normalizedCumulativeWeights == NULL)
	{
		fprintf(stderr, "UxHwDoubleDistFromWeightedSamples: malloc failed. Returning NAN...");

		return NAN;
	}

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		weightSum += samples[ii].sampleWeight;
		normalizedCumulativeWeights[ii] = weightSum;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	probability = gsl_ran_flat(gGSLr, 0.0, 1.0);

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		normalizedCumulativeWeights[ii] /= weightSum;

		if (probability <= normalizedCumulativeWeights[ii])
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
	float   probability;
	size_t  index;

	normalizedCumulativeWeights = (float *) malloc(weightedSampleCount * sizeof(float));

	if (normalizedCumulativeWeights == NULL)
	{
		fprintf(stderr, "UxHwFloatDistFromWeightedSamples: malloc failed. Returning NAN...");

		return NAN;
	}

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		weightSum += samples[ii].sampleWeight;
		normalizedCumulativeWeights[ii] = weightSum;
	}

	if (!gGSLr)
	{
		initializeGenerators();
	}

	probability = gsl_ran_flat(gGSLr, 0.0, 1.0);

	for (size_t ii = 0; ii < weightedSampleCount; ii++)
	{
		normalizedCumulativeWeights[ii] /= weightSum;

		if (probability <= normalizedCumulativeWeights[ii])
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
	double * *  castSamples = (double * *) samples;
	size_t      randomIndex;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	randomIndex = randomFromRange(sampleCount);

	for (size_t ii = 0; ii < sampleCardinality; ii++)
	{
		destinationArray[ii] = castSamples[randomIndex][ii];
	}

	return;
}

void
UxHwFloatDistFromMultidimensionalSamples(float * destinationArray, void *  samples, size_t sampleCount, size_t sampleCardinality)
{
	float * *   castSamples = (float * *) samples;
	size_t      randomIndex;

	if (!gGSLr)
	{
		initializeGenerators();
	}

	randomIndex = randomFromRange(sampleCount);

	for (size_t ii = 0; ii < sampleCardinality; ii++)
	{
		destinationArray[ii] = castSamples[randomIndex][ii];
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
	return (float) UxHwDoubleExponentialDist((double) mu);
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
	return (float) UxHwDoubleGumbel1Dist((double) mu, (double) beta);
}

double
UxHwDoubleUniformDist(double a, double b)
{
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

	if (!gGSLr)
	{
		initializeGenerators();
	}

	while (sample >= xMax)
	{
		sample = gsl_ran_pareto(gGSLr, alpha, xMin);
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

	if (random() < aScale * RAND_MAX)
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

	if (random() < aScale * RAND_MAX)
	{
		return a;
	}

	return b;
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
	for (size_t ii = 0; ii < numberOfDistributions; ii++)
	{
		destDistArray[ii] = srcDistArray[ii];
	}

	return;
}

void
UxHwFloatGetIndependentJointMultidimensionalCopy(float * srcDistArray, float *  destDistArray, size_t numberOfDistributions)
{
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
	return value;
}

float
UxHwFloatQuantile(float value, float probability)
{
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
	fprintf(stderr, "Warning: UxHwDoubleBayesLaplace is not supported in native execution mode! Returning prior...");

	return prior;
}

float
UxHwFloatBayesLaplace(float ( * statisticalModel )(void *, float), void * modelParams, float prior, float observedData, size_t numberOfObservations)
{
	fprintf(stderr, "Warning: UxHwFloatBayesLaplace is not supported in native execution mode! Returning prior...");

	return prior;
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
		minimizingArgumentInstanceArray);
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
		minimizingArgumentInstanceArray);
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
	 *	Sanity checking	of input arguments
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
	 *	Sanity checking	of input arguments
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
	 *	Sanity checking	of input arguments
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

			return;
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
		for (size_t step = 0; step < maxPathLength; step++)
		{
			if (stateGeneratorFuncPtr(parameterStructArray[ii], paths, step, numberOfStateVariables))
			{
				/*
				 *	Function is indicating it has converged
				 */
				resultArray[ii] = paths[step];

				break;
			}

			if (step == maxPathLength - 1)
			{
				/*
				 *	Function has failed to converge - return the final step
				 */
				resultArray[ii] = paths[step];
			}
		}
	}

	for (size_t step = 0; step < maxPathLength; step++)
	{
		free(paths[step]);
	}

	free(paths);
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
	 *	Sanity checking	of input arguments
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

			return;
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
		for (size_t step = 0; step < maxPathLength; step++)
		{
			if (stateGeneratorFuncPtr(parameterStructArray[ii], paths, step, numberOfStateVariables))
			{
				/*
				 *	Function is indicating it has converged
				 */
				for (size_t jj = 0; jj < numberOfStateVariables; jj++)
				{
					resultArray[ii][jj] = paths[step][jj];
				}

				for (size_t step = 0; step < maxPathLength; step++)
				{
					free(paths[step]);
				}

				free(paths);

				break;
			}

			if (step == maxPathLength - 1)
			{
				/*
				 *	Function has failed to converge - return the final step
				 */
				for (size_t jj = 0; jj < numberOfStateVariables; jj++)
				{
					resultArray[ii][jj] = paths[step][jj];
				}
			}
		}
	}

	for (size_t step = 0; step < maxPathLength; step++)
	{
		free(paths[step]);
	}

	free(paths);
}

void
UxHwFloatPropagateFunction(
	void ( * functionPtr )(void * args, float *  input, size_t sizeOfInput, float *  output, size_t sizeOfOutput),
	void *      args,
	float *     input,
	size_t sizeOfInput,
	float *     output,
	size_t sizeOfOutput)
{
	functionPtr(args, input, sizeOfInput, output, sizeOfOutput);

	return;
}

void
UxHwDoublePropagateFunction(
	void ( * functionPtr )(void * args, double *  input, size_t sizeOfInput, double *  output, size_t sizeOfOutput),
	void *      args,
	double *    input,
	size_t sizeOfInput,
	double *    output,
	size_t sizeOfOutput)
{
	functionPtr(args, input, sizeOfInput, output, sizeOfOutput);

	return;
}
