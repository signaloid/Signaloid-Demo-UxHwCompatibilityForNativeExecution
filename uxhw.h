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
#include <stddef.h>
#include <stdint.h>

/*
 *	Link to UxHw API documentation: https://docs.signaloid.io/docs/hardware-api/.
 */

#define SignaloidParticleModifier	""

/*
 *	Data structures for specifying a weighted sample
 */
typedef struct
{
	float		sample;
	float		sampleWeight;
} WeightedFloatSample;

typedef struct
{
	double		sample;
	double		sampleWeight;
} WeightedDoubleSample;

#ifdef __cplusplus
extern "C"
{
#endif

double	UxHwDoubleSample(double value);
float	UxHwFloatSample(float value);

void	UxHwDoubleSampleBatch(double value, double *  destSampleArray, size_t numberOfRandomSamples);
void	UxHwFloatSampleBatch(float value, float *  destSampleArray, size_t numberOfRandomSamples);

double	UxHwDoubleDistFromSamples(double *  samples, size_t sampleCount);
float	UxHwFloatDistFromSamples(float *  samples, size_t samplesCount);

double	UxHwDoubleDistFromWeightedSamples(WeightedDoubleSample *  samples, size_t sampleCount, size_t unweightedSampleCount);
float	UxHwFloatDistFromWeightedSamples(WeightedFloatSample *  samples, size_t sampleCount, size_t unweightedSampleCount);

void	UxHwDoubleDistFromMultidimensionalSamples(double *  destinationArray, void *  samples, size_t sampleCount, size_t sampleCardinality);
void	UxHwFloatDistFromMultidimensionalSamples(float *  destinationArray, void *  samples, size_t sampleCount, size_t sampleCardinality);

double	UxHwDoubleExponentialDist(double mu);
float	UxHwFloatExponentialDist(float mu);

double	UxHwDoubleUniformDist(double min, double max);
float	UxHwFloatUniformDist(float min, float max);

double	UxHwDoubleGaussDist(double mean, double stddev);
float	UxHwFloatGaussDist(float mean, float stddev);

double	UxHwDoubleGumbel1Dist(double mu, double beta);
float	UxHwFloatGumbel1Dist(float mu, float beta);

double	UxHwDoubleLogisticDist(double location, double scale);
float	UxHwFloatLogisticDist(float location, float scale);

double	UxHwDoubleLaplaceDist(double mu, double b);
float	UxHwFloatLaplaceDist(float mu, float b);

double	UxHwDoubleWeibullDist(double k, double lambda);
float	UxHwFloatWeibullDist(float k, float lambda);

double	UxHwDoubleLognormalDist(double mu, double sigma);
float	UxHwFloatLognormalDist(float mu, float sigma);

double	UxHwDoubleBoundedparetoDist(double alpha, double min, double max);
float	UxHwFloatBoundedparetoDist(float alpha, float min, float max);

double	UxHwDoubleMixture(double dist1, double dist2, double p);
float	UxHwFloatMixture(float dist1, float dist2, float p);

double	UxHwDoubleNthMoment(double value, size_t n);
float	UxHwFloatNthMoment(float value, size_t n);

double	UxHwDoubleNthMode(double value, size_t n);
float	UxHwFloatNthMode(float value, size_t n);

double	UxHwDoubleSupportMin(double value);
float	UxHwFloatSupportMin(float value);

double	UxHwDoubleSupportMax(double value);
float	UxHwFloatSupportMax(float value);

double	UxHwDoubleProbabilityGT(double value, double cutoff);
float	UxHwFloatProbabilityGT(float value, float cutoff);

double	UxHwDoubleCopyDistShape(double value);
float	UxHwFloatCopyDistShape(float value);

double	UxHwDoubleLimitDistributionSupport(double value, double supportMin, double supportMax);
float	UxHwFloatLimitDistributionSupport(float value, float supportMin, float supportMax);

double	UxHwDoubleQuantile(double value, double quantileProbability);
float	UxHwFloatQuantile(float value, float quantileProbability);

double	UxHwDoubleBayesLaplace(double (*likelihood)(double), double prior, double evidence);
float	UxHwFloatBayesLaplace(float (*likelihood)(float), float prior, float evidence);

#define min(a, b)	((a) < (b) ? (a) : (b))
#define max(a, b)	((a) > (b) ? (a) : (b))

#ifdef __cplusplus
}
#endif
