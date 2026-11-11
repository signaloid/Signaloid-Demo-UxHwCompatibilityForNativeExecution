# UxHwCompatibilityForNativeExecution
Stubs for UxHw® functions so that developers can compile against UxHw API[^1] on non-uncertainty-tracking hardware.
The compatibility version of the UxHw API implements all available UxHw functions and uses GNU Scientific Library (GSL)[^2] to generate random samples
from parametric distributions. It also uses GSL to draw samples from empirical distributions provided as user input.
You can use `UxHwCompatibilityForNativeExecution` to run native Monte Carlo simulations that quantify the benefit of running an application
on Signaloid's platform compared to running the same application on conventional architectures, without the need to modify the application source code.

## Usage
UxHw Compatibility API is used as a submodule in the Signaloid public application examples.
In order to compile and run natively an application that uses `UxHwCompatibilityForNativeExecution`, you need to:

0. Install dependencies (e.g., on Linux)
```
sudo apt-get install libgsl-dev libgslcblas0
```
1. Use the source files and header files provided by `UxHwCompatibilityForNativeExecution` to compile your program
```
gcc -I/path/to/Signaloid-Demo-UxHwCompatibilityForNativeExecution /path/to/Signaloid-Demo-UxHwCompatibilityForNativeExecution/uxhw.c main.c -o main.bin -lgsl -lgslcblas -lm
```
Do not forget to link with `-lgsl` and `-lgslcblas`.

2. Run the application natively
```
./main.bin ...
```

## Random number generation
By default the GSL generator is seeded from the current time, so each run gives different samples.
To make runs repeatable, set the `GSL_RNG_SEED` environment variable.
```
GSL_RNG_SEED=1 ./main.bin ...
```
You can also choose the generator with `GSL_RNG_TYPE`, for example `GSL_RNG_TYPE=ranlxd2`. The default is `mt19937`.
When either variable is set, GSL prints its value to stderr.

## Behavior in native execution
Native execution does not track uncertainty. Each variable holds a single value, and the functions behave as follows.
- Distribution functions, such as `UxHwDoubleGaussDist`, return one random sample. Invalid parameters return NaN.
- `UxHwDoubleDistFromSamples` returns one of the input samples, picked at random. `UxHwDoubleDistFromWeightedSamples` picks by weight.
- `UxHwDoubleMixture(a, b, aScale)` returns `a` with probability `aScale`, and `b` otherwise.
- `UxHwDoubleSample`, `UxHwDoubleSupportMin`, `UxHwDoubleSupportMax`, `UxHwDoubleGetIndependentCopy`, `UxHwDoubleNthMode` and `UxHwDoubleQuantile` return the value itself. `UxHwDoubleSampleBatch` fills the array with the value.
- `UxHwDoubleNthMoment` returns 1 for n = 0, the value for n = 1, and 0 for higher n.
- `UxHwDoubleProbabilityGT` returns 1 if the value is greater than the cutoff, and 0 otherwise.
- `UxHwDoubleEvaluatePDF` returns infinity at the value, and 0 elsewhere.
- `UxHwDoubleLimitDistributionSupport` returns the value if it lies within the limits, and NaN otherwise.
- `UxHwDoubleBayesLaplace` returns the prior and prints a warning.
- `UxHwDoubleDerivative` returns NaN.

Each function has a `Float` version that behaves the same way.

---

[^1]: [UxHw API Documentation](https://docs.signaloid.io/docs/hardware-api/)
[^2]: [GNU Scientific Library](https://www.gnu.org/software/gsl/)