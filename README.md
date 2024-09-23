# UxHwCompatibilityForNativeExecution
Stubs for UxHw functions so that developers can compile against UxHw API[^1] on non-uncertainty-tracking hardware.
The compatibility version of the UxHw API implements all available UxHw functions and uses GNU Scientific Library (GSL)[^2] to generate random samples
from parametric distributions. It also uses GSL to draw samples from empirical distributions provided as user input.
You can use `UxHwCompatibilityForNativeExecution` to implement native Monte Carlo implementations in order to quantify the benefit of running an application 
on Signaloid's platform compared to running the same application on conventional architectures, without the need to modify the application source code.

## Usage
UxHw Compatibility API is used as a submodule in many of the Signaloid public application examples.
In order to compile and run natively an application that uses `UxHwCompatibilityForNativeExecution`, you need to:

0. Install dependencies (e.g., on Linux)
```
sudo apt-get install libgsl-dev libgslcblas0
```
1. Use the source files and header files provided by `UxHwCompatibilityForNativeExecution` to compile your program
```
gcc -I/path/to/Signaloid-Demo-UxHwCompatibilityForNativeExecution /path/to/Signaloid-Demo-UxHwCompatibilityForNativeExecution/uxhw.c main.c -o main.bin -lgsl -lgslcblas -lm
```
Don't forget to link with `-lgsl` and `-lgslcblas`.

2. Run the application natively
```
./main.bin ...
```

---

[^1]: [UxHw API Documentation](https://docs.signaloid.io/docs/hardware-api/)
[^2]: [GNU Scientific Library](https://www.gnu.org/software/gsl/)