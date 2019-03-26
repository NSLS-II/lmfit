# lmfit

This is the home page of **lmfit**,
a self-contained C library for Levenberg-Marquardt
least-squares minimization and curve fitting.

# Background information

## History

The core algorithm of lmfit has been invented by K Levenberg (1944) and D W Marquardt (1963).
The seminal FORTRAN implementation of MINPACK
is due to J M Moré, B S Garbow, and K E Hillstrom (1980).
An automatic translation into C was first published by S Moshier.
The present library lmfit started from a similar automatic translation.
With time, every single line was hand edited to convert the code into genuine, readable C.
The API was modified; examples, man pages, and build scripts were added.
Corrections and refinements were contributed by users (see CHANGELOG).

## Mathematics

For a gentle introduction to the Levenberg-Marquardt algorithm,
see e.g. K Madsen, H B Nielsen, O Tingleff (2004):
*Methods for non-linear least squares problems.*

## License

[FreeBSD License](http://opensource.org/licenses/BSD-2-Clause).

## Citation

When using lmfit as a tool in scientific work, no acknowledgement is needed:
If fit results are robust, it does not matter by which software they have been obtained.
If results are not robust, they should not be published anyway.
In methodological publications, e g when describing software that interacts with lmfit,
the preferred form of citation is:

Joachim Wuttke: lmfit – a C library for Levenberg-Marquardt least-squares minimization and curve fitting. Version […], retrieved on […] from https://jugit.fz-juelich.de/mlz/lmfit.

# Usage

* Concise documentation is available in form of manual pages:
  * [lmfit(7)](http://apps.jcns.fz-juelich.de/man/lmfit.html): Summary.
  * [lmmin2(3)](http://apps.jcns.fz-juelich.de/man/lmmin2.html):
      Generic routine for minimizing whatever sum of squares.
  * [lmmin(3)](http://apps.jcns.fz-juelich.de/man/lmmin.html):
      Simpler legacy API without error estimates.
  * [lmcurve(3)](http://apps.jcns.fz-juelich.de/man/lmcurve.html):
      Simpler interface for fitting a parametric function f(x;p) to a data set y(x).
* Sample code:
  * [curve fitting with lmcurve()]()
  * [surface fitting as example for minimization with lmmin()]()
  * [nonlinear equations solving with lmmin()]()
* Application questions:
  * Constraining parameters
  * Weighing data points
* Language bindings:
  * https://github.com/jvail/lmfit.js by Jan Vaillant

# Installation

## From source

Build&install are based on CMake. Out-of-source build is enforced. After unpacking the source, go to the source directory and do:

```
mkdir build
cd build
cmake ..
ctest
make
make install
```