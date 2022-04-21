### Usage
Compile the code by
```bash
make clean; make
```
Once the compilation is done, an executable named `test_Ylm` will be created. Run the test by
```bash
./test_Ylm <Lmax> <n>
```
`Lmax` is the maximum degree number for `l` (note that the current tests work only for `Lmax <= 6` since the reference answers for `Lmax > 6` are not available. However, the `sph_harmonics` routine works for any non-negative integer `Lmax`). `n` is the number of random coordinates to test. For example, the following comand will generate `100000` random coordinates to test the spherical harmonics routines and compare the results.
```bash
./test_Ylm 4 100000
```
Result of the above command gives:
```bash
l =  0, m =  0, error_lm = 0.000e+00
l =  1, m = -1, error_lm = 3.331e-16
l =  1, m =  0, error_lm = 4.441e-16
l =  1, m =  1, error_lm = 2.776e-16
l =  2, m = -2, error_lm = 4.441e-16
l =  2, m = -1, error_lm = 4.996e-16
l =  2, m =  0, error_lm = 5.375e-16
l =  2, m =  1, error_lm = 4.441e-16
l =  2, m =  2, error_lm = 1.617e-16
Success! All tests passed!
===============
= Timing info =
===============
Run-time of sph_harmonics: 22.867 ms
Run-time of SPARC routine: 2.121 ms
Total run-time of the test: 41.810 ms
```