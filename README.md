### Usage
Compile the code by
```bash
make clean; make
```
Once the compilation is done, an executable named `test_Ylm` will be created. Run the test by
```bash
./test_Ylm <Lmax> <n>
```
`Lmax` is the maximum degree number for `l`, `n` is the number of random coordinates to test. 

> Note that the current tests work only for `Lmax <= 6` since the reference answers for `Lmax > 6` are not available. However, the `sph_harmonics` routine works for any non-negative integer `Lmax`.

For example, the following comand will generate `100000` random coordinates to test the spherical harmonics routines for `l` up to `3` and compare the results.

```bash
./test_Ylm 3 100000
```

The above command gives the following results:

```python
l =  0, m =  0, error_lm = 0.000e+00
l =  1, m = -1, error_lm = 3.331e-16
l =  1, m =  0, error_lm = 4.441e-16
l =  1, m =  1, error_lm = 2.776e-16
l =  2, m = -2, error_lm = 4.441e-16
l =  2, m = -1, error_lm = 4.996e-16
l =  2, m =  0, error_lm = 5.375e-16
l =  2, m =  1, error_lm = 4.441e-16
l =  2, m =  2, error_lm = 1.617e-16
l =  3, m = -3, error_lm = 4.718e-16
l =  3, m = -2, error_lm = 6.661e-16
l =  3, m = -1, error_lm = 9.159e-16
l =  3, m =  0, error_lm = 5.551e-16
l =  3, m =  1, error_lm = 8.604e-16
l =  3, m =  2, error_lm = 2.470e-16
l =  3, m =  3, error_lm = 6.661e-16
Success! All tests passed!
===============
= Timing info =
===============
Run-time of sph_harmonics: 23.172 ms
Run-time of SPARC routine: 4.169 ms
Total run-time of the test: 45.258 ms
```