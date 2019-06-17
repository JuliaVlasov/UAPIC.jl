# UAPIC.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->


Uniformly Accurate scheme implementation for 4d Vlasov-Poisson PIC simulation

**NOTE: This package is still very much under development and is not fully tested.**

- UA scheme for 4d Vlasov-Poisson in Fluid-scaling with b(x). Update b(x(tn)) every step.


- Fortran reference program `bupdate` in `fortran` directory
```
time ./bupdate

real	0m33.632s
user	0m33.222s
sys	0m0.399s
```

- Julia program
```
julia -O3 --check-bounds=no test/bupdate.jl

108.887904 seconds (30.72 M allocations: 2.646 GiB, 0.81% gc time)
```
