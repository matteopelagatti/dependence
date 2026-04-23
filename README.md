# An independence test for any kind of variables

The main function of the package is `indeptest(x, y, ...)` that implements two similar tests of independence.
The vectors `x` and `y` contain the values of the two variables for which we want to test the independence.
`x` and `y` can be numerical or categorial (factor, in R terminology). The tests use the same rationality
regardless of the nature of `x` and `y`.

The theoretical details of the test will be pubished later. 

Here are a couple of examples.

Numerical vs. numerical
```
> indeptest(rnorm(100), rnorm(100))

Pelagatti-Monti independence test
  Variable types  : numeric vs. numeric 
  Basis type      : poly 
  Basis dimensions: p = 2, q = 2
  Sample size     : n = 100
  Pillai stat.    : 2.9207  (p-value = 0.5712)
  Bartlett stat.  : 2.8563  (p-value = 0.5822)
```

Factor vs. numerical
```
> indeptest(sample(factor(letters[1:5]), 100, T), rnorm(100))

Pelagatti-Monti independence test
  Variable types  : factor vs. numeric 
  Basis type      : poly(factor) poly(numeric) 
  Basis dimensions: p = 4, q = 2
  Sample size     : n = 100
  Pillai stat.    : 4.0111  (p-value = 0.8561)
  Bartlett stat.  : 3.9011  (p-value = 0.8659)
```

Factor vs. factor
```
> indeptest(sample(factor(letters[1:5]), 100, T), sample(factor(letters[1:5]), 100, T))

Pelagatti-Monti independence test
  Variable types  : factor vs. factor 
  Basis type      : poly 
  Basis dimensions: p = 4, q = 4
  Sample size     : n = 100
  Pillai stat.    : 13.3019  (p-value = 0.6506)
  Bartlett stat.  : 13.0464  (p-value = 0.6694)
```
