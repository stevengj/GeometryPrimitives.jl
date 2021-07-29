# GeometryPrimitives

[![Build status](https://github.com/stevengj/GeometryPrimitives.jl/workflows/CI/badge.svg)](https://github.com/stevengj/GeometryPrimitives.jl/actions)
[![codecov.io](http://codecov.io/github/stevengj/GeometryPrimitives.jl/coverage.svg?branch=master)](http://codecov.io/github/stevengj/GeometryPrimitives.jl?branch=master)

This package provides a set of geometric primitive types (spheres, boxes,
cylinders, and so on) and operations on them designed to enable piecewise
definition of functions, especially for finite-difference and finite-element
simulations, in the Julia language.

For example, suppose that you are discretizing a PDE like the Poisson
equation ∇⋅c∇u = f, and you want to provide a simple user interface
for the user to specify the function `c(x)`.  In many applications,
`c` will be piecewise constant, and you want to be able to specify
`c = 1` in one box, `c = 2` in some cylinders, etcetera.   The
GeometryPrimitives package allows the user to provide a list of
shapes with associated data (in this case, the value of `c`) to
define such a `c(x)`.

Furthermore, the application to discretized simulations imposes a couple
of additional requirements:

* One needs to be able to evaluate `c(x)` a huge number of times (once
  for every point on a grid).  So, we provide a fast O(log n) K-D tree
  data structure for rapid searching of shapes.

* Often, one wants to compute the *average* of `c(x)` over a voxel,
  so we provide routines for rapid *approximate* voxel averages.

* Often, one needs not only the value `c(x)` but the normal vector
  to the nearest shape, so we provide normal-vector computation.

This package was inspired by the geometry utilities in my
[Libctl](http://ab-initio.mit.edu/wiki/index.php/Libctl) package.
