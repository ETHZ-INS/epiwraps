# Inject (insert) values at positions in a vector

Inject (insert) values at positions in a vector

## Usage

``` r
inject(what, inWhat, at)
```

## Arguments

- what:

  The values to inject

- inWhat:

  The vector in which to inject them

- at:

  The positions in the vector \*after\* which to inject the values.

## Value

A vector of same mode at \`inWhat\`, and of length equal to the sum of
the lengths of \`what\` and \`inWhat\`.

## Examples

``` r
inWhat <- 1:10
inject(c(21,22,23), inWhat, at=as.integer(c(0,5,10)))
#>  [1] 21  1  2  3  4  5 22  6  7  8  9 10 23
```
