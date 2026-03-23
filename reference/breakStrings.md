# breakStrings

breaks a string of long-enough words (or vector thereof) into two lines.

## Usage

``` r
breakStrings(x, minSizeForBreak = 20, lb = "\n")
```

## Arguments

- x:

  A character (or factor) vector.

- minSizeForBreak:

  The minimum number of characters to be broken into two lines.

- lb:

  The separation character (e.g. line break).

## Value

A character vector of length=length(x).

## Examples

``` r
breakStrings("this is too long for practical purposes")
#>    this is too long for practical purposes 
#> "this is too long for\npractical purposes" 
```
