# salmonIPM

This is the development repo for `salmonIPM`, an `R` package for fitting integrated population models to salmon data. Click [here](https://ebuhle.github.io/salmonIPM) for the full vignette.

## Installation

You can install the development version using `devtools`.

```{r }
if(!require("devtools")) {
  install.packages("devtools")
  library("devtools")
}
devtools::install_github("ebuhle/salmonIPM")
```
