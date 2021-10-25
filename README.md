# WrightFisherSimulation
A R-package that simulates the behavior of population genetics according to the [Wright-Fisher Model](https://en.wikipedia.org/wiki/Genetic_drift#Wright%E2%80%93Fisher_model). It assumes a diploid population with two possible alleles. Using that simulation, grafics can be generated. They show the evolution of different realizations from a given starting point. One focuses on the number of the alleles, another one on the variance between the realizations, and the third one the heterozygosity.


## Installation
You can install the package in R via
```devtools::install_git("https://github.com/johschnee/WrightFisherSimulation.git")```.

## Example of usage
```
library(WrightFisherSimulation)
init(N = 50, X0 = 50, realizations = 50, gens = 400)
p <- simulation()
graphic_individuals(p)
graphic_variance(p)
graphic_heterozygosity(p)
```

## Reference
V. Hösel, C. Kuttler, J. Müller. _Mathematical Population Genetics and Evolution of Bacterial Cooperation_, Word Scientific Press, Singapore, 2020, pp. 26-36.
