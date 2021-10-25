#' Initialization
#'
#' Sets some values, which are needed for all other functions. Always run this
#' function in the beginning before using any other function of the package!
#'
#' @param N Number of individuals in the population
#' @param X0 Number of occurrences for that allele you are looking at
#' @param realizations Number of independent runs of the simulaion
#' @param gens Number of generations which should be simulated
#' @examples
#' init(N = 50, X0 = 50, realizations = 50, gens = 400)
#' @export
init <- function(N, X0, realizations, gens){
  X0 <<- X0
  N <<- N
  realizations <<- realizations
  gens <<- gens
  H0 <<- X0/N*(1-(X0-1)/(2*N-1))
}

#' Simulation
#'
#' Creates an array, which contains all necessary information about the genetic
#' development. This array is needed to create the graphics. It can get quite large.
#'
#' All necessary parameters need to be set by using \code{init()}.
#'
#' @examples
#' p <- simulation()
#' @export
simulation <- function(){
  pop <- array(0, dim = c(2*N, realizations, gens))
  start <- rep(c(1,0), c(X0,2*N-X0))
  for (i in 1:realizations) {
    pop[,i,1] <- sample(start, replace = TRUE)
  }
  for (g in 2:gens) {
    for (i in 1:realizations) {
      pop[,i,g] <- sample(pop[,i,g-1], replace = TRUE)
    }
  }
  pop
}

#' Heterozygosity
#'
#' Creates an array, which contains a \code{1} for every heterozygous and a \code{0} for
#' every homozygous individual. The array has a similar structure like that after the
#' simulation, but the first dimension is only the number of individuals.
#'
#' @param pop The array containing the simulation from \code{simulation()}
#' @return An array that indicates whether an individual is heterozygous
#' @examples
#' p <- simulation()
#' h <- heterozygous_individuals(pop = p)
#' @export
heterozygous_individuals <- function(pop) {
  het <- array(0, dim = c(N, realizations, gens))
  for (n in 1:N) {
    for (i in 1:realizations) {
      for (g in 1:gens) {
        if(pop[2*n-1, i, g] != pop[2*n, i, g]) het[n,i,g] <- 1
      }
    }
  }
  het
}

#' Graphic realizations
#'
#' Creates a graphic that shows the evolution of the simulated populations. The mean of
#' the number of alleles is added.
#'
#' @param pop Data from the function \code{simulation()}
#' @param limit The number of realizations that is included in the graphic as a maximum.
#' The mean is always built over all realizations. Default: -1, no limit.
#' @return The described graphic
#' @examples
#' p <- simulation()
#' h <- graphic_individuals(pop = p)
#' @import ggplot2
#' @export
graphic_individuals <- function(pop, limit=-1) {
  evolution_array <- apply(pop, MARGIN = c(3,2), sum)
  evo_mean <- apply(evolution_array, MARGIN = 1, mean)
  evolution <- tibble::as_tibble(evolution_array, .name_repair = "universal")
  cols <- names(evolution)
  evolution <- dplyr::mutate(.data = evolution, gen = 1:dim(evolution_array)[1], mean = evo_mean, .before=1)
  graph <- ggplot(data = evolution)
  l <- 0
  for (c in cols) {
    if(l == limit) break
    graph <- graph + geom_line(mapping = aes(x=gen, y=!!sym(c)))
    l <- l+1
  }
  graph <- graph + geom_line(mapping = aes(x=gen, y=mean, colour="red")) +
    labs(x = "generation n", y = expression("Number of alles X"[n]),
         title = "Time Course of the Simulation",
         subtitle = "Number of alleles same the the starting alleles for (maybe just some) realizations.\nBesides that, the mean over all realizations is shown in red.") +
    theme(legend.title = element_blank()) +
    scale_color_hue(labels = c("mean"))
  graph
}

#' Graphic variance
#'
#' Creates a graphic that shows the variance of the number of alleles over generations.
#'
#' @param pop Data from the function \code{simulation()}
#' @return The described graphic
#' @examples
#' p <- simulation()
#' h <- graphic_variance(pop = p)
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @export
graphic_variance <- function(pop) {
  var_gens <- tibble::as_tibble(apply(apply(pop, MARGIN = c(3,2), sum), MARGIN = 1, var), .name_repair = "universal") %>%
    dplyr::mutate(gen = 1:gens)
  ggplot(data = var_gens) +
    stat_function(fun = function(.x) X0*(2*N-X0)) +
    geom_point(mapping = aes(x=gen, y=value)) +
    labs(x = "generation n", y = expression("Variance of X"[n]),
         title = expression("Variance of X"[n]),
         subtitle = "The points show the Variance in every generation. The line marks the calculated limit.")
}

#' Graphic heterozygousity
#'
#' Creates a graphic that shows the heterozygosity of the number of alleles over
#' generations.
#'
#' @param pop Data from the function \code{simulation()}
#' @return The described graphic
#' @examples
#' p <- simulation()
#' h <- graphic_heterozygosity(pop = p)
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @export
graphic_heterozygosity <- function(pop) {
  het <- heterozygous_individuals(pop)
  het_mean <- tibble::as_tibble(apply(het, MARGIN = 3, sum)/(N*realizations), .name_repair = "universal") %>%
    dplyr::mutate(gen = 1:gens)
  ggplot(data.frame(x = c(0,gens+1)), aes(x)) +
    stat_function(fun = function(.x) H0*(1-1/(2*N))^(.x)) +
    geom_point(data = het_mean, mapping = aes(x = gen, y=value)) +
    labs(x = "generation n", y = expression("Heterozygosity in generation n"),
         title = "Heterozygosity in the course of time",
         subtitle = "The line shows the theoretical calculated heterozygosity. The points visualize\nthe real heterozygosity in the simulation.")
}
