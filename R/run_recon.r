#' Standardized density function for generalized normal (Subbotin) distribution
#'
#' Standardized version of \code{dgnorm}.
#' 
#' @param \code{data} data frame with the following \code{colnames}:
#' \describe{
#' \item{\code{brood.yr}}{The year the fish were spawned.}
#' \item{\code{pop}}{Population ID.}
#' \item{\code{p3...pN}}{Multiple columns of the proportion of age-3, age-4, ..., age-N fish.}
#' \item{\code{pHOS}}{The proportion of Hatchery Origin Spawners.}
#' \item{\code{nS}}{The total number of wild + hatchery-origin spawners.}
#' \item{\code{wild.broodstk}}{The number of fish taken for hatchery broodstock.}
#' \item{\code{hrate.w}}{The harvest rate (0-1) of wild fish.}
#' }
#' 
#' @export
run_recon <- function(data)
{
  year <- data$brood.yr
  pop <- data$pop
  q <- sweep(data[,c("p3","p4","p5")], 1, rowSums(data[,c("p3","p4","p5")]), "/")
  pHOS <- data$pHOS
  S_W <- data$nS*(1 - pHOS)
  S_H <- data$nS*pHOS
  B_take <- data$wild.broodstk
  F_rate <- data$hrate.w
  R <- matrix(NA, nrow(data), 3)
  for(i in 1:nrow(data))
    for(age in 3:5)
    {
      if(year[i] + age <= max(year[pop==pop[i]]))
      {
        B_rate <- ifelse(age==3, 0, B_take[i+age]/(S_W[i+age]*sum(q[i+age,-1]) + B_take[i+age]))
        F_rate23 <- ifelse(age==3, 0, F_rate[i+age])
        R[i,age-2] <- S_W[i+age]*q[i+age,age-2]/((1 - B_rate)*(1 - F_rate23))
      }
    }
  return(list(R = R, R_tot = rowSums(R)))
}
