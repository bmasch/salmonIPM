#' Perform run reconstruction on brood table data.
#'
#' @param fish_data Data frame that includes the following \code{colnames}, in no particular order except where noted:
#' \describe{
#' \item{\code{pop}}{Numeric or character population ID.}
#' \item{\code{year}}{Numeric variable giving the year the fish spawned (i.e., the brood year).}
#' \item{\code{A}}{Spawning habitat size (either stream length or area). Will usually be time-invariant within a population, but need not be.}
#' \item{\code{S_tot_obs}}{Total number (not density) of wild and hatchery-origin spawners.}
#' \item{\code{n_age_minAge...n_age_maxAge}}{Multiple columns of observed spawner age frequencies (i.e., counts), where minAge (maxAge) is the numeral age in years (total, not ocean age) of the youngest (oldest) spawners.}
#' \item{\code{n_W_obs}}{Observed frequency of natural-origin spawners.}
#' \item{\code{n_H_obs}}{Observed frequency of hatchery-origin spawners.}
#' \item{\code{F_rate}}{Total harvest rate (proportion) of natural-origin fish.}
#' \item{\code{B_take_obs}}{Number of adults taken for hatchery broodstock.}
#' }
#' 
#' @return A data frame with the following \code{colnames}, some of which are simply replicated from fish_data:
#' \describe{
#' \item{\code{pop}}{See above.}
#' \item{\code{year}}{See above.}
#' \item{\code{A}}{See above.}
#' \item{\code{S}}{Same as S_tot_obs above.}
#' \item{\code{p_age_minAge...p_age_maxAge}}{Multiple columns of spawner age \emph{relative} frequencies corresponding to the frequencies in fish_data.}
#' \item{\code{p_HOS}}{Proportion of hatchery-origin spawners.}
#' \item{\code{F_rate}}{See above.}
#' \item{\code{B_take_obs}}{See above.}
#' #' \item{\code{R}}{Total natural-origin recruits from the brood year in each row.}
#' }
#' 
#' @export
run_recon <- function(fish_data)
{
  with(fish_data, {
    age_range <- as.numeric(substring(names(fish_data)[grep("n_age", names(fish_data))], 6, 6))
    nn <- fish_data[,grep("n_age", names(fish_data))]
    qq <- sweep(nn, 1, rowSums(nn), "/")
    for(i in unique(pop))
      if(any(pop==i & rowSums(nn)==0))
        qq[pop==i & rowSums(nn)==0,] <- rep(colMeans(qq[pop==i & rowSums(nn) > 0,]),
                                            each = sum(pop==i & rowSums(nn)==0))
    substr(names(qq), 1, 1) <- "p" 
    p_HOS <- n_H_obs/(n_H_obs + n_W_obs)
    recon_dat <- cbind(pop = pop, A, year = year, S = S_tot_obs, qq, p_HOS = p_HOS, 
                       F_rate = F_rate, B_take_obs = B_take_obs, R = NA)
    R <- matrix(NA, nrow(recon_dat), length(age_range))
    
    for(i in 1:nrow(recon_dat))
    {
      for(j in 1:length(age_range))
      {
        a <- age_range[j]
        if(year[i] + a <= max(year[pop==pop[i]]))
        {
          B_rate <- ifelse(j==1, 0, B_take_obs[i+a]/(S_tot_obs[i+a]*(1 - p_HOS[i+a])*(1 - qq[i+a,1]) + B_take_obs[i+a]))
          F_eff <- ifelse(j==1, 0, F_rate[i+a])
          R[i,j] <- S_tot_obs[i+a]*(1 - p_HOS[i+a])*qq[i+a,j]/((1 - B_rate)*(1 - F_eff))
        }
      }
    }
    
    recon_dat$R <- rowSums(R)
    return(recon_dat)
  })
}
