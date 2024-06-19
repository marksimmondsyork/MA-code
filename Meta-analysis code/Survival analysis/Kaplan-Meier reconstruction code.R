#' Augment a raw number at risk table with the necessary information to run
#' the reconstruction algorithm.
#'
#' @param raw_NAR A data frame with the columns 'time' and NAR' at least.
#' @param raw_surv A data frame with the columns 'time' and 'survival' at least.
#' @param tau End of follow-up time, defaults to last time in NAR table.
#'
#' @return An augmented tab that can be used as input in KM_reconstruct().
#'
#' @export
#'
#' @examples
#' data(TTfields_pfs_trt_NAR)
#' data(TTfields_pfs_trt_clicks)
#' augmented_NAR <- format_NAR_tab(rawNAR=TTfields_pfs_trt_NAR, rawSurv=TTfields_pfs_trt_clicks)
#'
#'
#'
format_raw_tabs <- function(raw_NAR, raw_surv, tau=NULL) {
  
  # check clicks tab has correct columns
  has_col <- length(which(colnames(raw_surv) %in% c('time', 'survival')))
  if (has_col != 2) { stop('raw_surv must have columns named time and survival exactly') }
  
  # subset and order clicks
  raw_surv <- dplyr::select(raw_surv, time, survival) %>%
    dplyr::arrange(time)
  
  # make sure survival is non-increasing, starts with (0,1), ends with an event
  if (is.unsorted(rev(raw_surv$survival))) {
    stop('survival must be non-increasing in time') }
  if (raw_surv[1, 1] != 0 | raw_surv[1,2] != 1) {
    stop('Your raw_clicks_tab did not have a t=0,S=1 row, check your work') }
  last_click_row <- nrow(raw_surv)
  last_click_t <- raw_surv$time[last_click_row]
  last_surv <- raw_surv$survival[last_click_row]
  if (last_click_t <= raw_surv$time[last_click_row-1] |
      last_surv >= raw_surv$survival[last_click_row-1]) {
    stop('Your last click should have been at the end of a vertical (not
         horizontal) segment')
  }
  
  # check NAR tab has correct columns
  has_col <- length(which(colnames(raw_NAR) %in% c('time', 'NAR')))
  if (has_col != 2) { stop('raw_NAR must have columns named time and NAR exactly') }
  
  # subset and order NAR
  raw_NAR <- dplyr::select(raw_NAR, time, NAR) %>%
    dplyr::arrange(time)
  
  # make sure NAR is non-increasing
  if (is.unsorted(rev(raw_NAR$NAR))) {
    stop('NAR must be non-increasing in time') }
  
  # follow-up end is the last NAR time unless otherwise specified (e.g. surv goes to 0)
  if (is.null(tau)) {tau=max(raw_NAR$time)}
  
  # match NAR intervals with raw_clicks rows - remember last row done manually
  ints <- data.frame(lower=rep(NA, nrow(raw_NAR)-1), upper=NA)
  for (int_idx in 1:nrow(ints)) {
    temp_rows <- which(raw_surv$time >= raw_NAR$time[int_idx] &
                         raw_surv$time < raw_NAR$time[int_idx+1])
    if (length(temp_rows) == 0) {
      next
    } else {
      ints$lower[int_idx] <- min(temp_rows)
      ints$upper[int_idx] <- max(temp_rows)
    }
  }
  
  # augment NAR, remove NA rows
  aug_NAR <- dplyr::bind_cols(raw_NAR[-nrow(raw_NAR), ], ints) %>%
    dplyr::filter(!is.na(lower))
  
  # manually add last row to NAR and clicks tables
  last_NAR_row <- data.frame(time=max(raw_NAR$time),
                             NAR=min(raw_NAR$NAR), lower=aug_NAR$upper[nrow(aug_NAR)]+1,
                             upper=aug_NAR$upper[nrow(aug_NAR)]+1)
  last_surv_row <- data.frame(time=tau, survival=last_surv)
  
  aug_NAR <- dplyr::bind_rows(aug_NAR, last_NAR_row)
  aug_surv <- dplyr::bind_rows(raw_surv, last_surv_row)
  
  return(list(aug_NAR=aug_NAR, aug_surv=aug_surv))
  
  }







#' Reconstruct individual-level data from augmented survival table and
#' NAR table, with augmentation performed by format_raw_tabs().
#'
#' @param aug_NAR A data frame processed through format_raw_tabs().
#' @param aug_surv A data frame processed through format_raw_tabs().
#'
#' @return A list including IPD_time, IPD_event, n_hat=n_hat,
#' KM_hat, n_cen, n_event, int_censor
#'
#' @export
#'
#' @examples
#' data(TTfields_pfs_trt_NAR)
#' data(TTfields_pfs_trt_clicks)
#' augmented_NAR <- format_NAR_tab(rawNAR=TTfields_pfs_trt_NAR, rawSurv=TTfields_pfs_trt_clicks)
#' KM_reconstruct(aug_NAR=augmented_NAR$aug_NAR, aug_surv=augmented_NAR$aug_surv)
#'
#'
#'
#'
KM_reconstruct <- function(aug_NAR, aug_surv) {

    # info from NAR table
    TAR <- aug_NAR$time
    NAR <- aug_NAR$NAR
    lower <- aug_NAR$lower
    upper <- aug_NAR$upper

    # make sure the time/survival is nonincreasing
    t_surv <- aug_surv$time
    surv <- aug_surv$surv
    if ( is.unsorted(t_surv) | is.unsorted(rev(surv)) ) {stop('aug_surv unsorted')}

    # number of intervals
    total_ints <- length(NAR)
    # number of event times
    total_e_times <- upper[total_ints]

    # number censored on each interval (except last)
    int_censor <- rep(0, total_ints-1)
    # last value of t where we had an event
    last_event <- rep(1, total_ints)

    # estimated subjects remaining at each k
    n_hat <- rep(NAR[1]+1, total_e_times)
    # number censored at each k
    n_cen <- rep(0, total_e_times)
    # number of events at each k
    n_event <- rep(0, total_e_times)
    # S(t) at each k
    KM_hat <- rep(1, total_e_times)

    # loop through intervals
    for (int_idx in 1:(total_ints-1)) {

        # it's possible that surv[lower[int_idx]] = 0 if the KMC goes to 0
        if (surv[lower[int_idx]] == 0) {
            int_censor[int_idx] <- 0
        } else {
            # first approximation of no. censored on interval int_idx
            int_censor[int_idx] <- round(NAR[int_idx] * surv[lower[int_idx+1]] /
                                             surv[lower[int_idx]] - NAR[int_idx+1])
        }

        # adjust int_censor[int_idx] until n_hat = NAR at the start of the next interval
        # if we have too many events, then just add more censoring
        # if too few events and no. censored > 0, then remove censoring
        # if too few events and no.censored <=0, then stuck, just move on
        while ( n_hat[lower[int_idx+1]] > NAR[int_idx+1] |
                (n_hat[lower[int_idx+1]] < NAR[int_idx+1]&&int_censor[int_idx]>0) ) {

            # can't have negative censoring
            if (int_censor[int_idx] <= 0) {
                n_cen[lower[int_idx]:upper[int_idx]] <- 0
                int_censor[int_idx] <- 0
            } else {

                # evenly distribute censoring times
                cen_times <- t_surv[lower[int_idx]] + (1:int_censor[int_idx])*
                    (t_surv[lower[int_idx+1]] - t_surv[lower[int_idx]]) / (int_censor[int_idx]+1)
                n_cen[lower[int_idx]:upper[int_idx]] <- hist(cen_times,
                    breaks=t_surv[lower[int_idx]:lower[int_idx+1]], plot=F)$counts
            }

            # now account for all events in the interval
            n_hat[lower[int_idx]] <- NAR[int_idx]
            last <- last_event[int_idx]
            for (click_idx in lower[int_idx]:upper[int_idx]) {
                # initial row
                if (click_idx == 1) {
                    n_event[click_idx] <- 0
                    KM_hat[click_idx] <- 1
                } else {
                    # have to check if our KMC goes to zero
                    if (KM_hat[last] == 0) {
                        n_event[click_idx] <- 0
                        KM_hat[click_idx] <- 0
                    } else {
                        # KM_hat and S are ideally the same, but since we are rounding/estimating,
                        # there will be small differences
                        n_event[click_idx] <- round(n_hat[click_idx] * (1-(surv[click_idx] / KM_hat[last])))
                        KM_hat[click_idx] <- KM_hat[last] * (1-(n_event[click_idx] / n_hat[click_idx]))
                    }
                }

                # fill in next n_hat
                n_hat[click_idx+1] <- n_hat[click_idx] - n_event[click_idx] - n_cen[click_idx]
                # update last
                if (n_event[click_idx] != 0) {last <- click_idx}
            }

            # update amount of censoring we need
            int_censor[int_idx] <- int_censor[int_idx] + (n_hat[lower[int_idx+1]] - NAR[int_idx+1])
        } # end while loop through one interval

        # if ended the interval with fewer estimated at risk than published number,
        # that means there int_censor[int_idx] was not positive, so nobody else to redistribute,
        # so we need to change the number at risk from the published number to our estimated
        # number before continuing the estimation
        if (n_hat[lower[int_idx+1]] < NAR[int_idx+1]) {NAR[int_idx+1] <- n_hat[lower[int_idx+1]]}

        last_event[int_idx+1] <- last
    } # end looping through intervals

    # record events and event times
    IPD_event <- rep(0, NAR[1])
    IPD_event[1:sum(n_event)] <- 1
    IPD_time <- c()
    for (click_idx in 1:total_e_times) {
        IPD_time <- c(IPD_time, rep(t_surv[click_idx], n_event[click_idx]))
    }

    # record censoring times
    for (click_idx in 1:(total_e_times-1)) {
        IPD_time <- c(IPD_time, rep((t_surv[click_idx]+t_surv[click_idx+1])/2, n_cen[click_idx]))
    }

    # fill rest of events as censored at max(t_surv)
    ppl_remain <- length(IPD_event) - length(IPD_time)
    if (ppl_remain < 0) {
        stop('Algorithm failed, ended up with too many people')
    } else {
        IPD_time <- c(IPD_time, rep(max(t_surv), ppl_remain))
    }

    # return separated results
    return( list(IPD_time=IPD_time, IPD_event=IPD_event, n_hat=n_hat,
           KM_hat=KM_hat, n_cen=n_cen, n_event=n_event, int_censor=int_censor) )
}
