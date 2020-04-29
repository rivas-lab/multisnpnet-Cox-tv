#' Fit a proportional hazard model with time-dependent covariates
#'
#' @param phe a data frame or data table with the following columns:
#'            y: time-to-event;
#'            status: a binary variable, 1 if event  happened, 0 if censored;
#'            time independent covariates specified in ti_vars
#'            ID: the ID of each individual
#' @param ti_vars the names of time-independent predictors
#' @param tv_list a list of data frame or data table, each corresponding to
#'                one time-varying covariate. The elements of this list must
#'                have the following columns:
#'                ID: same ID in phe
#'                time: time when the measurement taken relative to t0, where
#'                time-to-event is defined as t0 -> t1
#'                value: the value of the measurement
#'                Note that the user is responsible to ensure that each ID must
#'                have a measurement before the first event time
#' @param lambda a vector of regularization parameters
#' @param B0 the initial parameter, zero if not provided
#' @param p.fac penalty factor for each variable, uniform if not provided
#' @return The lasso solution
#' @export
coxtv = function(phe, tv_list, ti_vars, lambda, B0=NULL, p.fac=NULL, info=NULL)
{
  if(is.null(info))
  {
    info = get_info(phe, tv_list)
  }
  X = as.matrix(dplyr::select(phe, ti_vars))
  if(any(is.na(X))){
    stop("NA in predictors detected.")
  }
  pti = ncol(X)
  ptv = ncol(info$V)

  if(is.null(B0)){
    B0 = numeric(pti+ptv)
  }

  if(is.null(p.fac)){
    p.fac = numeric(pti+ptv) + 1.0
  }
  result = fit(X,
               info$V,
               info$Vind,
               info$order-1L,
               info$erank,
               info$eind,
               info$te_count,
               1.0,
               B0,
               lambda,
               p.fac,
               2000,
               1.1,
               1e-5)

  return(lapply(result, function(x){names(x)= c(info$names, ti_vars);x}))
}

#' Compute the concordance index of the fit
#' @param B The fitted parameter, all the other variables are the same as above
#' @return The C-index
#' @export
cindex_tv = function(phe, tv_list, ti_vars, B, info=NULL)
{
  if(is.null(info))
  {
    info = get_info(phe, tv_list)
  }
  X = as.matrix(dplyr::select(phe, ti_vars))
  if(any(is.na(X))){
    stop("NA in predictors detected.")
  }
  pti = ncol(X)
  ptv = ncol(info$V)

  eta_td = info$V %*% B[1:ptv]
  eta = X %*% B[(ptv+1):(ptv+pti)]
  eta = eta[info$order]
  CCIndex(eta, eta_td, info$Vind, info$rank_all, info$erank_wt, info$eranke, info$te_count)
}


#' Draw a Kaplan Meier curve based on a time-varying covariate
#' @param tv The time varying covariate. Must satisfy the same condition on the
#'           elements of tv_list above
#' @return The KM-curve
#' @export
KM_curve_tv = function(phe, tv, ngroup=2)
{
  info = get_info(phe, list(tv))
  if(all(tv$value %in% c(0, 1))){
    bound = c(-1, 0.5, 2)
    ngroup = 2
  } else{
    q = seq(0, 1, length.out = ngroup+1)
    q = q[-c(1, length(q))]
    bound = c(-Inf, quantile(tv$value, q), Inf)
  }
  result = matrix(nrow=(1+length(info$erank)), ncol=ngroup)
  result[1,] = 1
  for(i in 1:length(info$erank)){
    start = info$Vind[i]+1
    end = info$Vind[i+1] # start and end of risk set

    values = info$V[start:end, 1]
    event_values = values[1:(info$te_count[i])]
    denoms = sapply(1:ngroup, function(j){sum( (values > bound[j]) & (values <= bound[j+1]))})
    numerators =sapply(1:ngroup, function(j){sum( (event_values > bound[j]) & (event_values <= bound[j+1]))})
    result[(i+1),] = result[i,] *(1 - (numerators/denoms))
  }
  result = data.frame(t=c(0, info$event_time), result)
  colnames(result) = c("t", paste("group", 1:ngroup, sep=""))
  result = data.frame(result[1], stack(result[2:(ngroup+1)]))
  g = ggplot2::ggplot() + ggplot2::geom_step(data=result, mapping=ggplot2::aes(x=t, y=values, color=ind), size=1) +
    ggplot2::xlab("Time") + ggplot2::ylab("Survival probability")+ ggplot2::theme(legend.position="top")

  return(g)
}

#' @export
get_info = function(phe, tv_list){
  ## Get Vind first
  o = order(phe$y, -phe$status, phe$ID) # Use lexicographic order(event time, event, ID)
  eind = which(as.logical(phe$status)) - 1L
  # These two are sorted
  y = phe$y[o]
  status = phe$status[o]
  if(any(is.na(y))){
    stop("NA in time-to-event detected")
  }
  if(any(is.na(status))){
    stop("NA in status detected")
  }
  rank_all = rank(y, ties.method = 'min')
  erank_wt = rank_all[as.logical(status)] # event rank that includes tied events
  te_count = as.numeric(get_te_count(erank_wt))
  event_rank = unique(erank_wt)
  event_time = y[event_rank]
  vind = numeric(length(event_rank)+1)
  vind[1] = 0
  N = length(y)
  for(i in 1:length(event_rank)){
    n = N - event_rank[i] + 1 # size of the risk set
    vind[i+1] = vind[i] + n
  }
  event_rank = event_rank - 1L
  erank_wt = erank_wt - 1L
  rank_all = rank_all - 1L
  vind = as.integer(vind)

  V = matrix(nrow = vind[length(vind)], ncol=length(tv_list))
  for(i in (1:length(tv_list))){
    tv = tv_list[[i]]
    if(!all(phe$ID %in% tv$ID)){
      stop("Some people does not have any time varying covariate measured.")
    }
    #tv$time = tv$age - phe$t0[match(tv$ID, phe$ID)] # time when the measurements are taken, relative to MI
    tv = dplyr::filter(tv, ID %in% phe$ID)
    tv$etime = phe$y[match(tv$ID, phe$ID)] # need this for the ordering
    tv$status = phe$status[match(tv$ID, phe$ID)] # need this for the ordering too
    tv_o = order(tv$etime, -tv$status, tv$ID, tv$time) #lexicalgraphic ordering (event time, event, ID, measuremen time)
    tv = tv[tv_o,]

    tmp = unique(tv$ID)
    tmp = match(tmp, tv$ID)
    if(!all(tv$time[tmp] < min(event_time)))
    {
      num_violation = sum((tv$time[tmp] >= min(event_time)))
      warning(paste(num_violation,
                    "people do not have time-varying covariates measured before",
                    "the first event.",
                    "The most recent measurement after the event is used."))
    }
    tmp =  tmp - 1L
    cumu_measurement = c(tmp, nrow(tv))
    v = numeric(vind[length(vind)])
    get_V(v, tv$value, tv$time, event_time, cumu_measurement, event_rank, vind)
    V[,i] = v
  }
  info = list()
  info$order = o
  info$V =V
  info$Vind = vind
  info$rank_all = rank_all
  info$erank = event_rank
  info$erank_wt = erank_wt
  info$eind = eind
  info$te_count = te_count
  info$eranke = match(erank_wt, event_rank) - 1L
  info$event_time = event_time
  info$names = names(tv_list)
  if(is.null(info$names))
  {
    info$names = paste0("TV", 1:length(tv_list))
  }
  return(info)
}

#' @export
cox_residual = function(phe, ti_vars, info, beta)
{
  X = as.matrix(dplyr::select(phe, ti_vars))
  get_residual(X,
               info$V,
               info$Vind,
               info$order-1L,
               info$erank,
               info$eind,
               info$te_count,
               beta)
}
