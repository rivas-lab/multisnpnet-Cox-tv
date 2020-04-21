## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(survival)
library(dplyr)
N = nrow(jasa)
jasa$ID = 1:N
jasa$y = pmax(0.5, as.numeric(jasa$fu.date - jasa$accept.dt))
jasa$age = jasa$age - 48
jasa$year = as.numeric(jasa$accept.dt - as.Date("1967-10-01"))/365.25

jasa$txtime= with(jasa, ifelse(tx.date== fu.date,
                               (tx.date -accept.dt) -.5,
                               (tx.date - accept.dt)))
jasa$status = jasa$fustat

## ---- message=FALSE-----------------------------------------------------------
# at the beginning no-one has transplant
tv_trans = data.frame(ID=1:N, time=0, value = 0)
tmp = jasa %>% select(c(ID, txtime)) %>% filter(!is.na(txtime))
tmp$time = tmp$txtime
tmp$value = 1
tmp = select(tmp, c(ID, time, value))
# At the transplant, the transplant indicator becomes 1
tv_trans = rbind(tv_trans, tmp)

# The interaction  between age and transplant
tv_inter = tv_trans
tv_inter$value = tv_inter$value * jasa$age[match(tv_inter$ID, jasa$ID)]
tv_list = list(transplant = tv_trans, transplant_age = tv_inter)


## ----setup--------------------------------------------------------------------
library(coxtv)
result = coxtv(jasa, tv_list,c('surgery', 'year', 'age'),c(0.01, 0.0))
result

## -----------------------------------------------------------------------------
cindex_tv(jasa, tv_list,c('surgery', 'year', 'age'), result[[1]])
cindex_tv(jasa, tv_list,c('surgery', 'year', 'age'), result[[2]])

## -----------------------------------------------------------------------------
KM_curve_tv(jasa, tv_list[[1]], ngroup=2)

