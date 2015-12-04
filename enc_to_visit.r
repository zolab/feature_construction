####################################################
### Multiple histogram function
####################################################

'
Purpose: Turn distinct physician encounters into contiguous visits. Can be applied for: (1) ED encs, (2) OP encs, and (3) IP encs.
Creator: MJP
Date: 12/4/2015
Note: Slow b/c goes row-by-row.
To do:
[ ] improve speed
'

#############################################################
# turn enc --> visits (need to customize fn if port to new context)
enc.to.visit <- function(dt, adm_date = adm_date, disch_date = disch_date) {
  
  # logical substitutions
  setkey(dt, empi, adm_date, disch_date)
  dt <- dt[, disch_date := max(disch_date), by = c("empi", "adm_date")]  # associating each adm to its most logical (latest) discharge date
  dt <- unique(dt)  # remove duplicates
  dt <- dt[, adm_date := min(adm_date), by = c("empi", "disch_date")]  # associating each disch to its most logical (earliest) admission date
  dt <- unique(dt)
  
  # ! ensure data in right order -- CRUCIAL
  setorder(dt, empi, adm_date, disch_date)
  
  # clean: (1) incorporate adjacent dates into one visit and (2) delete overlapping ones with redundant info
  dt$discard <- 0
  
  
  ## [ ] works, BUT SUPER SLOW --> apply? DT?
  for(i in 2:nrow(dt)) {
    samepat <- dt$empi[i]==dt$empi[i-1]
    curovl <- min(dt$disch_date[i], dt$disch_date[i-1]) - dt$adm_date[i] + 1
    if(curovl == 0 & samepat) {
      dt$adm_date[i] <- dt$adm_date[i-1]  # combine adm_date if adjacent dates but no overlap
      dt$discard[i-1] <- 1  # then mark old (adjacent) for deletion
    }
    if(curovl > 0 & samepat) dt$discard[i] <- 1
  }
  
  # discard overlapping and incorporated periods
  dt <- dt[discard != 1]
  dt[, discard := NULL]
  
  # set key
  setkey(dt, empi, adm_date, disch_date)
  
  return(dt)
  
}
#############################################################
