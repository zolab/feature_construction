# ##########################################
# DESCRIPTION: This file defines a function that classifies ed_ccs into categories 
# specified in a csv file (ed_cc_class) - a file which also contains 
# 'cleaning instructions' - and returns a single ed cc per encoutner based upon a 
# user input hierarchy (second argument to the function)

#NOTE: Below an example of how the function is to be called from a main code file
# CALL: 
## ed_cc_class_list <- ed_cc_classification(
## ed_encounters_cohort_control_long_ed_cc,ed_cc_class_hierarchy)
# Where 
# (a) the ed_encounters_cohort_control_long_ed_cc file contains ed cc in a 
# long-format, one row per day-patient (multiple ED visits per day are combined 
# into one row - linked by '---') & an encounter number (row-based)
# (b) the ed_cc_class_hierarchy file contains the ed cc categories ranked in the 
# order in which they are to be imposed 

# OUTPUT
# ed_encounters_cohort_control_long_ed_cc <- 
# ed_encounters_cohort_control_long_ed_cc[ed_cc_class_list$ed_cc_class_short, 
# on=c("encounter_number")]
# This call merges the ed_cc categories back into the original 
# ed_encounters_cohort_control_long_ed_cc file - at the patient-day level

# ###########################################
# Creator: Clara Marquardt
# Date: 23rd October
# ##########################################
# Language: R
# ##########################################
# TO-DO-LIST 
# [] Figure out how to load regular expressions without the current workarounds
# [] Validate results against the Stata code

# NOTE
# This code is to be superseded by a targex based implementation - at some point...
################################################################################

ed_cc_classification <- function(ed_cc_table_input, ed_cc_hierarchy) {

# load the relevant libraries and functions
source("/data/zolab/pants/helpers/30_day_cancer_chemo/pants_cancer_chemo_functions.R")
source("/data/zolab/pants/helpers/30_day_cancer_chemo/pants_cancer_chemo_libraries.R")

# make a copy of the input table so as to leave that unchanged
ed_cc_table <- copy(ed_cc_table_input)

# load and format the csv file with the regular expressions that underlie the 
# categorisation (can be changeg)
ed_cc_class <- fread(
  "/data/zolab/pants/helpers/30_day_cancer_chemo/ed_cc_classification.txt",sep=";")
ed_cc_class <- ed_cc_class[,setdiff(names(ed_cc_class), "gap_1"),with=F]

regx_var <- names(ed_cc_class)

for (var in regx_var) {
  mapply(function(datatable_variable,regx_var_variable,pattern_old_variable, 
  pattern_new_variable) regx_specialcharacters(datatable_variable,regx_var_variable,
  pattern_old_variable, pattern_new_variable),regx_var_variable=var,
  pattern_new_variable=pattern_regx, pattern_old_variable=pattern_input,
  MoreArgs = list(datatable_variable=ed_cc_class))
}

################################################################################
### Replace (conditionally) dirty with clean words (words & separators)

# generate a list of the conditions (3 conditions max per replacement)
ed_cc_class[is.na(word_condition_pos_1), word_condition_pos_1:=""]
ed_cc_class[is.na(word_condition_pos_2), word_condition_pos_2:=""]
ed_cc_class[is.na(word_condition_pos_3), word_condition_pos_3:=""]
ed_cc_class[is.na(word_condition_neg_1)|word_condition_neg_1=="", 
  word_condition_neg_1:="XXX$$$XXX"]
ed_cc_class[is.na(word_condition_neg_2)|word_condition_neg_2=="", 
  word_condition_neg_2:="XXX$$$XXX"]
ed_cc_class[is.na(word_condition_neg_3)|word_condition_neg_3=="", 
  word_condition_neg_3:="XXX$$$XXX"]

word_conditions_pos <- strsplit(paste(ed_cc_class$word_condition_pos_1, 
  ed_cc_class$word_condition_pos_2, ed_cc_class$word_condition_pos_3, sep="---"), 
  "---")
word_conditions_neg <- strsplit(paste(ed_cc_class$word_condition_neg_1, 
  ed_cc_class$word_condition_neg_2, ed_cc_class$word_condition_neg_3, sep="---"), 
  "---")

# generate list of the dirty and clean words
dirty_words <- as.list(ed_cc_class$dirty_words)
clean_words <- as.list(ed_cc_class$cleaned_words)

# define the function that updates the ed_cc column by replacing (conditionally) 
# the dirty with the clean words
update_dt_replace <- function(datatable,pattern_old, pattern_new, condition) {
 datatable[condition,cc_dx:=gsub(pattern_old,pattern_new,cc_dx)]
}

# execute the function by looping over the pattern_old, pattern_new and condition 
# list (joining the conditions using &  & then combining both results)
word_conditions_combined_pos <- lapply(word_conditions_pos, function(x) 
  Reduce(`&`, lapply(x, grepl, ed_cc_table$cc_dx)))
word_conditions_combined_neg <- lapply(word_conditions_neg, function(x) 
  Reduce(`&`, lapply(x, Negate(grepl), ed_cc_table$cc_dx)))
word_conditions_combined <- mapply('+',word_conditions_combined_pos,
  word_conditions_combined_neg,SIMPLIFY=F)

# turn the combined conditions back into a logical class list format
a <- as.relistable(word_conditions_combined)
u <- unlist(a)
u[u==2] <- 100
u[u!=2 & u!=100] <- 0
u[u==100] <- 1
u <- as.logical(u)
word_conditions_combined <- relist(u, a)

mapply(function(datatable_variable,dirty_words_variable,clean_words_variable,
  word_conditions_variable) update_dt_replace(datatable_variable,
  dirty_words_variable,clean_words_variable,word_conditions_variable),
  dirty_words_variable=dirty_words,
  clean_words_variable=clean_words,
  word_conditions_variable=word_conditions_combined,
  MoreArgs = list(datatable_variable=ed_cc_table))

################################################################################
### Split the ed_cc based on the / separator

n <- max(str_count(ed_cc_table$cc_dx,"/"))+1
ed_cc_table[, paste("cc_dx",1:n,sep=""):=tstrsplit(cc_dx,"/")][,cc_dx:=NULL]
ed_cc_table <- as.data.table(gather(ed_cc_table, ed_cc_number, cc_dx, 
  -(empi:ed_dx_code_type)))[!is.na(cc_dx)][, ed_cc_count:=.N, 
  by=c("encounter_number")]

################################################################################
### Generate ed_cc dummies

# generate a list of the conditions (3 conditions max per replacement)
ed_cc_class[is.na(ed_cc_condition_pos_1), ed_cc_condition_pos_1:=""]
ed_cc_class[is.na(ed_cc_condition_pos_2), ed_cc_condition_pos_2:=""]
ed_cc_class[is.na(ed_cc_condition_pos_3), ed_cc_condition_pos_3:=""]
ed_cc_class[is.na(ed_cc_condition_neg_1)|ed_cc_condition_neg_1=="", 
  ed_cc_condition_neg_1:="XXX$$$XXX"]
ed_cc_class[is.na(ed_cc_condition_neg_2)|ed_cc_condition_neg_2==""
  , ed_cc_condition_neg_2:="XXX$$$XXX"]
ed_cc_class[is.na(ed_cc_condition_neg_3)|ed_cc_condition_neg_3=="", 
  ed_cc_condition_neg_3:="XXX$$$XXX"]

ed_cc_conditions_pos <- strsplit(paste(ed_cc_class$ed_cc_condition_pos_1, 
  ed_cc_class$ed_cc_condition_pos_2, ed_cc_class$ed_cc_condition_pos_3, 
  sep="---"), "---")
ed_cc_conditions_neg <- strsplit(paste(ed_cc_class$ed_cc_condition_neg_1, 
  ed_cc_class$ed_cc_condition_neg_2,ed_cc_class$ed_cc_condition_neg_3, 
  sep="---"), "---")

# generate list of ed_cc categories & generate dummies
ed_cc_cat_names <-gsub("^\\s+|\\s+$","",ed_cc_class$ed_cc_cat[
  !ed_cc_class$ed_cc_cat==""])
ed_cc_cat_names_unique <- unique(ed_cc_cat_names)
ed_cc_cat_fill <- rep(0,nrow(ed_cc_table))
ed_cc_table[, c(ed_cc_cat_names_unique) := as.list(rep(list(ed_cc_cat_fill),
  length(ed_cc_cat_names_unique)))]

# define the function that updates(conditionally) the ed_cc_cat dummy 
update_dt_ed_cc_cat <- function(datatable,ed_cc_cat, condition) {
  ed_cc_cat_col <- parse(text=ed_cc_cat)
  datatable[condition,eval(ed_cc_cat_col):=1]
}

# execute the function by looping over the ed_cc_cat and the conditions list 
# (joining the pos and neg conditions using & then combining both results)
ed_cc_conditions_combined_pos <- lapply(ed_cc_conditions_pos, function(x) 
  Reduce(`&`, lapply(x, grepl, ed_cc_table$cc_dx)))
ed_cc_conditions_combined_neg <- lapply(ed_cc_conditions_neg, function(x) 
  Reduce(`&`, lapply(x, Negate(grepl), ed_cc_table$cc_dx)))
ed_cc_conditions_combined <- mapply('+',ed_cc_conditions_combined_pos,
  ed_cc_conditions_combined_neg, SIMPLIFY=F)

# turn the combined conditions back into a logical class list format
a <- as.relistable(ed_cc_conditions_combined)
u <- unlist(a)
u[u==2] <- 100
u[u!=2 & u!=100] <- 0
u[u==100] <- 1
u <- as.logical(u)
ed_cc_conditions_combined <- relist(u, a)

mapply(function(datatable_variable,ed_cc_cat_variable,
  ed_cc_conditions_variable) update_dt_ed_cc_cat(datatable_variable,
  ed_cc_cat_variable,
  ed_cc_conditions_variable),
  ed_cc_cat_variable=ed_cc_cat_names,
  ed_cc_conditions_variable=ed_cc_conditions_combined[1:length(ed_cc_cat_names)],
  MoreArgs = list(datatable_variable=ed_cc_table))

ed_cc_table[,"ed_cc_cat_sum":=rowSums(.SD),with=F,.SDcols=c(ed_cc_cat_names_unique)][
  ed_cc_cat_sum==0, ed_cc_cat_other:=1][is.na(ed_cc_cat_other),ed_cc_cat_other:=0] 
ed_cc_table_temp <- copy(ed_cc_table)


################################################################################
### Impose the hierarchy - at the level of individual ed_cc complaints rather 
### than the ed_enc level

ed_cc_table[,ed_cc_number_row:=.I]
ed_cc_table <- as.data.table(gather(ed_cc_table, ed_cc_cat,ed_cc_cat_dummy,
   -(empi:ed_cc_count),-ed_cc_number_row))
ed_cc_table <- unique(ed_cc_table[ed_cc_hierarchy, nomatch=NA, on=c("ed_cc_cat")][
  order(ed_cc_number,ed_cc_rank)][!ed_cc_cat_dummy==0], by=c("ed_cc_number_row"))

ed_cc_class_long <- ed_cc_table

# SAVE
write.csv(ed_cc_class_long, "ed_cc_class_long.csv", row.names=F)

################################################################################
### Impose the hierarchy - at the level of encounters
ed_cc_table <- unique(ed_cc_table[order(encounter_number,ed_cc_rank)]
, by=c("encounter_number"))

ed_cc_class_short <- ed_cc_table

# SAVE
write.csv(ed_cc_class, "ed_cc_class.csv", row.names=F)


################################################################################
### Output - Two Matrixes - (1) Long - At the level of individual ed_cc complaints 
### and (2) At the level of encounters (patient-day)
ed_cc_class_list<- list("ed_cc_class_long"=ed_cc_class_long, 
  "ed_cc_class_short"=ed_cc_class_short)
return(ed_cc_class_list)

}
################################################################################
#################################   END    #####################################
################################################################################
################################################################################

