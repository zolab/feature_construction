 ##########################################
# DESCRIPTION: 
# This file 
# (a) classifies a subset of microbiology results (Blood specimen based 
# results & a subset or urine based results (LEGIONELLA ANTIGEN, 
# STREPTOCOCCUS PNEUMONIAE ANTIGEN)) - a susbet of results which serves as the 
# 'gold standard' for infections. The code classifies the results into a set of 
# granular categories (invalid/mabigious, unambigiously negative, unambigiously 
# positive - tier 3, unambigiously positive - tier 2),unambigiously positive - 
 # tier 1) which are in a final step 
# aggregated into two higher-level categories (inavlid/ambigious/negative and 
# positive). The division into positive tier 1, 2 and 3 results is 
# determined by clinical considerations - all are unmabigiously 'positive' 
# results but only tier 1 test results are deemed to be unmabigious evidence of 
# an infection (tier 2,3 results are therefore in the final step classfied as 
# negative). (b) bacterial cultures identified as positive tier 1 results are in 
# a second step classified by the genus of the identified organisms by drawing on 
# the taxize R package (with some modifications to account for missidentifications)
# The code is designed to be modular in that (i) both granular and aggregated 
# categoreis are reported allowing for an easy reclassification and (ii) the 
# regular expressions underlying the code are contained in a cvs file that 
# can easily be modified
# ##########################################
# Creator: Clara Marquardt
# Date: 15th December
# ##########################################
# Language: R
# ##########################################
# TO-DO-LIST 
# [] Add better/more informative in line comments
# [] Validate - organism classification & micro classification
# [] Optimise
# [] FIgure out how to read in regx expressions without the workarounds
# [] Resolve last few instances of regular expressions living in the code rather
# than the csv file
# [] Ensure that the code is in fact entirely modular with respect to the 
# taxonomy level that microorganisms are to be identified at
# [] Re-think organism classification 
# [] Apply to other test results in Mic file - braoden coade to also classify 
# these

################################################################################
######################  (A) SET-UP - 'Control File'  ###########################
################################################################################
### operational specifications
# (a) working directory
setwd("/data/zolab/pants") 
modified_data_folder <- "modified_data/"
output_data_folder <- "output/"
temp_data_folder <- "temp/"

# (b) libraries and functions
source("/data/zolab/pants/helpers/30_day_cancer_chemo/pants_cancer_chemo_functions.R")
source("/data/zolab/pants/helpers/30_day_cancer_chemo/pants_cancer_chemo_libraries.R")

# (c) raw data files
micro_filename <- "/data/zolab/bwh_ed/dock/rpdr_2010_2012_Mic.txt"

# (d) regular expression related files
org_exceptions_filename <- 
  "/data/zolab/pants/helpers/30_day_cancer_chemo/spec_org_name.csv"
classification_regx_filename <- 
  "/data/zolab/pants/helpers/30_day_cancer_chemo/bc_regx_classification.txt"

### inputs
#(a) level - taxonmy hierarchy - at which to classify microbiology organisms
### NOTE: These have to correspond to the levels represented in the texonmy 
# databases that the taxise package draws on
tex_level <- "genus"

################################################################################
######################  (A) SET-UP - Load Files ################################
################################################################################
# (a) load raw micro files & format
mic_raw <- setkey(fread(micro_filename , 
  stringsAsFactors=F, sep="$", header=T),empi)[,microbiology_time:=
  format(strptime(sapply(strsplit(microbiology_time, " "), "[", 2), "%H:%M:%S"),
  "%H:%M")][, ':='(empi=as.integer(empi))]

# (b) load regx related files & format
# NOTE: The formatting here ensures that the regx expressions are correctly read 
# in (this step currently involves a number of rather convoluted workarounds which
# will need to be resolved)
spec_org_name_raw <-setkey(fread(org_exceptions_filename),spec_org_name)
spec_org_name <- paste(ifelse(spec_org_name_raw$spec_org_name_2!="",
  sprintf("(?<!no )%s,%s",spec_org_name_raw$spec_org_name, 
  spec_org_name_raw$spec_org_name_2),sprintf("(?<!no )%s",
  spec_org_name_raw$spec_org_name)),collapse="|")
spec_org_name <- gsub("_s_p_a_c_e_", " ", spec_org_name)

bc_regx_raw <- fread(classification_regx_filename,sep=";")
regx_var <- names(bc_regx_raw)
for (var in regx_var) {
  mapply(function(datatable_variable,regx_var_variable,pattern_old_variable, 
    pattern_new_variable) regx_specialcharacters(datatable_variable,
    regx_var_variable,pattern_old_variable, pattern_new_variable),
    regx_var_variable=var,pattern_new_variable=pattern_regx, 
    pattern_old_variable=pattern_input,
    MoreArgs = list(datatable_variable=bc_regx_raw))
}

col_number <- ncol(bc_regx_raw)
for (i in 1:col_number) {
  temp <- paste(bc_regx_raw[!(bc_regx_raw[[i]]=="")][[i]], collapse="|")
  category_name <- names(bc_regx_raw)[i]
  assign(category_name, temp)
  write.csv(category_name, 
  sprintf("/data/zolab/pants/helpers/30_day_cancer_chemo/%s_mod.csv", 
  	category_name), row.names=F)
}


################################################################################
###  (B) EXTRACT THE RELEVANT RESULTS & FORMAT THE MIC FILE FOR THE ANALYSIS ###
################################################################################
################################################################################
### (a) extract all valid blood and urine-specimen based results 

# extract all blood/subset of urine specimen based results
mic_inf <- setkey(mic_raw[specimen_type %like% "BLOOD" | 
  specimen_type %like% "URINE" & test_name %in% c("URINE FOR LEGIONELLA ANTIGEN", 
  "STREPTOCOCCUS PNEUMONIAE ANTIGEN (URINE)", "Legionella Urine AG", 
  "Legionella Urinary Antigen","LEGIONELLA ANTIGEN")], empi)

# omit non-interestign results (TEL/TELE results - report reporting of results to
# clinicians, empty results (no empi) & duplicates (on all variables other than 
# organism text) - assume that with the omission of Tel/Tele results and all 100% 
# duplicates each row now presents one unique result (pending elimiantio of 
# invalid, etc. results] 
mic_inf <- unique(mic_inf[empi!="." & !(organism_code %like% "TEL")],
  by=names(mic_inf)[names(mic_inf)!="organism_text"])

#################################################################################
### (b) format file for the regular expression analysis (focus on: specimen type, 
### organism_name, organism_comment, test_name, test_comments)

# define the variables relevant to the regx analysis
relevant_vars <- c("specimen_type", "organism_name","organism_comment",
  "test_name","test_comments")

# transform all the relevant variables to lower case 
regx_formatting(mic_inf, relevant_vars)

# SAVE
write.csv(mic_inf, paste0(temp_data_folder, "mic_inf_temp_", "1", ".csv"),
	row.names=F) 

################################################################################
#########################  (B) CLASSIFY THE RESULTS ############################
################################################################################
################################################################################
### (a) identify invalid and/or ambigious_a results

# main analysis
mic_inf[, valid := ifelse(grepl(invalid_ambigious_results_a,organism_name,perl=T) 
  | grepl(invalid_ambigious_results_a,organism_comment,perl=T) | 
  grepl(invalid_ambigious_results_a,test_comments,perl=T),"invalid/ambigious_a",
  "valid")]

# exceptions
invalid_ambigious_results_a_exception_i <- "this is a corrected report"
invalid_ambigious_results_a_exception_ii <- 
  "^positive|^negative|^antibody|^no antibody|^borderline|^mid"

mic_inf[valid=="invalid/ambigious_a" & 
  (grepl(invalid_ambigious_results_a_exception_i ,organism_name,perl=T) | 
  grepl(invalid_ambigious_results_a_exception_i,organism_comment,perl=T) | 
  grepl(invalid_ambigious_results_a_exception_i, test_comments,perl=T)) & 
  (grepl(invalid_ambigious_results_a_exception_ii ,organism_name,perl=T) | 
  grepl(invalid_ambigious_results_a_exception_ii,organism_comment,perl=T) | 
  grepl(invalid_ambigious_results_a_exception_ii,test_comments,perl=T)), 
  valid := "valid"]

# SAVE
write.csv(mic_inf, paste0(temp_data_folder, "mic_inf_temp_", "2", ".csv"),
	row.names=F) 

################################################################################
### (b) identify unambigiously pos vs. unambigiously negative vs. non-neg results

# analysis - unambigiously negative results (main analysis & exception)
mic_inf[,result := ifelse(!(grepl(neg_results,organism_name,perl=T) ) &
  !(grepl("negative",organism_name,perl=T) & !(grepl(spec_org_name,
  organism_name,perl=T))) &
  !(grepl(neg_results,test_comments,perl=T)) &
  !(grepl("negative",test_comments,perl=T) & !(grepl(spec_org_name,
  test_comments,perl=T) & !(grepl("tests suggested: ebv panel negative",
  test_comments,perl=T) ))),"non-neg","neg")]

# analysis - unambigiously positive (main analysis)
mic_inf[grepl(pos_results,test_comments,perl=T),result:="pos"]

# SAVE
write.csv(mic_inf, paste0(temp_data_folder, "mic_inf_temp_", "3", ".csv"),
	row.names=F) 

################################################################################
### (c) identify non_neg vs. unambigiously negative (i.e. pos_neg)

# analysis - classify observations with no organism name as pos_neg 
mic_inf[organism_name=="" & test_comments=="" & test_name!="" & result=="non-neg",
  result:="posneg"]

# SAVE
write.csv(mic_inf, paste0(temp_data_folder, "mic_inf_temp_", "4", ".csv"),
	row.names=F) 

#################################################################################
### (d) identify invalid and/or ambigious_b vs. pos & update 
### results for all invalid/ambigious_a results

# main analysis -  invalid and/or ambigious_b
mic_inf[organism_name=="" & test_comments=="" & test_name=="" & result=="non-neg",
  result:="invalid/ambigious_b"]

mic_inf[(grepl(invalid_ambigious_results_b,test_comments,perl=T) | 
  grepl(invalid_ambigious_results_b,organism_name,perl=T)) & 
  result=="non-neg",result:="invalid/ambigious_b"]

# exceptions -  invalid and/or ambigious_b
mic_inf[(grepl("(?<!viable for )susceptibility testing",organism_name,perl=T) & 
  result=="non-neg"),result:="invalid/ambigious_b"]
mic_inf[organism_name=="consult lab if further identification required" & 
  result=="non-neg",result:="invalid/ambigious_b"]

mic_inf[result=="invalid/ambigious_b" & grepl("resistance due",organism_name), 
  result := "non-neg"]
mic_inf[result=="invalid/ambigious_b" & grepl("seen",organism_name) & 
  !grepl("failed to grow",organism_name), result := "non-neg"]
mic_inf[result=="invalid/ambigious_b" & grepl(": carbapenem non-susceptible",
  organism_name), result := "non-neg"]
mic_inf[result=="invalid/ambigious_b" & grepl("speciation",organism_name), 
  result := "non-neg"]
mic_inf[result=="invalid/ambigious_b" & grepl("organism_isolated",organism_name),  
  result := "non-neg"]
mic_inf[result=="invalid/ambigious_b" & grepl("yeast|salmonella|streptococcus|
  gram negative bacilli|gram variable pleo rods|corynebacterium group d-2|
  gram positive rods|gram negative rods",organism_name),result := "non-neg"]

# main analysis - unmabigiously positive
mic_inf[result=="non-neg",result:="pos"]

# update result for invalid results
mic_inf[valid=="invalid/ambigious_a", result := "invalid/ambigious_a"]

# SAVE
write.csv(mic_inf, paste0(temp_data_folder, "mic_inf_temp_", "5", ".csv"),
	row.names=F) 

################################################################################
### (e) categorise all positive results into three tiers 

# tier 1
mic_inf[result=="pos" & test_name=="" & test_comments=="" | result=="pos" &  
  specimen_type %like% "urine" & test_name %in% c("urine for legionella antigen", 
  "streptococcus pneumoniae antigen (urine)", "legionella urine ag", 
  "legionella urinary antigen","legionella antigen"),result:="pos_tier_1"]

# tier 2 - note: no longer include pcr as separate category as all 
# results that contain pcr also contain viral load - check this
mic_inf[result=="pos" & grepl(pos_tier_2_results, test_name,perl=T),
  result:="pos_tier_2"]

# tier 3 - other positive
mic_inf[result=="pos",result:="pos_other"]

# SAVE
write.csv(mic_inf, paste0(temp_data_folder, "mic_inf_temp_", "6", ".csv"),
	row.names=F) 

################################################################################
####  (D) CLASSIFY THE ORGANISM NAMES OF ALL BACTERIAL CULTURES IN TIER 1 ######
################################################################################
################################################################################
### (a) classify the tier 2 results 
mic_inf[result=="pos_tier_2" & test_name %like% "igm", result_mod_tier_pos:="igm"]
mic_inf[result=="pos_tier_2" & test_name %like% "pcr|viral load", 
  result_mod_tier_pos:="pcr & viral load"]
mic_inf[result=="pos_tier_2" & test_name %like% "beta-d-glucan|cytomegalovirus", 
  result_mod_tier_pos:="fungal markers"]

### (b) classify the urine results in tier 1 (group into two tests)
mic_inf[result=="pos_tier_1" & specimen_type %like% "urine", result_mod_tier_pos:=
  ifelse(test_name %like% "legionella", "urine legionella antigen", 
  "urine streptococcus pneumoniae antigen")]

### (c) CLASSIFY THE BLOOD CULTURES (TIER 1) - level of classifcation genus
# extract organism names
organism_name_trimmed <- mic_inf[result=="pos_tier_1" & specimen_type %like% 
  "blood"]$organism_name    

# generate table of organism names (split first three words into words)
n <- sapply(gregexpr("\\S+",organism_name_trimmed ), length)
for (i in 1:length(organism_name_trimmed)) {
  if (n[i]>2)  organism_name_trimmed[i] <- word(organism_name_trimmed[i],1,3)
  if (n[i]==2) organism_name_trimmed[i] <- word(organism_name_trimmed[i],1,2)
  if (n[i]==1) organism_name_trimmed[i] <- word(organism_name_trimmed[i],1,1)
}

organism_name_trimmed <- strsplit(organism_name_trimmed, " ")

max.len <- max(sapply(organism_name_trimmed, length))
organism_name_trimmed <- lapply(organism_name_trimmed, function(x) {c(x, rep(NA, 
  max.len - length(x)))})

temp_names <- c(unlist(organism_name_trimmed))
organism_name_trimmed <- as.data.table(matrix(temp_names,
  ncol=3,nrow=length(temp_names)/3,byrow=T))
organism_name_trimmed[,entry_number := .I]

# SAVE
write.csv(organism_name_trimmed, paste0(temp_data_folder, 
  "organism_name_trimmed_temp_", "1", ".csv"), row.names=F) 

# generate table of unique organism names 
organism_name_trimmed_unique <- unique(organism_name_trimmed, 
  by=c("V1", "V2", "V3"))

# SAVE
write.csv(organism_name_trimmed_unique, paste0(temp_data_folder, 
  "organism_name_trimmed_unique_temp_", "1", ".csv"), row.names=F) 

# convert word table (unique) to long word list
organism_name_trimmed_cut <- as.data.table(gather(organism_name_trimmed_unique,
  word_number, word, V1,V2,V3))
organism_name_trimmed_cut[!is.na(word), word := gsub("[^[:alnum:]]","",word)]

# look up the individual words in the ncbi database using the taxise R package
genus <- tryCatch(tax_name(query=organism_name_trimmed_cut$word,get = tax_level, 
  db = "ncbi",ask=F))[[3]]

# SAVE
write.csv(genus, paste0(temp_data_folder, "genus_temp_", "1", ".csv"), 
  row.names=F) 

# merge genus results with organism words
organism_name_trimmed_cut$genus <- genus
organism_name_trimmed_cut[, number_of_classifications:=sum(!is.na(genus)), 
  by=c("entry_number")]
organism_name_trimmed_cut[is.na(genus), genus:=""]
organism_name_trimmed_cut[, ':='(word_number=NULL,word=NULL)]
organism_name_trimmed_cut[,genus_final:=gsub("^\\s+|\\s+$","",paste(genus,
  collapse=" ")),by=c("entry_number")]
organism_name_trimmed_cut <- unique(organism_name_trimmed_cut[,genus:=NULL], 
  by=c("entry_number"))
organism_name_trimmed_cut[number_of_classifications>1 & gsub("^(.*?) (.*)",  
  "\\1",genus_final) !=gsub("^(.* )(.*)",  "\\2",genus_final), genus_final:=""]
organism_name_trimmed_cut[number_of_classifications>1 & gsub("^(.*?) (.*)",  
  "\\1",genus_final) ==gsub("^(.* )(.*)",  "\\2",genus_final), 
  genus_final:=gsub("^(.*?) (.*)", "\\1",genus_final)]

organism_name_trimmed_unique <- organism_name_trimmed_unique[
  organism_name_trimmed_cut, .(V1,V2,V3,genus_final, entry_number), 
  on=c("entry_number")]

# extract entries with missing information
genus_missing_info <- organism_name_trimmed_unique[genus_final==""]

# missing genus information - exception#1 - account for plurals (most common)
plurals <- c("cocci", "bacilli","bacteria")
singulars <- c("coccus", "bacillus","bacterium")
var_words <- c("V1","V2","V3")

for (var in var_words) {
  mapply(function(datatable_variable,var_words_variable,pattern_old_variable, 
    pattern_new_variable) pattern_replacement(datatable_variable,var_words_variable,
    pattern_old_variable, pattern_new_variable),var_words_variable=var,
    pattern_new_variable=singulars, pattern_old_variable=plurals,
    MoreArgs = list(datatable_variable=organism_name_trimmed_unique))
}

# missing genus infomation - main exceptions - (i) gram and rods (correct - impose 
# some grouping upon range of gram variable types and eliminate words such as short) 
# and (ii) acid - acid fast bacilli
organism_genus_exceptions_1 <- "gram|rods"
organism_genus_exceptions_1_old <- c("suspected",",","short","(neg\\s)|(neg$)", 
  "(pos\\s)|(pos$)","non(\\s)*fermenting")
organism_genus_exceptions_1_new <-c("", "", "" ,"negative", "positive",
  "non-fermenting")
organism_genus_exceptions_2 <- "acid\\s|acid$"

organism_name_trimmed_unique[genus_final=="" & paste(V1,V2,V3,sep=" ") %like% 
  organism_genus_exceptions_1, genus_final := sprintf("non-genus: %s %s %s",
  V1,V2,V3)]
temp <- organism_name_trimmed_unique[paste(V1,V2,V3,sep=" ") %like% 
    organism_genus_exceptions_1]
mapply(function(datatable_variable,var_variable,pattern_old_variable, 
    pattern_new_variable) pattern_replacement(datatable_variable,var_variable,
    pattern_old_variable, pattern_new_variable),var_variable="genus_final",
    pattern_new_variable=organism_genus_exceptions_1_new, 
    pattern_old_variable=organism_genus_exceptions_1_old,
    MoreArgs = list(datatable_variable=temp))
organism_name_trimmed_unique[paste(V1,V2,V3,sep=" ") %like% 
    organism_genus_exceptions_1, genus_final:=temp$genus_final]

organism_name_trimmed_unique[genus_final=="" & paste(V1,V2,V3,sep=" ") %like% 
  organism_genus_exceptions_2, genus_final := "non-genus: acid fast bacilli"]

# re-run - account for plurals and exceptions
organism_name_trimmed_cut_temp <- as.data.table(gather(
  organism_name_trimmed_unique[genus_final==""], 
  word_number, word, V1,V2,V3))
organism_name_trimmed_cut_temp[!is.na(word), word := gsub("[^[:alnum:]]","",
  word)]

genus_temp_1 <- tryCatch(tax_name(query=organism_name_trimmed_cut_temp$word,
  get = tax_level, db = "ncbi",ask=F))[[3]]

# SAVE
write.csv(genus_temp_1, paste0(temp_data_folder, "genus_temp_1_temp_", "1", 
  ".csv"), row.names=F)

# merge genus _temp_1 results with organism words
organism_name_trimmed_cut_temp$genus <- genus_temp_1
organism_name_trimmed_cut_temp[, number_of_classifications:=sum(!is.na(genus)), 
  by=c("entry_number")]
organism_name_trimmed_cut_temp[is.na(genus), genus:=""]
organism_name_trimmed_cut_temp[, ':='(word_number=NULL,word=NULL)]
organism_name_trimmed_cut_temp[,genus_final:=gsub("^\\s+|\\s+$","",paste(genus,
  collapse=" ")),by=c("entry_number")]
organism_name_trimmed_cut_temp <- unique(organism_name_trimmed_cut_temp[,
  genus:=NULL], by=c("entry_number"))[order(entry_number)]
organism_name_trimmed_cut_temp[number_of_classifications>1 & gsub("^(.*?) (.*)",  
  "\\1",genus_final) !=gsub("^(.* )(.*)",  "\\2",genus_final), genus_final:=""]
organism_name_trimmed_cut_temp[number_of_classifications>1 & gsub("^(.*?) (.*)",  
  "\\1",genus_final) ==gsub("^(.* )(.*)",  "\\2",genus_final), genus_final:=
  gsub("^(.*?) (.*)",  "\\1",genus_final)]

organism_name_trimmed_unique <- organism_name_trimmed_unique[order(entry_number)]
organism_name_trimmed_unique[genus_final=="",genus_final :=
  organism_name_trimmed_cut_temp$genus_final ]

# re-run - account for still missing observations by considering multimatch 
# results (force taxize package to report back first result ) - report results as 
# 'uncertain'
organism_name_trimmed_cut_temp <- as.data.table(gather(organism_name_trimmed_unique[
  genus_final==""], word_number, word, V1,V2,V3))
organism_name_trimmed_cut_temp[!is.na(word), word := gsub("[^[:alnum:]]","",word)]

genus_temp_2 <- tryCatch(tax_name(query=organism_name_trimmed_cut_temp$word,get = 
	tax_level, db = "ncbi",rows=1))[[3]]

# SAVE
write.csv(genus_temp_2, paste0(temp_data_folder, "genus_temp_2_temp_", "1", 
  ".csv"), row.names=F)

# merge genus _temp_2 results with organism words
organism_name_trimmed_cut_temp$genus <- genus_temp_2
organism_name_trimmed_cut_temp[, number_of_classifications:=sum(!is.na(genus)), 
  by=c("entry_number")]
organism_name_trimmed_cut_temp[is.na(genus), genus:=""]
organism_name_trimmed_cut_temp[, ':='(word_number=NULL,word=NULL)]
organism_name_trimmed_cut_temp[,genus_final:=gsub("^\\s+|\\s+$","",paste(genus,
  collapse=" ")),by=c("entry_number")]
organism_name_trimmed_cut_temp <- unique(organism_name_trimmed_cut_temp[,
  genus:=NULL], by=c("entry_number"))[order(entry_number)]
organism_name_trimmed_cut_temp[number_of_classifications>1 & gsub("^(.*?) (.*)",  
  "\\1",genus_final) !=gsub("^(.* )(.*)",  "\\2",genus_final), genus_final:=""]
organism_name_trimmed_cut_temp[number_of_classifications>1 & gsub("^(.*?) (.*)",  
  "\\1",genus_final) ==gsub("^(.* )(.*)",  "\\2",genus_final), genus_final:=
  gsub("^(.*?) (.*)",  "\\1", genus_final)]

organism_name_trimmed_unique <- organism_name_trimmed_unique[order(entry_number)]
organism_name_trimmed_unique[genus_final=="",genus_final :=
  paste("uncertain result:", organism_name_trimmed_cut_temp$genus_final, sep=" ")][
  genus_final=="uncertain result: ", genus_final:=""]

# eliminate non-sensical categories (hard-coded exceptions)
organism_name_trimmed_unique[genus_final %in% c("This"),genus_final := ""]

# mark missing results as non-classified  
organism_name_trimmed_unique[genus_final=="",genus_final :=
  "non-classifiable organism"]

# merge results with mic_inf - step by step 
organism_name_trimmed_unique[genus_final=="",genus_final :=
  "non-classifiable organism"]

# SAVE
write.csv(organism_name_trimmed_unique, paste0(temp_data_folder, 
  "organism_name_trimmed_unique_temp_", "2", ".csv"), row.names=F) 

organism_name_trimmed_unique[, ':='(V1=unique(organism_name_trimmed, by=c("V1", 
  "V2", "V3"))$V1, V2=unique(organism_name_trimmed, by=c("V1", "V2", "V3"))$V2, 
  V3=unique(organism_name_trimmed, by=c("V1", "V2", "V3"))$V3)]

organism_name_trimmed <- organism_name_trimmed[organism_name_trimmed_unique,
  .(V1,V2,V3,genus_final, entry_number),on=c("V1", "V2", "V3")][
  order(entry_number)]

# SAVE
write.csv(organism_name_trimmed, paste0(temp_data_folder, 
  "organism_name_trimmed_temp_", "2", ".csv"), row.names=F) 

mic_inf[result=="pos_tier_1" & specimen_type %like% "blood", 
  result_mod_tier_pos:=organism_name_trimmed$genus_final]

# SAVE
write.csv(mic_inf, paste0(temp_data_folder, "mic_inf_temp_", "7", ".csv"),
	row.names=F) 

################################################################################
######### (E) IMPOSE THE FINAL CATEGORISATION (GRANULAR & AGGREGATED) ##########
################################################################################
################################################################################
# (a) aggregated categories (2)
mic_inf[result %in% c("invalid/ambigious_a", "invalid/ambigious_b",
  "pos_other", "posneg","neg","pos_tier_2"), result_mod := 
  "invalid/ambigious/neg/no result"][result %in% c("pos_tier_1"), result_mod := 
  "positive"]

# (b) granular categories (invalid/ambigious, positive tier 1, positive tier 2 
# and negative)
mic_inf[result %in% c("invalid/ambigious_a", "invalid/ambigious_b", "posneg"), 
  result_mod_granular := "invalid/ambigious"][result %in% c("pos_tier_1"),
  result_mod_granular := "positive_tier_1"][result %in% c("pos_tier_2"),
  result_mod_granular := "positive_tier_2"][result %in% c("pos_other","neg"), 
  result_mod_granular := "negative"]

# SAVE
write.csv(mic_inf, paste0(modified_data_folder, "mic_inf.csv"), row.names=F) 

################################################################################
#################  (F) GRAPHIC/TABLE OUTPUT - PROVIDE OVERVIEW #################
################################################################################
################################################################################
# (a) raw tables of results and their classification 
mic_inf_results <- unique(mic_inf[,list(specimen_type,
  test_name, test_comments,organism_name, organism_comment,result,result_mod, 
  result_mod_granular)])

# SAVE
write.csv(mic_inf_results, paste0(modified_data_folder, "mic_inf_results.csv"), 
	row.names=F) 

# (b) frequency/proportion tables
mic_inf_results_table <- addmargins(table(mic_inf$result_mod_granular))

# SAVE
write.csv(mic_inf_results_table, paste0(output_data_folder, 
	"mic_inf_results_table.csv"), row.names=F) 

mic_inf_results_prop_table <- round((prop.table(table(mic_inf$result_mod_granular))
  )*100, digits=2)

# SAVE
write.csv(mic_inf_results_prop_table, paste0(output_data_folder, 
	"mic_inf_results_prop_table.csv"), row.names=F) 

################################################################################
################################################################################
############################      END      #####################################
################################################################################
################################################################################