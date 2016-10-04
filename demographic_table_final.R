#-------------------------------------------------------------------------------------------------#
# Create demographic table by family history
#-------------------------------------------------------------------------------------------------#

### libraries used
library(dplyr)

### functions used
cont.summary <- function(variable, digits) {
	cont.median <- tapply(variable,kidney$first_deg_primary,median, na.rm = TRUE)
	cont.25 <- tapply(variable,kidney$first_deg_primary,quantile,0.25, na.rm = TRUE)
	cont.75 <- tapply(variable,kidney$first_deg_primary,quantile,0.75, na.rm = TRUE)
	#p.value <- t.test(variable[hemo$contrast=="v"],variable[hemo$contrast=="i"])$p.value
	p.value.2 <- wilcox.test(variable ~ kidney$first_deg_primary)$p.value 
	p.value.2 <- round(p.value.2, 3)
	p.value.2 <- ifelse(p.value.2 == 0, "<0.001", p.value.2)
	cont.summary <- c(paste(round(cont.median[1], digits), " (", round(cont.25[1], digits), ", ", 
		round(cont.75[1], digits), ")", sep =""), paste(round(cont.median[2], digits),
			" (", round(cont.25[2], digits), ", ", 
			round(cont.75[2], digits), ")", sep ="") , p.value.2)
	return(cont.summary) }

cont.summary <- function(variable, byvar, digits) {
	cont.median <- tapply(variable, byvar, median, na.rm = TRUE)
	cont.25 <- tapply(variable,byvar,quantile,0.25, na.rm = TRUE)
	cont.75 <- tapply(variable,byvar,quantile,0.75, na.rm = TRUE)
	#p.value <- t.test(variable[hemo$contrast=="v"],variable[hemo$contrast=="i"])$p.value
	p.value.2 <- wilcox.test(variable ~ byvar)$p.value 
	p.value.2 <- round(p.value.2, 3)
	p.value.2 <- ifelse(p.value.2 == 0, "<0.001", p.value.2)
	cont.summary <- c(paste(round(cont.median[1], digits), " (", round(cont.25[1], digits), ", ", 
		round(cont.75[1], digits), ")", sep =""), paste(round(cont.median[2], digits),
			" (", round(cont.25[2], digits), ", ", 
			round(cont.75[2], digits), ")", sep ="") , p.value.2)
	return(cont.summary) }

binary.summary <- function(variable, byvar) {
	table.bin <- table(variable, byvar)
	prop.bin <- prop.table(table.bin,2)*100
	chisq.bin <- chisq.test(table.bin)$p.value
	chisq.exp <- chisq.test(table.bin)$expected
	if (min(chisq.exp)<=5) {
		chisq.bin <- fisher.test(table.bin)$p.value }
	chisq.bin <- round(chisq.bin, 3)
	chisq.bin <- ifelse(chisq.bin  == 0, "<0.001", chisq.bin)
	binary.summary <- c(paste(table.bin[2,1], " (", round(prop.bin[2,1], 1), ")", sep =""),
		paste(table.bin[2,2], " (", round(prop.bin[2,2], 1), ")", sep =""), chisq.bin)
	return(binary.summary) }

### set file direction
file.dir <- "C:\\Users\\bstvock\\Documents\\research\\surgery_transplant\\Matas, Arthur\\"
file.dir <- "C:\\Users\\David\\Google Drive\\research\\surgery_transplant\\Matas, Arthur\\"

### read-in data 
#kidney <- read.table(paste(file.dir, "egfr_cohort_v2.txt", sep=""), sep = "\t", header = T)
#kidney <- read.table(paste(file.dir, "egfr_cohort_all_whites_v2.txt", sep=""), sep = "\t", 
#	header = T)
kidney <- read.table(paste(file.dir, "egfr_cohort_all_post_donation_v3.txt", sep=""), sep = "\t", 
	header = T)
kidney <- read.table(paste(file.dir, "egfr_cohort_all_post_donation_all_whites_v3.txt", sep=""), sep = "\t", 
	header = T)

### define new variable which combines first degree relative with primary hx of kidney dx 
### for unrelated donors
table(kidney$rel_category, kidney$FAM_KID_DZ_SW)

kidney <- mutate(kidney, 
	first_deg_primary = ifelse(rel_category == 1 | FAM_KID_DZ_SW == "Y", 
		"Primary Hx", "No Hx"),
	female = (gender == "F")*1) 

demog_by_group <- cbind(
	cont.summary(kidney$age_at_donation, kidney$first_deg_primary, 1),
	binary.summary(kidney$female, kidney$first_deg_primary),
	cont.summary(kidney$donation_epi, kidney$first_deg_primary, 1),
	cont.summary(kidney$donation_cr, kidney$first_deg_primary, 2),
	cont.summary(kidney$donation_BMI, kidney$first_deg_primary, 1),
	cont.summary(kidney$donation_gluc, kidney$first_deg_primary, 0),
	cont.summary(kidney$donation_hgb, kidney$first_deg_primary, 1),
	cont.summary(kidney$SYS_donation, kidney$first_deg_primary, 0),
	binary.summary(kidney$smoke_donation, kidney$first_deg_primary),
	cont.summary(kidney$tx_year, kidney$first_deg_primary, 0),
	cont.summary(kidney$follow_up_time, kidney$first_deg_primary, 1),
	binary.summary(kidney$epi_N_one_year > 0, kidney$first_deg_primary),
	cont.summary(kidney$epi_N, kidney$first_deg_primary, 0),
	cont.summary(kidney$epi_N[kidney$epi_N_one_year > 0], 
		kidney$first_deg_primary[kidney$epi_N_one_year > 0], 0),
	cont.summary(kidney$epi_N[kidney$epi_N_one_year == 0], 
		kidney$first_deg_primary[kidney$epi_N_one_year == 0], 0)
	)

demog_by_group_df <- data.frame(No_Hx = demog_by_group[1, ], Primary_Hx = demog_by_group[2, ], 
	p_value = demog_by_group[3, ])

rownames(demog_by_group_df) <- c("Age", "Female", "eGFR", "Creatinine", "BMI", "Glucose", 
	 "Hemoglobin", "SBP", "Smoker", "Year of Donoation", "Follow-Up Length",
	"Follow-Up > 1 Year", 
	"Number of eGFR Measurements",
	"Number of eGFR Measurements Subjects > 1yr",
	"Number of eGFR Measurements Subjects < 1yr"
	)
print(demog_by_group_df)

