/* ----------------------------------------------------------------------------------------------*/
/* Longitudinal model of eGFR */

/* Analysis with all whites transplanted between January 1, 1990 and December 31, 2014*/

/* Read-in Cohort Data */ 

proc import file = "C:\Users\bstvock\Documents\research\surgery_transplant\Matas, Arthur\\egfr_cohort_all_post_donation_v3.txt" 
out = kidney dbms = tab replace;
run;


data kidney; set kidney;
kidney_cohort = 1;
run;

/* Read-in Longitudinal Data */ 

libname me "C:\Users\bstvock\Documents\research\surgery_transplant\Matas, Arthur\\";

data egfr; set me.req1266_lkd_egfr_2016jul22update;
run;

data egfr; set egfr;
drop tx_dt gender eGFR_60 eGFR_45 eGFR_30;
egfr_cohort = 1;
run;

/* Merge Cohort and Longitudinal Data */  
proc sort data = egfr;
by DONOR_ID; run;

proc sort data = kidney;
by DONOR_ID; run;

data egfr_total;
merge egfr kidney;
by DONOR_ID; run;

data egfr_total;
set egfr_total;
where kidney_cohort = 1 & egfr_cohort = 1; run;

data egfr_total; set egfr_total; 
*where yrs_tx_to_test > 1 & test_dt < '1NOV2015'd; 
where yrs_tx_to_test > 0.1 & test_dt le '30JUN2016'd; 
run;

data egfr_total; set egfr_total;
tx_year_2000 = max(tx_year - 2000, 0);
yrs_tx_to_test_1 = max(yrs_tx_to_test - 1, 0); 
yrs_tx_to_test_5 = max(yrs_tx_to_test - 5, 0);
yrs_tx_to_test_10 = max(yrs_tx_to_test - 10, 0);
yrs_tx_to_test_15 = max(yrs_tx_to_test - 15, 0);
yrs_tx_to_test_20 = max(yrs_tx_to_test - 20, 0);
age_at_donation_35 = max(age_at_donation - 35, 0); 
first_deg_primary = "No Hx"; 
if donation_eGFR = . then donation_eGFR = 87.89;
if donation_BMI = . then donation_BMI = 26;
if smoke_donation = . then smoke_donation = 0; 
if rel_category = 1 | FAM_KID_DZ_SW = "Y" then  first_deg_primary = "Primary Hx";
age_70_i = 0;
if age > 70 then age_70_i = (age - max(age_at_donation, 70));
age_60_i = 0;
if age > 60 then age_60_i = (age - max(age_at_donation, 60));
age_50_i = 0;
if age > 50 then age_50_i = (age - max(age_at_donation, 50));
age_at_donation_per_decade = age_at_donation/10;
age_at_donation_35_per_decade = age_at_donation_35/10;
donation_eGFR_per_10 = donation_eGFR/10;
donation_epi_per_10 = donation_epi/10;
donation_BMI_per_5 = donation_BMI/5;
donation_epi = 0;
if gender = 'F' & (0 < donation_cr le 0.7) then donation_epi = 144*((donation_cr/0.7)**(-0.329))*((0.993)**age_at_donation);
if gender = 'F' & donation_cr gt 0.7 then donation_epi = 144*((donation_cr/0.7)**(-1.209))*((0.993)**age_at_donation);
if gender = 'M' & (0 < donation_cr le 0.9) then donation_epi = 141*((donation_cr/0.9)**(-0.411))*((0.993)**age_at_donation);
if gender = 'M' & donation_cr gt 0.9 then donation_epi = 141*((donation_cr/0.9)**(-1.209))*((0.993)**age_at_donation);
run;


title "Parsimonious Model:  After 6 weeks"; 
proc mixed data = egfr_total method = ml;
class DONOR_ID first_deg_primary gender; 
model epi = age_at_donation age_at_donation_35 gender donation_epi
donation_BMI smoke_donation 
first_deg_primary tx_year  
yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20 
age_70_i  
first_deg_primary*yrs_tx_to_test 
/ s;
random intercept yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20/ type = un subject = DONOR_ID; 
ods output covparms=cov_model_final;
run;

title "Parsimonious Model:  After 6 weeks"; 
proc mixed data = egfr_total method = ml;
class DONOR_ID first_deg_primary gender; 
model epi = age_at_donation_per_decade age_at_donation_35_per_decade gender 
donation_epi_per_10
donation_BMI_per_5 
smoke_donation 
first_deg_primary 
tx_year 
yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20 
age_70_i
first_deg_primary*yrs_tx_to_test 
/ s;
random intercept yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20/ type = un subject = DONOR_ID; 
parms / parmsdata=cov_model_final;
estimate "slope 6 w -5 yrs post-donation fhx <70 years" yrs_tx_to_test 1 first_deg_primary*yrs_tx_to_test 0 1 ;
estimate "slope 5-10 yrs post-donation fhx <70 years" yrs_tx_to_test 1  yrs_tx_to_test_5 1 first_deg_primary*yrs_tx_to_test 0 1 ; 
estimate "slope 10-20 yrs post-donation fhx <70 years" yrs_tx_to_test 1 yrs_tx_to_test_5 1 yrs_tx_to_test_10 1 first_deg_primary*yrs_tx_to_test 0 1 ; 
estimate "slope 20+ yrs post-donation fhx <70 years" yrs_tx_to_test 1  yrs_tx_to_test_5 1 yrs_tx_to_test_10 1 yrs_tx_to_test_20 1 first_deg_primary*yrs_tx_to_test 0 1 ;
estimate "age at donation >35 per decade" age_at_donation_per_decade 1 age_at_donation_35_per_decade 1;
estimate "age 55 versus 25" age_at_donation_per_decade 3 age_at_donation_35_per_decade 2;
run;

/* Analysis with all whites */

proc import file = "C:\Users\bstvock\Documents\research\surgery_transplant\Matas, Arthur\\egfr_cohort_all_post_donation_all_whites_v3.txt" 
out = kidney_all dbms = tab replace;
run;

data kidney_all; set kidney_all;
kidney_cohort = 1;
run;

/* Merge Cohort and Longitudinal Data */  
proc sort data = egfr;
by DONOR_ID; run;

proc sort data = kidney_all;
by DONOR_ID; run;

data egfr_total_all;
merge egfr kidney_all;
by DONOR_ID; run;

data egfr_total_all;
set egfr_total_all;
where kidney_cohort = 1 & egfr_cohort = 1; run;

data egfr_total_all; set egfr_total_all; 
where yrs_tx_to_test > 0.1 & test_dt le '30JUN2016'd; 
run;

data egfr_total_all; set egfr_total_all;
tx_year_2000 = max(tx_year - 2000, 0);
yrs_tx_to_test_1 = max(yrs_tx_to_test - 1, 0); 
yrs_tx_to_test_5 = max(yrs_tx_to_test - 5, 0);
yrs_tx_to_test_10 = max(yrs_tx_to_test - 10, 0);
yrs_tx_to_test_15 = max(yrs_tx_to_test - 15, 0);
yrs_tx_to_test_20 = max(yrs_tx_to_test - 20, 0);
age_at_donation_35 = max(age_at_donation - 35, 0); 
first_deg_primary = "No Hx"; 
if donation_eGFR = . then donation_eGFR = 87.89;
if donation_BMI = . then donation_BMI = 26;
if smoke_donation = . then smoke_donation = 0; 
if rel_category = 1 | FAM_KID_DZ_SW = "Y" then  first_deg_primary = "Primary Hx";
age_70_i = 0;
if age > 70 then age_70_i = (age - max(age_at_donation, 70));
age_60_i = 0;
if age > 60 then age_60_i = (age - max(age_at_donation, 60));
age_50_i = 0;
if age > 50 then age_50_i = (age - max(age_at_donation, 50));
age_at_donation_per_decade = age_at_donation/10;
age_at_donation_35_per_decade = age_at_donation_35/10;
donation_eGFR_per_10 = donation_eGFR/10;
donation_BMI_per_5 = donation_BMI/5;
donation_epi = 0;
if gender = 'F' & (0 < donation_cr le 0.7) then donation_epi = 144*((donation_cr/0.7)**(-0.329))*((0.993)**age_at_donation);
if gender = 'F' & donation_cr gt 0.7 then donation_epi = 144*((donation_cr/0.7)**(-1.209))*((0.993)**age_at_donation);
if gender = 'M' & (0 < donation_cr le 0.9) then donation_epi = 141*((donation_cr/0.9)**(-0.411))*((0.993)**age_at_donation);
if gender = 'M' & donation_cr gt 0.9 then donation_epi = 141*((donation_cr/0.9)**(-1.209))*((0.993)**age_at_donation);
donation_epi_per_10 = donation_epi/10;
run;


title "Parsimonious Model:  After 6 weeks; all whites"; 
proc mixed data = egfr_total_all method = ml;
class DONOR_ID first_deg_primary gender; 
model epi = age_at_donation age_at_donation_35 gender donation_epi
donation_BMI smoke_donation 
first_deg_primary tx_year  
yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20 
age_70_i  
first_deg_primary*yrs_tx_to_test 
/ s;
random intercept yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20/ type = un subject = DONOR_ID; 
parms / parmsdata=cov_model_final;
ods output covparms=cov_model_final_all;
run;


title "Parsimonious Model:  After 6 weeks; all whites"; 
proc mixed data = egfr_total_all method = ml;
class DONOR_ID first_deg_primary gender; 
model epi = age_at_donation_per_decade age_at_donation_35_per_decade gender 
donation_epi_per_10
donation_BMI_per_5 
smoke_donation 
first_deg_primary 
tx_year 
yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20 
age_70_i
first_deg_primary*yrs_tx_to_test 
/ s;
random intercept yrs_tx_to_test yrs_tx_to_test_5 yrs_tx_to_test_10 yrs_tx_to_test_20/ type = un subject = DONOR_ID; 
parms / parmsdata=cov_model_final_all;
estimate "slope 6 w -5 yrs post-donation fhx <70 years" yrs_tx_to_test 1 first_deg_primary*yrs_tx_to_test 0 1 ;
estimate "slope 5-10 yrs post-donation fhx <70 years" yrs_tx_to_test 1  yrs_tx_to_test_5 1 first_deg_primary*yrs_tx_to_test 0 1 ; 
estimate "slope 10-20 yrs post-donation fhx <70 years" yrs_tx_to_test 1 yrs_tx_to_test_5 1 yrs_tx_to_test_10 1 first_deg_primary*yrs_tx_to_test 0 1 ; 
estimate "slope 20+ yrs post-donation fhx <70 years" yrs_tx_to_test 1  yrs_tx_to_test_5 1 yrs_tx_to_test_10 1 yrs_tx_to_test_20 1 first_deg_primary*yrs_tx_to_test 0 1 ;
estimate "age at donation >35 per decade" age_at_donation_per_decade 1 age_at_donation_35_per_decade 1;
estimate "age 55 versus 25" age_at_donation_per_decade 3 age_at_donation_35_per_decade 2;
run;
