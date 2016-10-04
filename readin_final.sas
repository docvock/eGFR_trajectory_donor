/* ----------------------------------------------------------------------------------------------*/
/* Create dataset of baseline information on white donors from January 1, 1990 - December 31, 2013
with at least one eGFR measurement one year post-donation */
/* Version 3 only examines measurement after 6 weeks */

 
OPTIONS nofmterr;
libname kidney "C:\Users\bstvock\Documents\research\surgery_transplant\Matas, Arthur"; 

%let all_white = "TRUE";

data kidney; set kidney.nov2015_sds_donor; run;

data kidney; set kidney;
rel_category = 0;
if donor_relationship = "Child" | donor_relationship = "Father" | donor_relationship = "Mother" |
	donor_relationship = "Full Sibling" |  donor_relationship = "Identical Twin" 
	then rel_category  = 1;
if donor_relationship = "Half Sibling" | donor_relationship = "Other Relative"  
	then rel_category = 2;
kidney_cohort = 1;
rename current_eGFR_dt = current_eGFR_dt_alt;
run;


/* Read-in longitudinal information to get date of last eGFR measurements and number of eGFR measurements */

data egfr; set kidney.req1266_lkd_egfr_2016jul22update; run;

proc sort data = egfr;
by DONOR_ID test_dt; run;

proc print data = egfr (obs = 20); run;

data egfr; set egfr;
by DONOR_ID; 
where test_dt le '30JUN2016'd;
last_egfr_measurement = last.DONOR_ID; run;

data egfr_current; set egfr;
where last_egfr_measurement  = 1;
rename test_dt = current_egfr_dt; run;

data egfr_current; set egfr_current(keep = DONOR_ID current_egfr_dt);
run;

proc tabulate data = egfr out = donor_obs;
where yrs_tx_to_test > 0.1;
  class DONOR_ID;
  var epi;
  tables DONOR_ID, epi*(N);
run;

proc tabulate data = egfr out = donor_obs_one_year;
where yrs_tx_to_test > 1;
  class DONOR_ID;
  var epi;
  tables DONOR_ID, epi*(N);
run;


data donor_obs; set donor_obs;
drop _TYPE_ _PAGE_ _TABLE_ ; run;

data donor_obs_one_year; set donor_obs_one_year;
drop _TYPE_ _PAGE_ _TABLE_ ; 
rename epi_N = epi_N_one_year;
run;


/* Merge Cohort and Longitudinal Count and Date Data */  
proc sort data = donor_obs;
by DONOR_ID; run;

proc sort data = donor_obs_one_year;
by DONOR_ID; run;

proc sort data = egfr_current;
by DONOR_ID; run;

proc sort data = kidney;
by DONOR_ID; run;

proc print data = kidney (obs = 20); run;

data kidney;
merge donor_obs donor_obs_one_year egfr_current kidney;
by DONOR_ID; run;

data kidney; set kidney;
if epi_N_one_year = . then epi_N_one_year = 0; 
where kidney_cohort = 1; run;
%macro dataset;

%if &all_white = "FALSE" %then %do; 
data kidney_short; set kidney(keep =  age_at_donation rel_category donation_eGFR donation_BMI 
	gender donation_gluc donation_cr donation_hgb SYS_donation smoke_donation tx_dt race 
	FAM_KID_DZ_SW donor_relationship donation_ACR current_eGFR_dt DONOR_ID epi_N epi_N_one_year ESRD);
where tx_dt le '31DEC2014'd & tx_dt ge '1JAN1990'd & race = "Caucasian / White" ; 
tx_time = (tx_dt - 1253)/365;
tx_year = year(tx_dt);
follow_up_time = (current_eGFR_dt - tx_dt)/365;
donation_epi = 0;
if gender = 'F' & (0 < donation_cr le 0.7) then donation_epi = 144*((donation_cr/0.7)**(-0.329))*((0.993)**age_at_donation);
if gender = 'F' & donation_cr gt 0.7 then donation_epi = 144*((donation_cr/0.7)**(-1.209))*((0.993)**age_at_donation);
if gender = 'M' & (0 < donation_cr le 0.9) then donation_epi = 141*((donation_cr/0.9)**(-0.411))*((0.993)**age_at_donation);
if gender = 'M' & donation_cr gt 0.9 then donation_epi = 141*((donation_cr/0.9)**(-1.209))*((0.993)**age_at_donation);
run;
%end; 

%if &all_white = "TRUE" %then %do; 
data kidney_short; set kidney(keep =  age_at_donation rel_category donation_eGFR donation_BMI 
	gender donation_gluc donation_cr donation_hgb SYS_donation smoke_donation tx_dt race 
	FAM_KID_DZ_SW donor_relationship donation_ACR current_eGFR_dt DONOR_ID epi_N epi_N_one_year ESRD);
where tx_dt le '31DEC2014'd & race = "Caucasian / White" ; 
tx_time = (tx_dt - 1253)/365;
tx_year = year(tx_dt);
follow_up_time = (current_eGFR_dt - tx_dt)/365;
donation_epi = 0;
if gender = 'F' & (0 < donation_cr le 0.7) then donation_epi = 144*((donation_cr/0.7)**(-0.329))*((0.993)**age_at_donation);
if gender = 'F' & donation_cr gt 0.7 then donation_epi = 144*((donation_cr/0.7)**(-1.209))*((0.993)**age_at_donation);
if gender = 'M' & (0 < donation_cr le 0.9) then donation_epi = 141*((donation_cr/0.9)**(-0.411))*((0.993)**age_at_donation);
if gender = 'M' & donation_cr gt 0.9 then donation_epi = 141*((donation_cr/0.9)**(-1.209))*((0.993)**age_at_donation);
run;
%end; 

data kidney_short; set kidney_short;
where current_eGFR_dt ne . & follow_up_time > 0.1;
kidney_cohort = 1;
drop race; run;

%if &all_white = "FALSE" %then %do; 
proc export outfile = 
	"C:\Users\bstvock\Documents\research\surgery_transplant\Matas, Arthur\egfr_cohort_all_post_donation_v3.txt" 
data = kidney_short dbms = tab replace;
run;
%end;

%if &all_white = "TRUE" %then %do; 
proc export outfile = 
	"C:\Users\bstvock\Documents\research\surgery_transplant\Matas, Arthur\egfr_cohort_all_post_donation_all_whites_v3.txt" 
data = kidney_short dbms = tab replace;
run;
%end;
%mend;

%dataset;
