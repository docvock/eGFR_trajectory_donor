#-----------------------------------------------------------------------------#
# Plot estimated post-donation eGFR trajectory by family history and age at donation based on 
# linear mixed-effect models fit in SAS
#-----------------------------------------------------------------------------#

## Coeficients from longitudinal model using data on white donors since 1990
cv <- c(-256.48,
	-0.4195,
	0.2120,
	1.4412,
	0.4722,
	-0.1332,
	1.4251,
	-0.4877,
	0.1460,
	1.1193,
	-0.8754,
	-0.1708,
	-0.2654,
	-0.3807,
	0.2002
)

## Coeficients from longitudinal model using data on all white donors
cv <- c(-126.30,
	-0.4994,
	0.2032,
	1.3829,
	0.3634,
	-0.1288,
	1.6413,
	-0.6790,
	0.08791,
	1.0476,
	-0.9376,
	-0.2048,
	-0.1959,
	-0.2657,
	0.2396
)


estimate_egfr <- function(tx_year, gender, age, yrs_tx, donation_eGFR, 
	donation_BMI, donation_smoke, donation_gluc, kdx_hx) {
	yrs_tx_1 <- pmax(yrs_tx - 1, 0)
	yrs_tx_5 <- pmax(yrs_tx - 5, 0)
	yrs_tx_10 <- pmax(yrs_tx - 10, 0)
	yrs_tx_20 <- pmax(yrs_tx - 20, 0)
	age_35 <- ifelse(age > 35, age-35 , 0)
	attained_age <- age+yrs_tx
	attained_age_70 <- ifelse(attained_age > 70, attained_age - pmax(age, 70),0)
	attained_age_50 <- ifelse(attained_age > 50, attained_age - pmax(age, 50),0)
	
	estimate <- colSums(rbind(1, age, age_35,  gender, 
		donation_eGFR, donation_BMI, donation_smoke,
		kdx_hx,
		tx_year,
		yrs_tx, yrs_tx_5, yrs_tx_10, yrs_tx_20,
		attained_age_70, 
		kdx_hx*yrs_tx)*cv)
	return(estimate)
}



tx_year_use <- 1995
est_25_hx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(25), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 0)
est_35_hx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(35), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 0)
est_45_hx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(45), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 0)
est_55_hx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(55), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 0)


est_25_nhx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(25), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 1)
est_35_nhx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(35), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 1)
est_45_nhx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(45), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 1)
est_55_nhx <- estimate_egfr(tx_year = tx_year_use, gender = 1, age = c(55), 
	yrs_tx = c(1:25), donation_eGFR = 100, donation_BMI = 25,
	donation_gluc = 90, donation_smoke = 0, kdx_hx = 1)

plot(26:50, est_25_nhx, col = "red", type = "l", ylim = c(min(c( 
	est_55_hx, est_55_nhx)), max(est_25_hx, est_25_nhx)), xlim = c(25, 80), lwd =2,
	xlab = "Attained Age (years)", ylab = "Mean eGFR (mL/min/1.73 m2)")
lines(26:50, est_25_hx, col = "red", lty = 2, lwd = 2)
lines(36:60, est_35_nhx, col = "orange", lwd = 2)
lines(36:60, est_35_hx, col = "orange", lty = 2, lwd = 2)
lines(46:70, est_45_nhx, col = "green", lwd = 2)
lines(46:70, est_45_hx, col = "green", lty = 2, lwd = 2)
lines(56:80, est_55_nhx, col = "blue", lwd = 2)
lines(56:80, est_55_hx, col = "blue", lty = 2, lwd = 2)

legend("topright", c("Age 25 No FHx", "Age 25 FHx",
	"Age 35 No FHx", "Age 35 FHx", "Age 45 No FHx", "Age 45 FHx",
	"Age 55 No FHx", "Age 55 FHx"), col = c("red", "red", "orange", "orange",
		"green", "green", "blue", "blue"), lty = c(1, 2, 1, 2, 1, 2, 1, 2))
