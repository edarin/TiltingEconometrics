/***************************************************************************************************/
/* AST NLSY79 dataset preparation 	 					   		         						   */
/* Bryan S. Graham, NYU (w/ Dan Egel and Cristine Pinto)			         			   		   */
/* bsg1@nyu.edu                 						         		   			   			   */
/* March 2011                               								      				   */
/***************************************************************************************************/

/***************************************************************************************************/
/* This do file and the accompanying Stata dictionary file report estimation results and figures   */
/* presented in the paper "Auxiliary-to-Study Tilting". The data and do                            */
/* file are provided "as is". I am unable to assist with their interpretation or use. However      */
/* please do feel free to e-mail me if you find any mistakes at bryan.graham@nyu.edu.              */
/***************************************************************************************************/

/* use a semicolon as the command delimiter */
#delimit ;

clear matrix;
clear;

set matsize 800;
set memory 100m;

/* Adjust the SOURCE_DATA directory to point to the location of the NLSY_BlkWhtGap.DCT dictionary file. Adjust the    */
/* WRITE_DATA and DO_FILES directorys to point to the location of where you would like to write any created files and */
/* where you have placed this do file respectively. */

global SOURCE_DATA "C:\Documents and Settings\bsg1\My Documents\BSG_WORK_19W4th\Research\AST_11Spr\EmpiricalApplication\Source_Data";
global WRITE_DATA "C:\Documents and Settings\bsg1\My Documents\BSG_WORK_19W4th\Research\AST_11Spr\EmpiricalApplication\Created_Data";
global DO_FILES "C:\Documents and Settings\bsg1\My Documents\BSG_WORK_19W4th\Research\AST_11Spr\EmpiricalApplication\Stata_Do";

/* read in source data (extract from April 30, 2008 release for NLSY79) */
infile using "$SOURCE_DATA\NLSY_BlkWhtGap.DCT";

g HHID_79 = R0000149;	/* household ID number (for `clustering') */						

/* parents years of completed schooling at baseline */
g DadSch_in_79r = R0007900 if R0007900>=0;
g MomSch_in_79r = R0006500 if  R0006500>=0;

/* Basic respondent demographics */
g usborn = (R0000700==1);
g mother_usborn = (R0006100==1);
g father_usborn = (R0007300==1);
g male = (R0214800==1);
g hispanic = (R0214700==1);
g black = (R0214700==2);
g born1962to1964 = (R0000500>=62);
g yearborn = R0000500;
g yearborn62 = (R0000500==62);
g yearborn63 = (R0000500==63);
g yearborn64 = (R0000500==64);

/* sample weights */
g core_sample = (R0173600<=8 | R0173600==10  | R0173600==11 | R0173600==13  | R0173600==14);
g male_blkwhthis_sample = (R0173600<=4 | R0173600==10 | R0173600==11);
g sample_wgts =  R0216100;

/* age in base survey year */
g AgeIn1979 = R0000600;
g Age13In1979 = (AgeIn1979==13);
g Age14In1979 = (AgeIn1979==14);
g Age15In1979 = (AgeIn1979==15);
g Age16In1979 = (AgeIn1979==16);
g Age17In1979 = (AgeIn1979==17);
g Age18In1979 = (AgeIn1979==18);
g Age19In1979 = (AgeIn1979==19);
g Age20In1979 = (AgeIn1979==20);
g Age21In1979 = (AgeIn1979==21);
g Age22In1979 = (AgeIn1979==22);

/* AFQT percentile */
g AFQT = R0618300 if R0618300>0;
g AFQT_NoProb = (R0614800==51);			/* AFQT score based on test with no reported "problems" */
g AFQT_Adj1 = AFQT if AFQT_NoProb==1;   /* AFQT scores, problem free only */
g AFQT_Adj2 = invnormal(AFQT_Adj1/100) if AFQT_Adj1~=.;	   /* transform to approximate normality */

/* Calculate real annual earnings 1990 to 1993 (1993 prices) */
/* CPI with 1982-84 = 100: 1990: 130.7, 1991: 136.2, 1992: 140.3, 1993: 144.5 */

g earnings90 = R3559001*(144.5/130.7) if R3559001>=0;
g earnings91 = R3897101*(144.5/136.2) if R3897101>=0;
g earnings92 = R4295101*(144.5/140.3) if R4295101>=0;
g earnings93 = R4982801               if R4982801>=0;

egen AvgEarnings_90to93 = rowmean(earnings90 earnings91 earnings92 earnings93);
g LogEarn = log(AvgEarnings_90to93);

/* Calculate average hourly wages */
g wages90 = R3127800*(144.5/130.7) if R3127800>=100 & R3127800<=7500;
g wages91 = R3523500*(144.5/136.2) if R3523500>=100 & R3523500<=7500;
g wages92 = R3728500*(144.5/140.3) if R3728500>=100 & R3728500<=7500;
g wages93 = R4416800 			   if R4416800>=100 & R4416800<=7500;

egen AvgHourlyWages_90to93 = rowmean(wages90 wages91 wages92 wages93);
g LogWage = log(AvgHourlyWages_90to93);

/* replication of Neal and Johnson (1996) sample (without hispanics) */
g NJ_target_sample = (male_blkwhthis_sample==1 & born1962to1964==1 & hispanic ~= 1 & hispanic ~= .);
g NJ_sample = (male_blkwhthis_sample==1 & born1962to1964==1 & hispanic ~= 1 & LogWage~=. & yearborn~=. & black ~= .  & hispanic ~= . & AFQT_Adj1~=.);
tab NJ_target_sample;
tab NJ_sample;

log using "$WRITE_DATA\NLSY79_Sample_Log", replace;
log on;
table black [pweight=sample_wgts], c(mean yearborn mean AFQT_Adj1 mean LogWage);

/* summary statistics on black and white differences */
reg AvgHourlyWages_90to93 black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
reg LogWage black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
reg AFQT_Adj1 black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
reg yearborn62 black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
reg yearborn63 black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
reg yearborn64 black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);

/* replicate Neal and Johnson (1996) basic finding */
reg LogWage yearborn62-yearborn64 black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79) nocons;
reg LogWage yearborn62-yearborn64 black AFQT_Adj1 [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79) nocons;
reg LogWage yearborn62-yearborn64 black AFQT_Adj2 [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79) nocons;

/* compute estimates of average log hourly wages by race */
reg LogWage       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg LogWage       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg LogWage black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);

/* compute estimates of CDFs of black and white wage distributions (and differences at various points) */
g t = (AvgHourlyWages_90to93<=500);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg t black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
replace t = (AvgHourlyWages_90to93<=750);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg t black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
replace t = (AvgHourlyWages_90to93<=1000);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg t black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
replace t = (AvgHourlyWages_90to93<=1250);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg t black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
replace t = (AvgHourlyWages_90to93<=1500);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg t black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
replace t = (AvgHourlyWages_90to93<=1750);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg t black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);
replace t = (AvgHourlyWages_90to93<=2000);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==1, cluster(HHID_79);
reg t       [pweight=sample_wgts] if NJ_sample==1 & black==0, cluster(HHID_79);
reg t black [pweight=sample_wgts] if NJ_sample==1, cluster(HHID_79);

log off;
log close;

keep if NJ_sample==1;
sort HHID_79;
outsheet sample_wgts HHID_79 AvgHourlyWages_90to93 LogWage yearborn black AFQT_Adj1 if NJ_sample==1 using "$WRITE_DATA\NLSY79_Sample.out", replace;
save "$WRITE_DATA\NLSY79_Sample.dta", replace;
