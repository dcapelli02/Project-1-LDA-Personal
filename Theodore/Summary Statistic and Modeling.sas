data alzheimer25;
	set '/home/u63619091/LDA HWs/HW1/alzheimer25.sas7bdat';
run;

/* reformat data to long */
data alzheimer_long;
    set alzheimer25;
    
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6] bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;
    
    do TIME = 0 to 6;
        CDRSB = cdrsb_arr[TIME];
        BPRS = bprs_arr[TIME];
        ABPET = abpet_arr[TIME];
        TAUPET = taupet_arr[TIME];
        
        if not missing(CDRSB) or not missing(BPRS) or 
           not missing(ABPET) or not missing(TAUPET) then output;
    end;

    keep TRIAL PATID SEX AGE EDU BMI INKOMEN JOB ADL WZC 
         TIME CDRSB BPRS ABPET TAUPET;
run;

/* long format with only baseline measurements for covariates that are only using baseline */
data alzheimer_long2;
    set alzheimer25;
    
    array bprs_arr[0:6] bprs0-bprs6;   /* only BPRS varies over time */
    
    /* baseline variables */
    CDRSB0 = cdrsb0;
    ABPET0 = abpet0;
    TAUPET0 = taupet0;
    
    do TIME = 0 to 6;
        BPRS = bprs_arr[TIME];
        
       
        if not missing(BPRS) then do;
            CDRSB = CDRSB0;   /* keep baseline */
            ABPET = ABPET0;
            TAUPET = TAUPET0;
            output;
        end;
    end;

    keep TRIAL PATID SEX AGE EDU BMI INKOMEN JOB ADL WZC 
         TIME BPRS CDRSB ABPET TAUPET;
run;

/* Summary Statistic: AUC */

/* Extract baseline BPRS (time = 0) */
data baseline;
    set alzheimer_long;
    where time = 0;
    keep patid bprs;
    rename bprs = baseline_bprs;
run;

/* Merge baseline BPRS into the AUC dataset */
proc sort data=baseline; by patid; run;
proc sort data=auc; by patid; run;

data auc;
    merge auc(in=a) baseline(in=b);
    by patid;
    if a; /* keep only subjects with AUC */
run;

proc glm data=auc;
    model auc = baseline_bprs;
run;


proc mixed data=alzheimer_long method=reml;
	class patid age sex edu bprs time;
	model bprs = age sex edu bprs time;
	repeated age / type=arh(1) subject = patid;
run;


proc mixed data=alzheimer_long2 method=ml;
  class patid time;
  model bprs = trial wzc time age cdrsb abpet inkomen sex age*time sex*time / solution;
  random intercept / subject=patid g gcorr v vcorr;
  repeated time / type=toep subject=patid r rcorr;
run;

/* Full Multivariate model (random intercept, no random effects) */
proc mixed data=alzheimer_long2 method=reml;
	class patid time wzc job edu;
	model bprs = time age job bmi adl wzc cdrsb abpet edu inkomen sex taupet age*time sex*time 
  				job*time bmi*time wzc*time cdrsb*time abpet*time inkomen*time adl*time/ solution;
  random intercept / type=un subject=patid g gcorr v vcorr;
  
/* Reduced Multivariate model - removing non-significant terms*/
proc mixed data=alzheimer_long2 method=reml;
	class patid time wzc job;
	model bprs = time age job bmi adl wzc cdrsb cdrsb*time /solution;
	random intercept / type =un subject=patid g gcorr v vcorr;
run;

/* Two-Stage Analysis attempt - idk man */
proc sort data=alzheimer_analysis; by PATID; run;

ods output ParameterEstimates=stage1_estimates;
proc reg data=alzheimer_analysis noprint;
    by PATID;
    model BPRS = TIME;
run;
quit;

/* Keep only slope estimates */
data slopes;
    set stage1_estimates;
    where Variable = "TIME";
    keep PATID Estimate StdErr;
    rename Estimate = Slope StdErr = Slope_SE;
run;

proc glm data=slopes;
    class SEX TRIAL;
    model Slope = AGE SEX TRIAL CDRSB0 ABPET0 TAUPET0;
    title "Stage 2: Predicting Individual BPRS Slopes from Baseline Covariates";
run;
quit;


/* Full Random effects model*/
proc mixed data=alzheimer_long2 method=ml;
  class patid sex job wzc;
  model bprs = time age job bmi adl wzc cdrsb abpet inkomen sex taupet age*time sex*time 
  				job*time bmi*time wzc*time cdrsb*time abpet*time inkomen*time adl*time/ solution;
  random intercept time / type=un subject=patid g gcorr v vcorr;
run;

/* Remove covariates with non-significant p-values */
/* removed ABPET and inkomen, and all interaction terms EXCEPT time with CDRSB and with adl */
/* time*adl is only slightly non-significant */

proc mixed data=alzheimer_long2 method=reml;
	class patid sex job wzc;
	model bprs = time age job bmi adl wzc cdrsb sex time*cdrsb time*adl /solution;
	random intercept time / type=un subject=patid g gcorr v vcorr;
run;


