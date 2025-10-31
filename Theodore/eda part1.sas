data alzheimer_analysis;
    merge alzheimer_long(in=a) 
          alzheimer_baseline(keep=PATID CDRSB0 ABPET0 TAUPET0);
    by PATID;
    if a;
run;

/* Overall mean trajectory of BPRS over time */
proc sgplot data=alzheimer_analysis;
    title "Mean BPRS Trajectory Over Time";
    vline TIME / response=BPRS stat=mean markers 
                 limitstat=stderr limits=both;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="Mean BPRS";
run;

/* Mean trajectories by TRIAL */
proc sgplot data=alzheimer_analysis;
    title "Mean BPRS Trajectories by Trial Group";
    vline TIME / response=BPRS stat=mean group=TRIAL markers;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="Mean BPRS";
run;

/* Mean trajectories by SEX */
proc sgplot data=alzheimer_analysis;
    title "Mean BPRS Trajectories by Sex";
    vline TIME / response=BPRS stat=mean group=SEX markers;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="Mean BPRS";
run;

/* Spaghetti plot - individual trajectories */
proc sgplot data=alzheimer_analysis;
    title "Individual BPRS Trajectories (Sample)";
    where PATID <= 50; /* Show subset for clarity */
    series x=TIME y=BPRS / group=PATID transparency=0.7;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="BPRS";
run;

/* Mean trajectories by baseline CDRSB tertiles */
proc rank data=alzheimer_analysis groups=3 out=alzheimer_tertiles;
    var CDRSB0;
    ranks CDRSB0_tertile;
run;

proc sgplot data=alzheimer_tertiles;
    title "Mean BPRS Trajectories by Baseline CDRSB Tertiles";
    vline TIME / response=BPRS stat=mean group=CDRSB0_tertile markers;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="Mean BPRS";
    keylegend / title="CDRSB0 Tertile";
run;


/* Variance of BPRS over time */
proc sgplot data=alzheimer_analysis;
    title "Standard Deviation of BPRS Over Time";
    vline TIME / response=BPRS stat=std markers;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="SD of BPRS";
run;

/* Box plots to show variance at each time point */
proc sgplot data=alzheimer_analysis;
    title "Distribution of BPRS at Each Time Point";
    vbox BPRS / category=TIME;
    xaxis label="Time (years)";
    yaxis label="BPRS";
run;

/* Coefficient of variation over time */
proc means data=alzheimer_analysis noprint nway;
    class TIME;
    var BPRS;
    output out=variance_stats mean=mean_bprs std=sd_bprs cv=cv_bprs;
run;

data variance_stats;
    set variance_stats;
    cv_bprs_pct = cv_bprs * 100;
run;

proc sgplot data=variance_stats;
    title "Coefficient of Variation of BPRS Over Time";
    scatter x=TIME y=cv_bprs_pct / markerattrs=(size=10);
    series x=TIME y=cv_bprs_pct / lineattrs=(thickness=2);
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="CV (%)";
run;


/* Reshape to wide format for correlation analysis */
proc sort data=alzheimer_analysis; by PATID TIME; run;

proc transpose data=alzheimer_analysis out=bprs_wide prefix=BPRS_T;
    by PATID;
    id TIME;
    var BPRS;
run;

/* Calculate correlations between time points */
proc corr data=bprs_wide outp=corr_matrix noprint;
    var BPRS_T0-BPRS_T6;
run;

/* Extract correlation coefficients */
data corr_data;
    set corr_matrix;
    where _TYPE_='CORR';
    TIME1 = input(substr(_NAME_,7), 1.);
    array bprs_vars[*] BPRS_T0-BPRS_T6;
    do i = 1 to dim(bprs_vars);
        TIME2 = i - 1;
        CORRELATION = bprs_vars[i];
        LAG = abs(TIME2 - TIME1);
        output;
    end;
    keep TIME1 TIME2 CORRELATION LAG;
run;

/* Heatmap of correlations */
proc sgplot data=corr_data;
    title "Correlation Heatmap: BPRS Between Time Points";
    heatmap x=TIME1 y=TIME2 / colorresponse=CORRELATION 
            colormodel=(blue white red);
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="Time (years)" values=(0 to 6 by 1);
run;

/* Correlation vs lag plot */
proc means data=corr_data noprint nway;
    where TIME1 < TIME2;
    class LAG;
    var CORRELATION;
    output out=lag_corr mean=mean_corr std=sd_corr;
run;

proc sgplot data=lag_corr;
    title "Mean Correlation by Time Lag";
    scatter x=LAG y=mean_corr / markerattrs=(size=12);
    series x=LAG y=mean_corr / lineattrs=(thickness=2);
    xaxis label="Time Lag (years)" values=(0 to 6 by 1);
    yaxis label="Mean Correlation";
    refline 0 / axis=y lineattrs=(pattern=dash);
run;


/* BPRS vs baseline continuous predictors */
proc sgscatter data=alzheimer_analysis;
    where TIME=0;
    title "BPRS at Baseline vs Continuous Predictors";
    plot BPRS*(AGE EDU BMI CDRSB0 ABPET0 TAUPET0);
run;

/* Panel plot: BPRS trajectories by categorical variables */
/* proc sgpanel data=alzheimer_analysis; */
/*     title "BPRS Trajectories by Key Categorical Variables"; */
/*     panelby TRIAL SEX / rows=2 columns=2; */
/*     vline TIME / response=BPRS stat=mean markers; */
/*     rowaxis label="Mean BPRS"; */
/*     colaxis label="Time (years)"; */
/* run; */

/* ============================================
   5. VARIOGRAM FOR CORRELATION STRUCTURE
   ============================================ */

/* Calculate empirical variogram */
proc sort data=alzheimer_analysis; by PATID TIME; run;

data variogram_prep;
    set alzheimer_analysis;
    by PATID;
    retain bprs_lag1-bprs_lag6 time_lag1-time_lag6;
    array bprs_lags[6] bprs_lag1-bprs_lag6;
    array time_lags[6] time_lag1-time_lag6;
    
    if first.PATID then do i=1 to 6;
        bprs_lags[i] = .;
        time_lags[i] = .;
    end;
    
    do i=6 to 2 by -1;
        bprs_lags[i] = bprs_lags[i-1];
        time_lags[i] = time_lags[i-1];
    end;
    
    bprs_lags[1] = BPRS;
    time_lags[1] = TIME;
    
    do i=2 to 6;
        if not missing(bprs_lags[i]) and not missing(BPRS) then do;
            LAG_TIME = TIME - time_lags[i];
            SEMI_VAR = 0.5 * (BPRS - bprs_lags[i])**2;
            output;
        end;
    end;
    keep PATID TIME LAG_TIME SEMI_VAR;
run;

proc means data=variogram_prep noprint nway;
    class LAG_TIME;
    var SEMI_VAR;
    output out=variogram mean=mean_semivar n=n_pairs;
run;

proc sgplot data=variogram;
    title "Empirical Variogram of BPRS";
    scatter x=LAG_TIME y=mean_semivar / markerattrs=(size=12);
    series x=LAG_TIME y=mean_semivar / lineattrs=(thickness=2);
    xaxis label="Time Lag (years)";
    yaxis label="Semivariance";
run;

/* Summary statistics table */
proc means data=alzheimer_analysis n mean std min max;
    title "Summary Statistics by Time Point";
    class TIME;
    var BPRS;
run;
