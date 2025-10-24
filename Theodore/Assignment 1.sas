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