/*

* from wide to long format;

libname File '/home/u64347574'; 
* create the library where the file is stored; 

data alzheimer25;
    set File.alzheimer25;   
    * create a temporary copy in WORK;
run;


data dati_long;
    set alzheimer25;
    
    array cdrsb_array[7] cdrsb0-cdrsb6;  
    * array for the first variable;
    array bprs_array[7] bprs0-bprs6;  
    * array for the second variable;
    array abpet_array[7] abpet0-abpet6;  
    * array for the third variable;
    array taupet_array[7] taupet0-taupet6;
    
    do time = 1 to 7;
        cdrsb = cdrsb_array[time];   
        * new value of variable v;
        bprs = bprs_array[time];   
        * new value of variable w;
        abpet = abpet_array[time];
        taupet = taupet_array[time];  
        * new value of variable x;
        output;               
        * create a new row;
    end;

    drop cdrsb0-cdrsb6 bprs0-bprs6 abpet0-abpet6 taupet0-taupet6;   
    * remove the original columns;
run;

* dataset reorganized! Now let's make some plots;

* Import the CSV file;
proc import out=File.DATI_LONG_csv; 
    * this is where it will be saved;
    datafile="/home/u64347574/DATI_LONG.csv"; 
    * this is where it will be taken from;
    dbms=csv replace;
    getnames=yes;
run;

proc contents data=work.DATI_LONG_csv;
run;

* CSV imported, try making some plots;

* NOW FORMAT OKAY, READY FOR ANALYSIS;

*/

* now start the analysis;

libname File '/home/u64347574'; /* create the library where the file is stored */

* this code creates a temporary file in WORK, not necessary, you can also work directly on FILE;
data dati_long_csv;
    set File.dati_long_csv;   /* create a temporary copy in File */
run;

/*
proc sgpanel data=dati_long_csv;
    panelby job ;
    series x=time y=cdrsb / group=patid lineattrs=(thickness=1 pattern=solid);
    scatter x=time y=cdrsb / group=patid;
    colaxis label="Time (years)";
    rowaxis label="CDRSB";
run;
*/

/* time evolution with respect to education level */

proc sort data=dati_long_csv nodupkey out=unique_patients;
    by edu patid;
run;

proc surveyselect data=unique_patients 
    out=sample_ids
    method=srs
    sampsize=(10 10 10 10 10)     /* 10 for each education level */
    seed=12345;
    strata edu;
    id patid;
run;

proc sql;
    create table dati_sample as
    select *
    from dati_long_csv
    where patid in (select patid from sample_ids);
quit;

proc sgpanel data=dati_sample;
    panelby edu ;
    series x=time y=cdrsb / group=patid lineattrs=(thickness=1 pattern=solid);
    colaxis label="Time (years)";
    rowaxis label="CDRSB";
run;
