
/* Full multivariate model with unstructured covariance structure - works */
proc mixed data=alzheimer_long2 method=reml;
	class patid time wzc job edu;
	model bprs = time age job bmi adl wzc cdrsb abpet edu inkomen sex taupet age*time sex*time 
  				job*time bmi*time wzc*time cdrsb*time abpet*time inkomen*time adl*time/ solution;
  repeated time / type=un subject=patid ;

  
