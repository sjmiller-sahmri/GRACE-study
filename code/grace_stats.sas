#FDR correction performed using GraphPad Prism
#SAS language


##Fig. 2E-H
#Create ranks
data grace; set grace;
  if total_rpkm<=889.3159 then total_rpkm=.;
run;
proc rank data=grace out=gracee groups=3;
  var total_rpkm;
  ranks RPKM_Ranks;
run;
data grace; set grace;
  RPKM_Ranks=RPKM_Ranks + 1;
run;
data grace; set grace;
  if RPKM_Ranks=. then RPKM_Ranks=0;
run;

proc rank data=grace out=grace groups=3;
  var days_since_most_recent_abx;
  ranks most_recent_abx_days_rank;
run;
data grace; set grace;
  if most_recent_abx_days_rank=. then most_recent_abx_days_rank=3;
  most_recent_abx_days_rank_char=put(most_recent_abx_days_rank, 1.);
run;

proc rank data=grace out=grace groups=3;
  var num_exposure_events;
  ranks num_exp_events_ranks;
run;
data grace; set grace;
  if num_exp_events_ranks=. then num_exp_events_ranks=3;
  num_exp_events_ranks_char=put(num_exp_events_ranks, 1.);
run;

proc rank data=grace out=grace groups=3;
  var num_unique_classes;
  ranks num_unique_classes_rank;
run;
data grace; set grace;
  if num_unique_classes_rank=. then num_unique_classes_rank=3;
  num_unique_classes_char=put(num_unique_classes_rank, 1.);
run;

proc rank data=grace out=grace groups=3;
  var total_exposed_days;
  ranks total_exposed_days_rank;
run;
data grace; set grace;
  if total_exposed_days_rank=. then total_exposed_days_rank=3;
  total_exposed_days_rank_char=put(total_exposed_days_rank, 1.);
run;

#Model
proc logistic data=grace;
  class num_unique_classes_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=num_unique_classes_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_unique_classes_char;
run;
proc logistic data=grace;
  class num_exp_events_ranks_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=num_exp_events_ranks_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_exp_events_ranks_char;
run;
proc logistic data=grace;
  class total_exposed_days_rank_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=total_exposed_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_exp_events_ranks_char;
run;
proc logistic data=grace;
  class most_recent_abx_days_rank_char (ref='3') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=most_recent_abx_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio most_recent_abx_days_rank_char;
run;


##Fig. 3
#Create ranks
data grace; set grace;
  if cephalosporin_sum_rpkm<=244.639 then cephalosporin_sum_rpkm=.;
run;
proc rank data=grace out=grace groups=3;
  var cephalosporin_sum_rpkm;
  ranks Ceph_RPKM_Ranks;
run;
data grace; set grace;
  Ceph_RPKM_Ranks=Ceph_RPKM_Ranks + 1;
  if Ceph_RPKM_Ranks=. then Ceph_RPKM_Ranks=0;
run;

data grace; set grace;
  if diaminopyrimidine_sum_rpkm<=6.3851 then diaminopyrimidine_sum_rpkm=.;
run;
proc rank data=grace out=grace groups=3;
  var diaminopyrimidine_sum_rpkm;
  ranks Diamino_RPKM_Ranks;
run;
data grace; set grace;
  Diamino_RPKM_Ranks=Diamino_RPKM_Ranks + 1;
  if Diamino_RPKM_Ranks=. then Diamino_RPKM_Ranks=0;
run;

data grace; set grace;
  if penam_sum_rpkm<=229.5551 then penam_sum_rpkm=.;
run;
proc rank data=grace out=grace groups=3;
  var penam_sum_rpkm;
  ranks Penam_RPKM_Ranks;
run;
data grace; set grace;
  Penam_RPKM_Ranks=Penam_RPKM_Ranks + 1;
  if Penam_RPKM_Ranks=. then Penam_RPKM_Ranks=0;
run;

data grace; set grace;
  if tetracycline_sum_rpkm<=136.6433 then tetracycline_sum_rpkm=.;
run;
proc rank data=grace out=grace groups=3;
  var tetracycline_sum_rpkm;
  ranks Tetra_RPKM_Ranks;
run;
data grace; set grace;
  Tetra_RPKM_Ranks=Tetra_RPKM_Ranks + 1;
  if Tetra_RPKM_Ranks=. then Tetra_RPKM_Ranks=0;
run;

#Model
/** Any cephalosporin use **/
/* vs total */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No');
  model RPKM_Ranks (descending)=any_cephalosporin_used;
  oddsratio any_cephalosporin_used;
run;
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs cephalosporin only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_RPKM_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_RPKM_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_RPKM_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_RPKM_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;

/** Any Diaminopyrimidine use **/
/* vs total */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No');
  model RPKM_Ranks (descending)=any_diaminopyrimidine_used;
  oddsratio any_diaminopyrimidine_used;
run;
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Cephalosporin only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_RPKM_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_RPKM_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_RPKM_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_RPKM_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;

/** Any Penicillin use **/
/* vs total */
proc logistic data=grace;
  class any_penam_used (ref='No');
  model RPKM_Ranks (descending)=any_penam_used;
  oddsratio any_penam_used;
run;
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Cephalosporin only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_RPKM_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_RPKM_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_RPKM_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_RPKM_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;

/** Any Tetracycline use **/
/* vs total */
proc logistic data=grace;
  class any_tetracycline_used (ref='No');
  model RPKM_Ranks (descending)=any_tetracycline_used;
  oddsratio any_tetracycline_used;
run;
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model RPKM_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Cephalosporin only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_RPKM_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_RPKM_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_RPKM_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_RPKM_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;


##Fig. 5
#Create ranks
data grace_concern; set grace_concern;
  if total_rpkm20<=1.95045 then total_rpkm20=.;
run;
proc rank data=grace_concern out=grace_concern groups=3;
  var total_rpkm20;
  ranks total_rpkm20_ranks;
run;
data grace_concern; set grace_concern;
  total_rpkm20_ranks=total_rpkm20_ranks2 + 1;
  if total_rpkm20_ranks=. then total_rpkm20_ranks2=0;
run;

proc rank data=grace_concern out=grace_concern groups=4;
	var Aminoglycoside;
	ranks Rank_Aminoglycoside;
run;

proc rank data=grace_concern out=grace_concern groups=4;
	var AmpC;
	ranks Rank_AmpC;
run;

proc rank data=grace_concern out=grace_concern groups=4;
	var ESBL;
	ranks Rank_ESBL;
run;

proc rank data=grace_concern out=grace_concern groups=4;
	var Betalactam;
	ranks Rank_Betalactam;
run;

proc rank data=grace_concern out=grace_concern groups=4;
	var Glycopeptide;
	ranks Rank_Glycopeptide;
run;

proc rank data=grace_concern out=grace_concern groups=4;
	var Linezolid;
	ranks Rank_Linezolid;
run;

proc rank data=grace_concern out=grace_concern groups=4;
	var Fluoroquinolone;
	ranks Rank_Fluoroquinolone;
run;

#Model
/* Total ranked abundance */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model total_RPKM20_Ranks=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model total_RPKM20_Ranks=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model total_RPKM20_Ranks=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model total_RPKM20_Ranks=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model total_RPKM20_Ranks=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;

/* AmpC */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_AmpC any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_AmpC=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_AmpC any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_AmpC=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_AmpC any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_AmpC=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_AmpC any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_AmpC=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_AmpC any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_AmpC=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;

/* ESBL */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_ESBL any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_ESBL=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_ESBL any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_ESBL=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex 
		days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_ESBL any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_ESBL=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_ESBL any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_ESBL=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_ESBL any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_ESBL=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;

/* Beta-lactam */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Betalactam any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Betalactam=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Betalactam any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Betalactam=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex 
		days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Betalactam any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Betalactam=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Betalactam any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Betalactam=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Betalactam any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Betalactam=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;

/* Linezolid */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Linezolid any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Linezolid=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Linezolid any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Linezolid=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex 
		days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Linezolid any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Linezolid=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Linezolid any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Linezolid=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Linezolid any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Linezolid=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;

/* Fluoroquinolone */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Fluoroquinolone any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Fluoroquinolone=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Fluoroquinolone any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Fluoroquinolone=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex 
		days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Fluoroquinolone any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Fluoroquinolone=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Fluoroquinolone any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Fluoroquinolone=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Fluoroquinolone any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Fluoroquinolone=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;

/* Glycopeptide */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Glycopeptide any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Glycopeptide=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Glycopeptide any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Glycopeptide=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex 
		days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Glycopeptide any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Glycopeptide=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Glycopeptide any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Glycopeptide=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class Rank_Glycopeptide any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Glycopeptide=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;

/* Aminoglycoside - ranked abundance */
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Aminoglycoside=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Aminoglycoside=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Aminoglycoside=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Aminoglycoside=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;
proc logistic data=grace_concern plots(only)=(effect oddsratio) alpha=0.05;
	class any_abx_yn (ref='Yes') depression gastro_ref_disease hypertension pain sex;
	model Rank_Aminoglycoside=any_abx_yn depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
run;


##Supplementary Fig. 3
/* TETRACYCLINE EXPOSURE */
%macro analyze_genes;
%do i = 1 %to 359;

    proc rank data=grace_individual out=grace_individual_ranked groups=3;
        var gene&i;
        ranks gene&i._Ranks;
    run;
    
    data grace_individual_ranked; set grace_individual_ranked;
		    gene&i._Ranks=gene&i._Ranks + 1;
		    if gene&i._Ranks=. then gene&i._Ranks=0;
  	run;

    proc logistic data=grace_individual_ranked;
        class any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
        model gene&i._Ranks = any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
        ods output ParameterEstimates=Gene&i._Estimates;
    run;
%end;
%mend;
%analyze_genes;


##Supplementary Fig. 4
/* TETRACYCLINE EXPOSURE */
%macro analyze_taxa;
%do i = 1 %to 238;

    proc rank data=grace_individual out=grace_individual_ranked groups=3;
        var taxa&i;
        ranks taxa&i._Ranks;
    run;
    
    data grace_individual_ranked; set grace_individual_ranked;
		    taxa&i._Ranks=taxa&i._Ranks + 1;
		    if taxa&i._Ranks=. then taxa&i._Ranks=0;
  	run;

    proc logistic data=grace_individual_ranked;
        class any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
        model taxa&i._Ranks = any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
        ods output ParameterEstimates=taxa&i._Estimates;
    run;
%end;
%mend;
%analyze_taxa;


##Supplementary Fig. 6
proc logistic data=grace_concern;
  class most_recent_abx_days_rank_char (ref='3') depression gastro_ref_disease hypertension pain sex site;
  model total_RPKM20_Ranks (descending)=most_recent_abx_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio most_recent_abx_days_rank_char;
run;

proc logistic data=grace_concern;
  class total_exposed_days_rank_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model total_RPKM20_Ranks (descending)=total_exposed_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio total_exposed_days_rank_char;
run;

proc logistic data=grace_concern;
  class num_exp_events_ranks_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model total_RPKM20_Ranks (descending)=num_exp_events_ranks_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_exp_events_ranks_char;
run;

proc logistic data=grace_concern;
  class num_unique_classes_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model total_RPKM20_Ranks (descending)=num_unique_classes_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_unique_classes_char;
run;


##Supplementary Fig. 7
/** Cephalosporin exposure **/
%macro analyze_genes;
%do i = 1 %to 9;

    proc rank data=grace_concern out=grace_concern_ranked groups=3;
        var concern_gene_&i;
        ranks concern_gene_&i._Ranks;
    run;
    
    data grace_concern_ranked; set grace_concern_ranked;
		    concern_gene_&i._Ranks=concern_gene_&i._Ranks + 1;
		    if concern_gene_&i._Ranks=. then concern_gene_&i._Ranks=0;
  	run;

    proc logistic data=grace_concern_ranked;
        class any_cephalosporin_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
        model concern_gene_&i._Ranks = any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
        ods output ParameterEstimates=concern_gene_&i._Estimates;
    run;
%end;
%mend;
%analyze_genes; /* Run the macro */

/** Diaminopyrimidine exposure **/
%macro analyze_genes;
%do i = 1 %to 9;

    proc rank data=grace_concern out=grace_concern_ranked groups=3;
        var concern_gene_&i;
        ranks concern_gene_&i._Ranks;
    run;
    
    data grace_concern_ranked; set grace_concern_ranked;
		    concern_gene_&i._Ranks=concern_gene_&i._Ranks + 1;
		    if concern_gene_&i._Ranks=. then concern_gene_&i._Ranks=0;
	  run;

    proc logistic data=grace_concern_ranked;
        class any_diaminopyrimidine_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
        model concern_gene_&i._Ranks = any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
        ods output ParameterEstimates=concern_gene_&i._Estimates;
    run;
%end;
%mend;
%analyze_genes; /* Run the macro */

/** Penicillin (penam) exposure **/
%macro analyze_genes;
%do i = 1 %to 9;

    proc rank data=grace_concern out=grace_concern_ranked groups=3;
        var concern_gene_&i;
        ranks concern_gene_&i._Ranks;
    run;
    
    data grace_concern_ranked; set grace_concern_ranked;
		    concern_gene_&i._Ranks=concern_gene_&i._Ranks + 1;
		    if concern_gene_&i._Ranks=. then concern_gene_&i._Ranks=0;
	  run;

    proc logistic data=grace_concern_ranked;
        class any_penam_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
        model concern_gene_&i._Ranks = any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
        ods output ParameterEstimates=concern_gene_&i._Estimates;
    run;
%end;
%mend;
%analyze_genes; /* Run the macro */

/** Tetracycline exposure **/
%macro analyze_genes;
%do i = 1 %to 9;

    proc rank data=grace_concern out=grace_concern groups=3;
        var concern_gene_&i;
        ranks concern_gene_&i._Ranks;
    run;
    
    data grace_concern_ranked; set grace_concern;
		    concern_gene_&i._Ranks=concern_gene_&i._Ranks + 1;
		    if concern_gene_&i._Ranks=. then concern_gene_&i._Ranks=0;
	  run;

    proc logistic data=grace_concern;
        class any_tetracycline_used (ref='Yes') depression gastro_ref_disease hypertension pain sex;
        model concern_gene_&i._Ranks = any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry;
        ods output ParameterEstimates=concern_Gene_&i._Estimates;
    run;
%end;
%mend;
%analyze_genes; 


##Supplementary Table 8
#Create ranks
proc rank data=grace out=grace groups=4;
var number_of_genes;
ranks ARG_Ranks;
run;

proc rank data=grace out=grace groups=4;
var Shannon;
ranks Shannon_Ranks;
run;

#Model
/* Richness */
proc logistic data=grace;
  class num_exp_events_ranks_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=num_exp_events_ranks_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_exp_events_ranks_char;
run;
proc logistic data=grace;
  class most_recent_abx_days_rank_char (ref='3') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=most_recent_abx_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio most_recent_abx_days_rank_char;
run;
proc logistic data=grace;
  class total_exposed_days_rank_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=total_exposed_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio total_exposed_days_rank_char;
run;
proc logistic data=grace;
  class num_unique_classes_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=num_unique_classes_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_unique_classes_char;
run;

/* Shannon's Diversity */
proc logistic data=grace;
  class num_exp_events_ranks_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model Shannon_Ranks (descending)=num_exp_events_ranks_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_exp_events_ranks_char;
run;
proc logistic data=grace;
  class most_recent_abx_days_rank_char (ref='3') depression gastro_ref_disease hypertension pain sex site;
  model Shannon_Ranks (descending)=most_recent_abx_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio most_recent_abx_days_rank_char;
run;
proc logistic data=grace;
  class total_exposed_days_rank_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model Shannon_Ranks (descending)=total_exposed_days_rank_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio total_exposed_days_rank_char;
run;
proc logistic data=grace;
  class num_unique_classes_char (ref='0') depression gastro_ref_disease hypertension pain sex site;
  model Shannon_Ranks (descending)=num_unique_classes_char depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio num_unique_classes_char;
run;


##Supplementary Table 9
#Create ranks
proc rank data=grace out=grace groups=4;
var cephalosporin_num_args;
ranks Ceph_ARG_Ranks;
run;
proc rank data=grace out=grace groups=4;
var diaminopyrimidine_num_args;
ranks Diamino_ARG_Ranks;
run;
proc rank data=grace out=grace groups=4;
var penam_num_args;
ranks Penam_ARG_Ranks;
run;
proc rank data=grace out=grace groups=4;
var tetracycline_num_args;
ranks Tetra_ARG_Ranks;
run;

#Model
/** Any cephalosporin use **/
/* vs total */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs cephalosporin only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_ARG_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_ARG_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_ARG_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_cephalosporin_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_ARG_Ranks (descending)=any_cephalosporin_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_cephalosporin_used;
run;

/** Any Diaminopyrimidine use **/
/* vs total */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Cephalosporin only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_ARG_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_ARG_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_ARG_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_diaminopyrimidine_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_ARG_Ranks (descending)=any_diaminopyrimidine_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_diaminopyrimidine_used;
run;

/** Any Penicillin use **/
/* vs total */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Cephalosporin only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_ARG_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_ARG_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_ARG_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_penam_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_ARG_Ranks (descending)=any_penam_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_penam_used;
run;

/** Any Tetracycline use **/
/* vs total */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model ARG_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Cephalosporin only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Ceph_ARG_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Diaminopyrimidine only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Diamino_ARG_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Penicillin only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Penam_ARG_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
/* vs Tetracycline only */
proc logistic data=grace;
  class any_tetracycline_used (ref='No') depression gastro_ref_disease hypertension pain sex site;
  model Tetra_ARG_Ranks (descending)=any_tetracycline_used depression gastro_ref_disease hypertension pain age_at_enrolment_years sex days_since_entry site;
  oddsratio any_tetracycline_used;
run;
