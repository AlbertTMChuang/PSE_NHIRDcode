%let bb=104;
%let ee=109;

proc sql noprint;
	select name into : vars separated by " "
	from dictionary.columns
	where LIBNAME = upcase("final")
	and MEMNAME = upcase("cohort")
	and prxmatch("/^med_/i", name)>0 ;
quit;
%let med_=&vars.;
%put &med_;

proc sql noprint;
	select name into : vars separated by " "
	from dictionary.columns
	where LIBNAME = upcase("final")
	and MEMNAME = upcase("cohort")
	and prxmatch("/^cov_/i", name)>0 ;
quit;
%let cov_=&vars.;

%fine_stratification (in_data= final.cohort, exposure= g_sglt2,
PS_provided= no , 
ps_var= ps, 

ps_class_var_list=cov: med: id_s age_g,
ps_cont_var_list=, 
interactions= , 
PSS_method=cohort, 
estimand= ATE, 

n_of_strata= 50, 
out_data= PS_FS,
id_var= id, 
effect_estimate= hr, 
outcome= , 
survival_time= , 
time_unit= days, 

out_excel=,
work_lib=no);


PROC LOGISTIC DATA=final.cohort;
  CLASS cov: med: id_s age_g  g_sglt2(ref="0")  ;
  MODEL g_sglt2 = cov: med: id_s age_g;
  OUTPUT OUT=bb PREDICTED=PS;
RUN;

proc logistic data= final.cohort desc;
	class g_sglt2(ref="0");
	model g_sglt2=;
	output out=s_sglt2 p=sglt2;
run;

proc sort data= final.cohort;
	by id;run;
proc sort data=bb;
	by id;run;
proc sort data=s_sglt2;
	by id;run;
data  aa_m;
	merge final.cohort bb s_sglt2;
	by id;
	if g_sglt2=1 then uw=1/ps; else if g_sglt2=0 then uw=1/(1-ps);
	if g_sglt2=1 then sw=sglt2/ps; else if g_sglt2=0 then sw=(1-sglt2)/(1-ps);
	keep id ps sw;
run;
proc sql;
	create table final.weight as
	select a.*, b.psweight, c.* from final.cohort as a
	left join ps_fs as b
	on a.id=b.id
	left join aa_m as c
	on a.id=c.id;
quit;


%table1 (in_for_table1= final.cohort, treatment_var= g_sglt2, 
	categorical_var_list=&cov_ &med_ id_s age_g , continuous_var_list=age ,  out_table1= Crude_T1);

%macro aa(outcome);
	%let g=%sysfunc(countw(&outcome));
	%do i=1 %to %sysfunc(countw(&outcome));
		%let name  = %scan(&tt, &i, %str( ));
		%let pp  = %substr(&name, 5, %length(&name));
			proc sql; create table fine_counts as 
				select drug, 
				round(sum (&name.*psweight)) as n_events, round(sum(tr_&name.*psweight)/365) as cumulative_pyears, 
				(calculated n_events/calculated cumulative_pyears)*100 as IR_100py,
				(cinv(0.025,2*calculated n_events)/(2*calculated cumulative_pyears))*100 as IR_lcl,
				(cinv(0.975,2*(1+calculated n_events))/(2*calculated cumulative_pyears))*100 as IR_ucl
				from final.weight
				group by drug;
			quit;

			data fine_counts;
				format method $8.;
				set fine_counts;
				method="FINE";
			run;

			data event_counts&i; set crude_counts match1_counts IPTW_counts fine_counts; run;
			data event_counts&i;
				format outcome $20.;
				set event_counts&i;
				outcome="&pp";
			run;

			ods output ParameterEstimates=iptw_beta; 
			  proc phreg data=final.weight;
			  class drug(ref="dpp4");
			     model tr_&name.*&name.(0) = drug/rl;
			     weight sw;
			     run;
			ods output clear; 

			ods output ParameterEstimates=fine_beta; 
			  proc phreg data=final.weight;
			  class drug(ref="dpp4");
			     model tr_&name.*&name.(0) = drug/rl;
			     weight psweight;
			     run;
			ods output clear; 


			data results&i; length method outcome $8;  set crude_beta (in=r) match1_beta(in=m) iptw_beta(in=n) fine_beta(in=k); 
				if r=1 then Method= 'crude';
				else if k=1 then Method="FINE"; 
				outcome="&pp";
				drop chisq probchisq classval0 df;
			run;
	%end;
	data incidence;
		set event_counts1-event_counts&g.;
	run;
	data hrs;
		set results1-results&g.;
	run;
%mend;

%let tt=outcome_c outcome_e outcome_comb;
%aa(&tt.);
