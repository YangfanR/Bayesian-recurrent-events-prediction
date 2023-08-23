# A Bayesian framework for event prediction in clinical trials with recurrent event endpoints and terminal events

We employ a Bayesian framework based on a joint frailty model for prediction of the timing of observing the desired number of total events. Patient enrollment and censoring of patients due to other reasons are also modeled in the Bayesian predictive framework. The proposed approach is illustrated by a simulated case study, where predictive quantities informative for trial monitoring and interim decision making are highlighted.

main_code.R - main file of the simulated case analysis
func_analysis.R - function of conducting Bayesian analysis using rjags
PH_case1.txt - rjags function
pred_rec_exp.R - function of prediction using the result from analysis
data.Rdata - simulated data
