# combo-efficacy-safety
Supporting MATLAB code for CCR submitted manuscript: Dose Optimization for Anticancer Drug Combinations: Maximizing Therapeutic Index via Simultaneous Clinical Exposure-Toxicity/Preclinical Efficacy modeling

do_everything.m
Does almost everything. This is a script file that generates synthetic mouse efficacy and clinical toxicity data and then performs appropriate analyses with these data corresponding to some of the manuscript figures. This file calls fit_loewe_combo_tox_data2.m and unilogifit.m functions, also included. 

TGIvsGRI.m
Generates the supplemental figure comparing the TGI and GRI xenograft response metrics

fig_toxmodel.m
Generates supplementary figure that shows which alpha level corresponds to non-overlapping tox
