# Bayesian_modelling_business_intelligence
This repository contains all codes and reports relevant to the project titled "Bayesian Modeling of Networks in Complex Business Intelligence Problems" done as a part of the course MTH422 : An introduction to Bayesian Analysis - IIT Kanpur. 

This project attempts to implement the algorithm proposed in - Bayesian modelling of networks in complex business intelligence problems D. Durante et. al.
The paper proposes a Bayesian hierarchical model to address complex business intelligence problems, specifically for optimising cross-selling strategies. The authors utilise data on customer product choices and co-subscription networks across multiple agencies. Their model clusters agencies with similar customer purchasing behaviors and then, within each cluster, employs a mixture of latent eigenmodels to understand customer preferences. This framework allows for the development of targeted cross-sell strategies and the evaluation of their potential effectiveness.

Through a pipeline of simulating data, recreating their MCMC algorithm and testing results we validate the claims of the authors and provide critical analysis of the strengths and weaknesses of the proposed algorithm.

# Guide
Sim.R : For dataset simulation
Algo.R : Contains code for the MCMC algorithm implemented. This also contains the code to reproduce the trace plot.
Visualizations.R : Produces the other plots and results in the report.
If variables seem to be missing while running Visualizations.R, please run Sim.R and Algo.R first
