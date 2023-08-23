%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Bayesian radiocarbon age modeling                   % 
%                     MatCalib v2.0                                       % 
%                     Copyright: Shiyong Yu                               %                                  
%                     E-mail: shiyong.yu@gmail.com                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
%% Reading the radiocarbon ages to be modeled
DataFile = 'archaeo_data.txt';
P = read_data(DataFile);
%% Specifying the relationship of every two neighboring periods (overlapping|contiguous|disjoint)
boundary_relationship = {'overlapping'};  
%% Setting up model parameters
delta = 5;                 % round ages to the nearest 5 years
year_scale = 'BCE/CE';     % scale of year to be reported (BCE/CE or BP)
nchains = 3;               % number of Markov Chains to be run
nsamples = 1000;           % number of random samples to be saved
burn_in = 1000;            % steps of MCMC run when samples will be saved  
thin = 20;                 % steps in every which samples will be saved 
%% Modeling the radiocarbon ages
[samples,pdfs,ALPHA,BETA,THETA,CAGES,A,B] = age_modeling(P,boundary_relationship,delta,year_scale,nchains,nsamples,burn_in,thin);
%% Plotting the pdfs of modeled ages
figure(1)
plot_ages(pdfs,CAGES,year_scale,A,B);
%% Calculating the differnce between two modeled ages and saving the results
x1 = samples(:,11);
x2 = samples(:,12);
bin_size = 5;
[diff_pdfs,AGE_DIFF] = age_difference(x1,x2,bin_size);
save_age_diff('diff_pdfs.txt',diff_pdfs);
%% Plotting the difference of two modeled ages
figure(2)
plot_difference(diff_pdfs);