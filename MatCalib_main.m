%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Bayesian radiocarbon age modeling                   % 
%                     MatCalib v2.0                                       % 
%                     Copyright: Shiyong Yu                               %                                  
%                     E-mail: shiyong.yu@gmail.com                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
%% Reading radiocarbon ages to be modeled
DataFile = 'data_file_template.txt';
P = read_data(DataFile);
%% Specifying the relationship of every two neighboring periods (overlapping|contiguous|disjoint)
boundary_relationship = {'overlapping','disjoint'};    
%% Setting up model parameters
delta = 5;                 % years to which the modeled ages to be rounded
year_scale = 'BCE/CE';     % scale of year to be reported (BCE/CE or BP)
nchains = 2;               % number of Markov Chains to be run
nsamples = 500;           % number of random samples to be saved
burn_in = 1000;            % steps of MCMC run when samples will be saved  
thin = 20;                 % steps in every which samples will be saved 
%% Modeling the radiocarbon ages
[samples,pdfs,ALPHA,BETA,THETA,CAGES,A,B] = age_modeling(P,boundary_relationship,delta,year_scale,nchains,nsamples,burn_in,thin);
%% Plotting the pdfs of modeled ages
figure(1)
plot_ages(pdfs,CAGES,year_scale,A,B);
%%
%% Calculating the differnce between two modeled ages and saving the results
x1 = samples(:,13);
x2 = samples(:,14);
bin_size = 5;
[diff_pdfs,AGE_DIFF] = age_difference(x1,x2,bin_size);
save_age_diff('diff_pdfs.txt',diff_pdfs);
%% Plotting the difference of two modeled ages
figure(2)
plot_difference(diff_pdfs);
%% Calculating the pooled mean of modeled ages in a period and saving the results
period_number = 3; % Period 3
[pooled_pdfs,POOLED_AGE] = pooled_mean(period_number,pdfs,P,delta);
save_pooled_mean('pooled_pdfs.txt',pooled_pdfs);
%% Plotting the pooled mean ages in a period
figure(3)
plot_pooled_mean(pooled_pdfs,POOLED_AGE,year_scale);