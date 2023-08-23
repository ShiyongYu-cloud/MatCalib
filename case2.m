%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Bayesian radiocarbon age modeling                   % 
%                     MatCalib v2.0                                       % 
%                     Copyright: Shiyong Yu                               %                                  
%                     E-mail: shiyong.yu@gmail.com                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
%% Reading the radiocarbon ages to be modeled
DataFile = 'Baltic_data.txt';
P = read_data(DataFile);
%% Specifying depth unit
depth_unit = 'cm'; 
%% Specifying the relationship of every two neighboring phases (overlapping|contiguous|disjoint)
boundary_relationship = {'disjoint'};  
%% Specifying depths of phase boundaries
alpha_depth = [-850,-575];
beta_depth = [-575,0];
%% Setting up model parameters
delta = 5;                 % round ages to the nearest 5 years
year_scale = 'BP';         % scale of year to be reported (BC/AD or BP)
nchains = 3;               % number of Markov Chains to be run
nsamples = 1000;           % number of random samples to be saved
burn_in = 1000;            % steps of MCMC run when samples will be saved  
thin = 20;                 % steps in every which samples will be saved 
%% Modeling the radiocarbon ages
[samples,pdfs,ALPHA,BETA,THETA,CAGES,A,B] = age_modeling(P,boundary_relationship,delta,year_scale,nchains,nsamples,burn_in,thin);
%% Plotting the modeled ages 
figure(1)
plot_age_depth(pdfs,CAGES,P,alpha_depth,beta_depth,depth_unit,year_scale,A,B);
%% Calculating the differnce between two modeled ages and saving the results
x1 = samples(:,5);
x2 = samples(:,6);
bin_size = 5;
[diff_pdfs,AGE_DIFF] = age_difference(x1,x2,bin_size);
save_age_diff('diff_pdfs.txt',diff_pdfs);
%% Plotting the difference of two modeled ages
figure(2)
plot_difference(diff_pdfs);