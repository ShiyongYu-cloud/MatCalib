function [samples,pdfs,ALPHA,BETA,THETA,CAGES,A,B] = age_modeling(P,boundary_relationship,delta,year_scale,nchains,nsamples,burn_in,thin)
%% function for Bayesian age modeling
%INPUT
%P: structure containing the information of all radiocarbon ages
%boundary_relationship: vector containing relationship of every two neighboring periods 
%delta: years to which the modeled ages to be rounded
%year_scale: scale of year to be reported (BCE/CE or BP)
%nchains: number of Markov Chains to be run
%nsamples: number of random samples to be saved
%burn_in: steps of MCMC run when samples will be saved 
%thin: steps in every which samples will be saved
%OUTPUT:
%samples: matrix containing randon samples of the model parameters
%pdfs: matrix containing the pdf of the model parameters
%ALPHA: structure containing the 68.2%, %95.4% hpd regions, and the median prabability age 
%BETA:  structure containing the 68.2%, %95.4% hpd regions, and the median prabability age 
%THETA: structure containing the 68.2%, %95.4% hpd regions, and the median prabability age 
%CAGES: structure wrapping up ALPHA, BETA, and THETA
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%%
M = length(P);             % number of periods in the sequence  
N = zeros(1,M);            % number of radiocarbon ages in each period
for i = 1:M
    N(i) = length(P(i).age); 
end
K = 2*M+sum(N);            % total number of model parameters
mcmcsamples = zeros(nsamples,K,nchains);
%% Converting 14C ages to F14C space and correcting for reservoir ages if any
for i = 1:M
    [P(i).f14C,P(i).ferr] = radiocarbon2f(P(i).age,P(i).err,P(i).rage,P(i).rerr);
end
%% Estimating the early and late boundary of the sequence
[A,B] = sequence_bound(P,year_scale);
%% Reading the calibration curve and convert it to F14C space
CalCurves = read_curves(P,A,B,year_scale);                     
%% Inferring model parameters and saving the results 
disp('Begin Bayesian age modeling...');
%parpool('local'); % start a local parallel pool
for i = 1:nchains
    disp(['Generating chain ' num2str(i) '...']);
    [theta,alpha,beta] = initialization(P,A,B,year_scale,CalCurves,boundary_relationship);
    mcmcsamples(:,:,i) = mcmc(theta,alpha,beta,P,A,B,year_scale,CalCurves,boundary_relationship,nsamples,thin,burn_in);
end
%delete(gcp('nocreate')); % close the pocal paralell pool
disp('Bayesian age modeling completed');
%% Monitoring convergence of the Markov chains for the parameters
R_hat = convergence(mcmcsamples);
%% Mixing the chains for each parameters and saving the results
samples = round(mean(mcmcsamples,3));
savemcmc('mcmc_run.txt',P,samples);
%% Estimating empirical pdf of the model parameters and saving the results
[pdfs,hpd68_2,hpd95_4] = mcmc2pdf(samples,A,B,delta,year_scale);
savepdf('epdf.txt',P,year_scale,pdfs);
savepdf('hpd68_2.txt',P,year_scale,hpd68_2);
savepdf('hpd95_4.txt',P,year_scale,hpd95_4);
%% Calculating 95.4% and 68.2% highest density regions and median probability age of the parameters
[ALPHA,BETA,THETA,CAGES] = pdf2hpd(pdfs,delta,P,year_scale); %results were saved in three structure arrays
%Table = struct2table(CAGES);
end