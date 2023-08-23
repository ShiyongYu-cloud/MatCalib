function mcmcsamples = mcmc(theta,alpha,beta,P,A,B,year_scale,CalCurves,boundary_relationship,nsamples,thin,burn_in)
%% function for inferring model parameters using the MCMC method 
%INPUT
%theta: structure containing the old state of theta
%alpha: vector containing the old state of alpha
%beta: vector containing the old state of beta
%P: structure containing laboratory data of the periods
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%year_scale: type of year scale to be reported ('BP' or 'BCE/CE')
%CalCurves: structure containing the calibration curves in the F14C space
%boundary_relationship: relationship of two neighboring periods (overlapping|contiguous|disjoint)
%nsamples: numbers of MCMC runs
%thin: interval in every which random samples to be kept
%burn_in: steps after which random samples began to be kept
%OUTPUT
%mcmcsamples: randomly generated calendar ages of the raddiocaron ages and boundaries
%%
M = length(P);      %number of periods
N = zeros(1,M);     %number of calendar ages in each period
for i = 1:M
    N(i) = length(P(i).age); 
end
K = 2*M+sum(N);     %total number of parameters
mcmcsamples = zeros(nsamples,K);
%% mcmc run
h = waitbar(0,'Please wait...','CreateCancelBtn','closereq', 'name',...
    'MatCalib v2.0');
for i = 1-burn_in:nsamples*thin
    waitbar(i/(nsamples*thin+burn_in),h)
    theta = update_theta(theta,alpha,beta,P,A,B,year_scale,CalCurves);
    [alpha,beta] = update_boundary(alpha,beta,theta,boundary_relationship,year_scale,A,B);
    S = P2S(theta,alpha,beta);     %conver a structure to a vector for saving the results
    if (i > 0) && (mod(i,thin) == 0)
        mcmcsamples(i/thin,:) = S;
    end
end
close(h);
return;
%%
function S = P2S(theta,alpha,beta)
%% function for converting a structure to a vector of model paramters
%INPUT
%theta: structure containing the calendar ages of the periods 
%alpha: vector containing the calendar ages of ealy boundary of the periods
%beta: vector containing the calendar ages of late boundary of the periods
%OUTPUT
%S: vector containing all model paramters
%%
M = length(theta); %number of periods
N_theta = zeros(1,M);   
for i = 1:M
    N_theta(i) = length(theta(i).age);
end
S = zeros(1,2*M+sum(N_theta)); %set up the dimention of the age sequence
ind_alpha = cumsum([1 N_theta(1:end-1)+2]); %index of alpha in the age sequence
ind_b_theta = cumsum([1 N_theta(1:end-1)+2]) + 1; %index of beginning of theta in the sequence 
ind_e_theta = cumsum([1 N_theta(1:end-1)+2]) + N_theta; %index of end of theta in the sequence
ind_beta = cumsum([1 N_theta(1:end-1)+2]) + N_theta + 1;  %index of beta in the age sequence
for i = 1:M
    S(ind_alpha(i)) = alpha(i);
    S(ind_b_theta(i):ind_e_theta(i)) = theta(i).age;
    S(ind_beta(i)) = beta(i);
end
return;