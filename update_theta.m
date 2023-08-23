function theta_out = update_theta(theta,alpha,beta,P,A,B,year_scale,CalCurves)
%% function for updating the calendar ages of radiocarbon ages in a sequence
%INPUT
%theta: structure containing the old status of the calendar ages in a sequence
%alpha: vector containing the calendar age of early boundary of all periods
%beta: vector containing the calendar age of late boundary of all periods
%P: structure containing the data of all periods
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%year_scale: scale of year (e.g. BCE/CE or BP)
%CalCurve: sturcture containing the calibration curves
%OUTPUT 
%theta_out: structure containing the new state of the calendar ages in a
%sequence
%%
theta_star = struct('age',[]);
theta_out = struct('age',[]);
M = length(theta); %number of periods
N = zeros(1,M);
for i = 1:M
    N(i) = length(theta(i).age);
end   
% propose a move
for i = 1:M
     if strcmpi(P(i).ordering,'coeval') == 1
        theta_star(i).age = round(theta(i).age + 1*(2*rand(1,1)-1)); 
     else 
        theta_star(i).age = round(theta(i).age + 1*(2*rand(1,N(i))-1));
     end   
end    
% calulate the acceptance probability
num = likelihood(P,theta_star,A,B,year_scale,CalCurves)+theta_prior(theta_star,alpha,beta,P,year_scale);
den = likelihood(P,theta,A,B,year_scale,CalCurves)+theta_prior(theta,alpha,beta,P,year_scale);
omega = num-den;
omega = min(omega,0);
rho = log(rand(1,1));
if omega > rho
    for i = 1:M
        theta_out(i).age = theta_star(i).age;
    end    
else
    for i = 1:M
        theta_out(i).age = theta(i).age;
    end   
end
end