function [alpha_out,beta_out] = update_boundary(alpha,beta,theta,boundary_relationship,year_scale,A,B)
%% function for updating the celendar ages of the boundaries of all periods
%INPUT
%alpha: vector containing the old state of alpha
%beta: vector containing the old state of beta
%theta: structure containing the old state of theta
%boundary_relationship: relationship of two neighboring periods (overlapping|contiguous|disjoint)
%year_scale: scale of year ('BCE/CE' or 'BP')
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%OUTPUT
%alpha_out: vector containing the new state of alpha
%beta_out: vector containing the new state of beta
%%  
M = length(theta); %number of periods
alpha_star = zeros(1,M);
beta_star = zeros(1,M);
q_alpha_num = zeros(1,M);
q_alpha_den = zeros(1,M);
q_beta_num = zeros(1,M);
q_beta_den = zeros(1,M);
%% propose a move for alpha
for i = 1:M
    if i == 1 % the first period
        if strcmpi(year_scale,'BCE/CE') == 1
            x = (A + alpha(i))/2;
            alpha_star(i) = trirnd(x,alpha(i),min(theta(i).age),1);
            q_alpha_num(i) = tripdf(alpha(i),x,alpha_star(i),min(theta(i).age));
            q_alpha_den(i) = tripdf(alpha_star(i),x,alpha(i),min(theta(i).age));
        elseif strcmpi(year_scale,'BP') == 1
            x = (A + alpha(i))/2;
            alpha_star(i) = trirnd(max(theta(i).age),alpha(i),x,1);
            q_alpha_num(i) = tripdf(alpha(i),max(theta(i).age),alpha_star(i),x);
            q_alpha_den(i) = tripdf(alpha_star(i),max(theta(i).age),alpha(i),x);
        end
    else %other periods
        if strcmpi(year_scale,'BCE/CE') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                id = theta(i-1).age < alpha(i);
                x = max(theta(i-1).age(id));
                alpha_star(i) = trirnd(x,alpha(i),min(theta(i).age),1);
                q_alpha_num(i) = tripdf(alpha(i),x,alpha_star(i),min(theta(i).age));
                q_alpha_den(i) = tripdf(alpha_star(i),x,alpha(i),min(theta(i).age));
            elseif strcmpi(boundary_relationship(i-1),'contiguous') == 1
                alpha_star(i) = trirnd(max(theta(i-1).age),alpha(i),min(theta(i).age),1);
                q_alpha_num(i) = tripdf(alpha(i),max(theta(i-1).age),alpha_star(i),min(theta(i).age));
                q_alpha_den(i) = tripdf(alpha_star(i),max(theta(i-1).age),alpha(i),min(theta(i).age));
            elseif strcmpi(boundary_relationship(i-1),'disjoint') == 1
                x = (beta(i-1) + alpha(i))/2;
                x = (x + alpha(i))/2;
                alpha_star(i) = trirnd(x,alpha(i),min(theta(i).age),1);
                q_alpha_num(i) = tripdf(alpha(i),x,alpha_star(i),min(theta(i).age));
                q_alpha_den(i) = tripdf(alpha_star(i),x,alpha(i),min(theta(i).age));
            end
        elseif strcmpi(year_scale,'BP') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                id = theta(i-1).age > alpha(i);
                x = min(theta(i-1).age(id));
                alpha_star(i) = trirnd(max(theta(i).age),alpha(i),x,1);
                q_alpha_num(i) = tripdf(alpha(i),max(theta(i).age),alpha_star(i),x);
                q_alpha_den(i) = tripdf(alpha_star(i),max(theta(i).age),alpha(i),x);
            elseif strcmpi(boundary_relationship(i-1),'contiguous') == 1
                alpha_star(i) = trirnd(max(theta(i).age),alpha(i),min(theta(i-1).age),1);
                q_alpha_num(i) = tripdf(alpha(i),max(theta(i).age),alpha_star(i),min(theta(i-1).age));
                q_alpha_den(i) = tripdf(alpha_star(i),max(theta(i).age),alpha(i),min(theta(i-1).age));
            elseif strcmpi(boundary_relationship(i-1),'disjoint') == 1
                x = (beta(i-1) + alpha(i))/2;
                x = (x + alpha(i))/2;
                alpha_star(i) = trirnd(max(theta(i).age),alpha(i),x,1);
                q_alpha_num(i) = tripdf(alpha(i),max(theta(i).age),alpha_star(i),x);
                q_alpha_den(i) = tripdf(alpha_star(i),max(theta(i).age),alpha(i),x);
            end
        end 
    end
    alpha_star(i) = round(alpha_star(i));
end   
%% propose a move for beta
for i = 1:M
    if i < M % first M-1 periods
        if strcmpi(year_scale,'BCE/CE') == 1
            if strcmpi(boundary_relationship(i),'overlapping') == 1
                id = theta(i+1).age > beta(i);
                x = min(theta(i+1).age(id));
                beta_star(i) = trirnd(max(theta(i).age),beta(i),x,1);
                q_beta_num(i) = tripdf(beta(i),max(theta(i).age),beta_star(i),x);
                q_beta_den(i) = tripdf(beta_star(i),max(theta(i).age),beta(i),x);
            elseif strcmpi(boundary_relationship(i),'contiguous') == 1
                beta_star(i) = alpha_star(i+1);
                q_beta_num(i) = tripdf(beta(i),max(theta(i).age),beta_star(i),min(theta(i+1).age));
                q_beta_den(i) = tripdf(beta_star(i),max(theta(i).age),beta(i),min(theta(i+1).age));
            elseif strcmpi(boundary_relationship(i),'disjoint') == 1
                x = (beta(i) + alpha(i+1))/2;
                x = (beta(i) + x)/2;
                beta_star(i) = trirnd(max(theta(i).age),beta(i),x,1);
                q_beta_num(i) = tripdf(beta(i),max(theta(i).age),beta_star(i),x);
                q_beta_den(i) = tripdf(beta_star(i),max(theta(i).age),beta(i),x);
            end
        elseif strcmpi(year_scale,'BP') == 1
            if strcmpi(boundary_relationship(i),'overlapping') == 1
                id = theta(i+1).age < beta(i);
                x = max(theta(i+1).age(id));
                beta_star(i) = trirnd(x,beta(i),min(theta(i).age),1);
                q_beta_num(i) = tripdf(beta(i),x,beta_star(i),min(theta(i).age));
                q_beta_den(i) = tripdf(beta_star(i),x,beta(i),min(theta(i).age));
            elseif strcmpi(boundary_relationship(i),'contiguous') == 1
                beta_star(i) = alpha_star(i+1);
                q_beta_num(i) = tripdf(beta(i),max(theta(i+1).age),beta_star(i),min(theta(i).age));
                q_beta_den(i) = tripdf(beta_star(i),max(theta(i+1).age),beta(i),min(theta(i).age));
            elseif strcmpi(boundary_relationship(i),'disjoint') == 1
                x = (beta(i) + alpha(i+1))/2;
                x = (beta(i) + x)/2;
                beta_star(i) = trirnd(x,beta(i),min(theta(i).age),1);
                q_beta_num(i) = tripdf(beta(i),x,beta_star(i),min(theta(i).age));
                q_beta_den(i) = tripdf(beta_star(i),x,beta(i),min(theta(i).age));
            end
        end 
    else % the last period
        if strcmpi(year_scale,'BCE/CE') == 1
            x = (beta(i) + B)/2;
            beta_star(i) = trirnd(max(theta(i).age),beta(i),x,1);
            q_beta_num(i) = tripdf(beta(i),max(theta(i).age),beta_star(i),x);
            q_beta_den(i) = tripdf(beta_star(i),max(theta(i).age),beta(i),x);
        elseif strcmpi(year_scale,'BP') == 1
            x = (beta(i) + B)/2;
            beta_star(i) = trirnd(x,beta(i),min(theta(i).age),1);
            q_beta_num(i) = tripdf(beta(i),x,beta_star(i),min(theta(i).age));
            q_beta_den(i) = tripdf(beta_star(i),x,beta(i),min(theta(i).age));
        end
    end
    beta_star(i) = round(beta_star(i));
end
%% calulate the acceptance probability
num1 = boundary_prior(alpha_star,beta_star,A,B,boundary_relationship,year_scale);
num2 = sum(log(q_alpha_num)) + sum(log(q_beta_num));
den1 = boundary_prior(alpha,beta,A,B,boundary_relationship,year_scale);
den2 = sum(log(q_alpha_den)) + sum(log(q_beta_den));
omega = (num1+num2) - (den1+den2);
omega = min(omega,0);
rho = log(rand(1,1));
if omega > rho
    alpha_out = alpha_star;
    beta_out = beta_star;
else
    alpha_out = alpha;
    beta_out = beta;
end
end