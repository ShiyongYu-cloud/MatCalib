function p = theta_prior(theta,alpha,beta,P,year_scale)
%% function for calculating the prior probability of theta in log scale
%INPUT
%theta: structure containing calendar ages of the sequence
%alpha: vector containing calendar ages of early boundary of all periods
%beta: vector containing calendar ages of late boundary of all periods
%P: structure containing laboratory data of all periods
%year_scale: scale of year (e.g. BCE/CE or BP)
%OUTPUT 
%p: probability in log scale
%%
M = length(theta);%number of periods
N = zeros(1,M);
for i = 1:M
    N(i) = length(theta(i).age); %number of ages in each period
end
F = zeros(1,M);
I = zeros(1,M);
for i = 1:M
     if strcmpi(year_scale,'BCE/CE') == 1
         if all(theta(i).age > alpha(i)) && all(theta(i).age < beta(i))
             F(i) = power(beta(i)-alpha(i),-N(i));
         else 
             F(i) = 0;
         end
     elseif strcmpi(year_scale,'BP') == 1
         if all(theta(i).age < alpha(i)) && all(theta(i).age > beta(i))
             F(i) = power(alpha(i)-beta(i),-N(i));
         else 
             F(i) = 0;
         end  
     end
    if strcmpi(P(i).ordering,'unordered') == 1 % the "unordered" case
         I(i) = 1; % odering does not matter
    elseif strcmpi(P(i).ordering,'coeval') == 1 % the "coeval" case
         if all(diff(theta(i).age) == 0)
             I(i) = 1;
         else
             I(i) = 0;
         end
    elseif strcmpi(P(i).ordering,'ordered') == 1 %the "ordered" case
         if strcmpi(year_scale,'BCE/CE') == 1
             if all(diff(theta(i).age) > 0)
                 I(i) = 1;
             else
                 I(i) = 0;
             end
         elseif strcmpi(year_scale,'BP') == 1
             if all(diff(theta(i).age) < 0)
                 I(i) = 1;
             else
                 I(i) = 0;
             end       
         end
    end
end    
p = sum(log(F)) + sum(log(I));
end