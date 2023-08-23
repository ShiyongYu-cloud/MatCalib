function p = boundary_prior(alpha,beta,A,B,boundary_relationship,year_scale)
%% function for calculating the prior probability of period boundaries in log scale
%INPUT
%alpha: vector containing the calendar age of early bundary of all periods
%beta: vector containing the calendar age of late bundary of all periods
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%boundary_relationship: relationship of two neighboring periods (overlapping|contiguous|disjoint)
%year_scale: scale of age (BCE/CE or BP)
%OUTPUT
%p: prabability in log scale
%%
M = length(alpha); %number of periods
F = zeros(1,M);
I = zeros(1,M);
%%
for i = 1:M
    if i == 1 || M % the first or last period
        if strcmpi(year_scale,'BCE/CE') == 1
            if A < alpha(i) && alpha(i) < B && A < beta(i) && beta(i) < B 
                F(i) = power(B-A,-2);
            else
                F(i) = 0;
            end
            if alpha(i) < beta(i)
                I(i) = 1;
            else
                I(i) = 0;
            end    
        elseif strcmpi(year_scale,'BP') == 1
            if A > alpha(i) && alpha(i) > B && A > beta(i) && beta(i) > B
                F(i) = power(A-B,-2);
            else
                F(i) = 0;
            end 
            if alpha(i) > beta(i)
                I(i) = 1;
            else
                I(i) = 0;
            end    
        end    
    else  % internal periods
        A = alpha(1);
        B = beta(M);
        if strcmpi(year_scale,'BCE/CE') == 1
            if A < alpha(i) && alpha(i) < B && A < beta(i) && beta(i) < B
                F(i) = power(B-A,-2);
            else
                F(i) = 0;
            end
            if strcmpi(boundary_relationship{i-1},'overlapping') == 1
                if alpha(i-1) < alpha(i) && alpha(i) < beta(i-1) && beta(i-1) < beta(i)
                    I(i) = 1;
                else
                    I(i) = 0;
                end            
            elseif strcmpi(boundary_relationship{i-1},'contiguous') == 1
                if alpha(i-1) < beta(i-1) && beta(i-1) == alpha(i) && alpha(i) < beta(i)
                    I(i) = 1;
                else
                    I(i) = 0;
                end            
            elseif strcmpi(boundary_relationship{i-1},'disjoint') == 1
                if alpha(i-1) < beta(i-1) && beta(i-1) < alpha(i) && alpha(i) < beta(i)
                    I(i) = 1;
                else
                    I(i) = 0;
                end
            end
        elseif strcmpi(year_scale,'BP') == 1
            if A > alpha(i) && alpha(i) > B && A > beta(i) && beta(i) > B
                F(i) = power(A-B,-2);
            else
                F(i) = 0;
            end
            if strcmpi(boundary_relationship{i-1},'overlapping') == 1
                if alpha(i-1) > alpha(i) && alpha(i) > beta(i-1) && beta(i-1) > beta(i)
                    I(i) = 1;
                else
                    I(i) = 0;
                end            
            elseif strcmpi(boundary_relationship{i-1},'contiguous') == 1
                if alpha(i-1) > beta(i-1) && beta(i-1) == alpha(i) && alpha(i) > beta(i)
                    I(i) = 1;
                else
                    I(i) = 0;
                end            
            elseif strcmpi(boundary_relationship{i-1},'disjoint') == 1
                if alpha(i-1) > beta(i-1) && beta(i-1) > alpha(i) && alpha(i) > beta(i)
                    I(i) = 1;
                else
                    I(i) = 0;
                end
            end
        end            
    end
end 
p = sum(log(F) + log(I));    
end