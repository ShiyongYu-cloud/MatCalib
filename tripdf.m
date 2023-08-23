function f = tripdf(x,a,c,b)
%% function for calculating the pdf of the triangular distribution supported on [a,b] with shape parameter c
%INPUT
%x: random number
%a: lower bound of random variable X
%c: Peak parameter
%b: upper bound of random variable X
%OUTPUT
%f: probability of x
%%
if x < a
   f = 0; 
elseif x >= a && x < c
   f = 2*(x-a)/((b-a)*(c-a));
elseif x == c 
   f = 2/(b-a); 
elseif x > c && x <= b
   f = 2*(b-x)/((b-a)*(b-c));
elseif x > b
   f = 0;     
end
end