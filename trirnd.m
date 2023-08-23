function x = trirnd(a,c,b,n)
%% function for generating n random numbers from the triangular distribution on[a,b] with shape parameter c 
%INPUT
%a: lower bound of random variable X
%c: peak parameter
%b: upper bound of random variable X
%n: number of samples to be drawn 
%OUTPUT
%X: random number of X
%%
x = zeros(n,1);
for i = 1:n
    z = rand(1,1);
    if (sqrt(z*(b-a)*(c-a)) + a) < c
        x(i) = sqrt(z*(b-a)*(c-a)) + a;
    else
        x(i) = b - sqrt((1-z)*(b-a)*(b-c));
    end
end 
end