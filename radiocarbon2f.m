function [F_14C,F_err] = radiocarbon2f(R_age,R_err,R_rage,R_rerr)
%% function for converting radiocarbon ages to the F14C space
%INPUT
%R_age: Laboratory radiocarbon ages
%R_err: 1 sigma standard deviation of radiocarbon ages
%R_rage: local reservoir age or reservoir age offset
%R_rerr: 1 sigma standard deviation of local reservoir ages or age offsets
%OUTPUT
%F_14C: 14C age in the F14C space
%F_err: one standard deviaton of F_14C 
%function usage example: [F_14C,F_err] = radiocarbon2f(3000,20,150,30)
%% correct for the reservoir effect
R_age = R_age - R_rage;
R_err = sqrt(R_err.^2 + R_rerr.^2);
%% convert to F14C space
F_14C = exp(R_age/-8033);
F_err = F_14C.*R_err/8033;
end