function [A,B] = sequence_bound(P,year_scale)
%% function for estimating the early and late boundary of the age sequence
%INPUT
%P: structure containing laboratory data of all periods of a sequence
%year_scale: scale of year to be reported (e.g. BCE/CE or BP)
%OUTPUT
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%% find the max of the laboratory radiocarbon ages in the earliest period
[max_age,max_ind] = max(P(1).age);
max_err = P(1).err(max_ind);
%% find the min of the laboratory radiocarbon ages in the latest period
[min_age,min_ind] = min(P(end).age);
min_err = P(end).err(min_ind);
%% find the reservoir age and error of the max and min radiocarbin ages
max_rage = P(1).rage(max_ind);
min_rage = P(end).rage(min_ind);
max_rerr = P(1).rerr(max_ind);
min_rerr = P(end).rerr(min_ind);
%% find the calibration curves applied on the max and min radiocarbon ages
max_curve = P(1).cal_curve(max_ind);
min_curve = P(end).cal_curve(min_ind);
%% correct for the reservoir effect, if any, and convert to the F14C space
[max_F_14C,max_F_err] = radiocarbon2f(max_age,max_err,max_rage,max_rerr);
[min_F_14C,min_F_err] = radiocarbon2f(min_age,min_err,min_rage,min_rerr);
%% find the max and min of the calibrated radiocarbon ages
[max_cal_age,max_prob] = calibration(max_F_14C,max_F_err,year_scale,max_curve{1});
[min_cal_age,min_prob] = calibration(min_F_14C,min_F_err,year_scale,min_curve{1});
%% round up A and B to the nearest 10 years; 
if strcmpi(year_scale,'BCE/CE') == 1
    A = min(max_cal_age(max_prob>0.00001));
    B = max(min_cal_age(min_prob>0.00001));
    A = round(A/10)*10;
    B = round(B/10)*10;
    if B > 1950
        B = 1950;
    end  
elseif strcmpi(year_scale,'BP') == 1   
    A = max(max_cal_age(max_prob>0.00001));
    B = min(min_cal_age(min_prob>0.00001));
    A = round(A/10)*10;
    B = round(B/10)*10;
    if B < 0
       B = 0;
    end   
end    
return;
%%
function [cal_age,prob] = calibration(F_14C_val,F_14C_err,year_scale,cal_curve)
%% Function for calibrating a 14C age using Bayesian higher posterior density analysis
%INPUT
% F_14C:    Lab 14C age in F 14C space
% F_err:    Lab 14C uncertainty in F space
% year_scale: scale of year to be reported (e.g. BCE/CE or BP)
% cal_curve: name of the calibration curve
%OUTPUT
% cal_age:  calibrated age in year_scale (BC/AD or BP)
% prob: probability of the calibrated ages
%%
%% open the calibration curve data file
headerlines = 11;
h = fopen([cal_curve,'.14c']);
cal_dataset = textscan(h,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(h);
curve_cal_age = flipud(cal_dataset{1});
curve_C14_age = flipud(cal_dataset{2});
curve_C14_err = flipud(cal_dataset{3});
%% convert radiocrbon ages in the calibration curve to F_14C space
curve_F14_val = exp(curve_C14_age/-8033); %ages to F_14C
curve_F14_err = curve_F14_val.*curve_C14_err/8033; %errors to F_14C
% interpolate the calendar age in the curve to annual resolution
interpres = 1;
hi_curve_cal_age = curve_cal_age(1):interpres:curve_cal_age(end);
hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
%% Calculate probability for every calendar year in the F14 space
cal_age = hi_curve_cal_age;
a = (F_14C_val - hi_curve_F14_val).^2;
b = 2*(F_14C_err^2 + hi_curve_F14_err.^2);
c = sqrt(F_14C_err^2 + hi_curve_F14_err.^2);
prob = exp(-a./b)./c;
%% normalize the probabilities to 1
prob = prob/sum(prob);
%% convert to BC/AD scale if needed
if strcmpi(year_scale,'BCE/CE') == 1
   cal_age = 1950 - cal_age;
end
return;