function p = likelihood(P,theta,A,B,year_scale,CalCurves)
%% function for calculating the likelihood of radiocabon ages given calendar ages in F14C space
%INPUT
%P: structure containing the F14C values and erros of radiocarbon ages in the sequence
%theta: structure containing the calendar ages of the periods
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%year_scale: scale of year ("BCE/CE" or "BP")
%CalCurves: sturcture containing the calibration curves
%OUTPUT
%p: likelihood value in log scale
%% find the radiocaron age of theta from the calibration curve and convert it to F14C
M = length(P); % number of periods
N = zeros(1,M);
F = zeros(1,M);
for i = 1:M
    N(i) = length(P(i).age); 
end 
F_hat = struct('C14',[],'err',[]);
for i = 1:M
    if strcmpi(year_scale,'BCE/CE') == 1
        if any(theta(i).age < A) || any(theta(i).age > B) 
            F(i) = log(0);
        else
            [F_hat(i).C14, F_hat(i).err] = calendar2f(theta(i).age,P(i).cal_curve,CalCurves,year_scale);
            a = (P(i).f14C - F_hat(i).C14).^2;
            b = 2*(P(i).ferr.^2 + F_hat(i).err.^2);
            c = gamma(N(i)/2);
            F(i) = log(c)-(N(i)/2)*log(sum(a./b));
            %c = sqrt(2*pi*(P(i).ferr.^2 + F_hat(i).err.^2));
            %F(i) = sum(-log(c))+sum(-a./b);
        end   
    elseif strcmpi(year_scale,'BP') == 1
        if any(theta(i).age > A) || any(theta(i).age < B) 
            F(i) = log(0);
        else
            [F_hat(i).C14, F_hat(i).err] = calendar2f(theta(i).age,P(i).cal_curve,CalCurves,year_scale);
            a = (P(i).f14C - F_hat(i).C14).^2;
            b = 2*(P(i).ferr.^2 + F_hat(i).err.^2);
            c = gamma(N(i)/2);
            F(i) = log(c)-(N(i)/2)*log(sum(a./b));
            %c = sqrt(2*pi*(P(i).ferr.^2 + F_hat(i).err.^2));
            %F(i) = sum(-log(c))+sum(-a./b);
        end
    end
end
p = sum(F);
return;
%%
function [F14_val, F14_err] = calendar2f(cal_age,cal_curves,CalCurves,year_scale)
%% function for mapping a calendar age to radiocarbon age and converting it to F14C
%INPUT
%cal_age: calendar ages in a period
%cal_curves: calibration curves associated with the corresponding 14C ages
%CalCurves: structure containing the extracted calibration curves
%year_scale: scale of years to be reported (BCE/CE or BP)
%OUTPUT
%F14_val: F14C values in the calibration curve
%F14_err: F14C errors in the calibration curve  
%% convert to BP scale if needed
if strcmpi(year_scale,'BCE/CE') == 1
   cal_age = 1950 - cal_age; 
end
%% find the corresponding F_14C values in the calibration curve
N = length(cal_age);
F14_val = zeros(1,N);
F14_err = zeros(1,N);
for i = 1:N
    % retrieve the calibration dataset
    if strcmpi(cal_curves(i),'intcal13') == 1
        curve_cal_age = CalCurves.intcal13(:,1);
        curve_F14_val = CalCurves.intcal13(:,2);
        curve_F14_err = CalCurves.intcal13(:,3);
        F14_val(i) = interp1(curve_cal_age, curve_F14_val, cal_age(i));
        F14_err(i) = interp1(curve_cal_age, curve_F14_err, cal_age(i));
    elseif strcmpi(cal_curves(i),'marine13') == 1
        curve_cal_age = CalCurves.marine13(:,1);
        curve_F14_val = CalCurves.marine13(:,2);
        curve_F14_err = CalCurves.marine13(:,3);
        F14_val(i) = interp1(curve_cal_age, curve_F14_val, cal_age(i));
        F14_err(i) = interp1(curve_cal_age, curve_F14_err, cal_age(i));
    elseif strcmpi(cal_curves(i),'shcal13') == 1
        curve_cal_age = CalCurves.shcal13(:,1);
        curve_F14_val = CalCurves.shcal13(:,2);
        curve_F14_err = CalCurves.shcal13(:,3);   
        F14_val(i) = interp1(curve_cal_age, curve_F14_val, cal_age(i));
        F14_err(i) = interp1(curve_cal_age, curve_F14_err, cal_age(i));
    elseif strcmpi(cal_curves(i),'intcal20') == 1
        curve_cal_age = CalCurves.intcal20(:,1);
        curve_F14_val = CalCurves.intcal20(:,2);
        curve_F14_err = CalCurves.intcal20(:,3);
        F14_val(i) = interp1(curve_cal_age, curve_F14_val, cal_age(i));
        F14_err(i) = interp1(curve_cal_age, curve_F14_err, cal_age(i));
    elseif strcmpi(cal_curves(i),'marine20') == 1
        curve_cal_age = CalCurves.marine20(:,1);
        curve_F14_val = CalCurves.marine20(:,2);
        curve_F14_err = CalCurves.marine20(:,3);
        F14_val(i) = interp1(curve_cal_age, curve_F14_val, cal_age(i));
        F14_err(i) = interp1(curve_cal_age, curve_F14_err, cal_age(i));
    elseif strcmpi(cal_curves(i),'shcal20') == 1
        curve_cal_age = CalCurves.shcal20(:,1);
        curve_F14_val = CalCurves.shcal20(:,2);
        curve_F14_err = CalCurves.shcal20(:,3);        
        F14_val(i) = interp1(curve_cal_age, curve_F14_val, cal_age(i));
        F14_err(i) = interp1(curve_cal_age, curve_F14_err, cal_age(i));
    end
end
return;