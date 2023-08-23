function CalCurves = read_curves(P,A,B,year_scale)
%% function for reading calibration curves only needed for age modeling
%INPUT:
%P: structure containing laboratory data of all periods of a sequence
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%year_scale: scale of year to be reported (e.g. BCE/CE or BP)
%OUTPUT
%CalCurve: sturcture containing the calibration curves
%%
M = length(P);
N = zeros(1,M);            
for i = 1:M
    N(i) = length(P(i).cal_curve); 
end
%% put calibration curve into a cell array
curves = cell(1,sum(N));
b = cumsum([1 N(1:end-1)]);
e = cumsum([1 N(1:end-1)])+N-1;
for i = 1:M
    curves(b(i):e(i)) = P(i).cal_curve;
end    
%% extract unique curves
curves = unique(curves);
%% read and extract the curves and out them into a structure
for i = 1:length(curves)
    if any(strcmpi(curves(i),'intcal13') == 1)
        CalCurves.intcal13 = extract_curve('intcal13',A,B,year_scale);
    elseif any(strcmpi(curves(i),'marine13') == 1)
        CalCurves.marine13 = extract_curve('marine13',A,B,year_scale);
    elseif any(strcmpi(curves(i),'shcal13') == 1)
        CalCurves.shcal13 = extract_curve('shcal13',A,B,year_scale);
    elseif any(strcmpi(curves(i),'intcal20') == 1)
        CalCurves.intcal20 = extract_curve('intcal20',A,B,year_scale);
    elseif any(strcmpi(curves(i),'marine20') == 1)
        CalCurves.marine20 = extract_curve('marine20',A,B,year_scale);
    elseif any(strcmpi(curves(i),'shcal20') == 1)
        CalCurves.shcal20 = extract_curve('shcal20',A,B,year_scale);
    end
end
return;
%%
function cal_dataset = extract_curve(cal_curve,A,B,year_scale)
%% function for reading the calibration curve and converting it to the F14C space
%INPUT
%cal_curve: name of the calibration curve
%A: early bound of the study period
%B: late bound of the study period
%year_scale: scale of year to be reported (e.g. BCE/CE or BP)
%OUTPUT
%cal_dataset: array containing the calibration curve bounded on [A.B]
%% open the calibration curve data
headerlines = 11;
h = fopen([cal_curve,'.14c']);
cal_dataset = textscan(h,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(h);
curve_cal_age = flipud(cal_dataset{1});
curve_C14_age = flipud(cal_dataset{2});
curve_C14_err = flipud(cal_dataset{3});
%% convert radiocrbon ages in the calibration curve to F_14C space
curve_F14_val = exp(curve_C14_age/-8033); %convert the radiocarbon ages to the F_14C space
curve_F14_err = curve_F14_val.*curve_C14_err/8033; %convert the error to the F_14C space
%% extract data
if strcmpi(year_scale,'BCE/CE') == 1
    A = 1950 - A;
    B = 1950 - B;
end
ind = (curve_cal_age >= B & curve_cal_age <= A);
cal_dataset = [curve_cal_age(ind) curve_F14_val(ind) curve_F14_err(ind)];
return;