function P = read_data(DataFile)
%% function for reading radiocarbon age info from a data file
%INPUT
% datafile: name of file containing radiocarbon age data
%OUTPUT
% P: a struct containing the data
%%
h = fopen(DataFile);
data = textscan(h,'%d %s %f %f %f %f %f %s %s','headerlines',1,'delimiter','\t');
fclose(h);
period_number = data{1};
lab_ID = data{2};
depth = data{3};
age = data{4};
err = data{5};
rage = data{6};
rerr = data{7};
calcurve = data{8}; 
temporal_order = data{9};
M = max(unique(period_number)); % number of periods in the sequence
P = struct('lab_code',[],'depth',[],'age',[],'err',[],'rage',[],'rerr',[],'cal_curve',[],'ordering',[]);
for i = 1:M
    ind = period_number == i;
    P(i).lab_code = lab_ID(ind)';
    P(i).depth = depth(ind)';
    P(i).age = age(ind)';
    P(i).err = err(ind)';
    P(i).rage = rage(ind)';
    P(i).rerr = rerr(ind)';
    P(i).cal_curve = calcurve(ind)';
    P(i).ordering = unique(temporal_order(ind));
end    
end