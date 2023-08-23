function [] = plot_pooled_mean(pooled_pdfs,POOLED_AGE,year_scale)
%% function for plotting the pdf and hpd regions of the pooled mean ages
%INPUT
%pooled_pdfs: matrix containing pooled pdf and its higher posterior probabilities
%POOLED_AGE: structure wrapping up the pooled mean age and hpd regions
%year_scale: scale of year (BCE/CE or BP)
%% plot pdf
cal_age = pooled_pdfs(:,1); 
prob = pooled_pdfs(:,2);
ind1 = prob > 0.000001;
if strcmpi(year_scale,'BCE/CE') == 1 
    early = min(cal_age(ind1));
    late = max(cal_age(ind1)); 
    ind2 = (cal_age >= early) & (cal_age <= late);
elseif strcmpi(year_scale,'BP') == 1
    early = max(cal_age(ind1));
    late = min(cal_age(ind1)); 
    ind2 = (cal_age <= early) & (cal_age >= late);
end        
calage = cal_age(ind2);
CAL_AGE = [calage(1); calage; calage(end); calage(1)];
PROB = [0; prob(ind2); 0; 0];
fill(CAL_AGE,PROB,[0.301 0.745 0.933]);  
hold on
%% plot the 95.4% hpd regions
p95_4 = POOLED_AGE.P95_4_regions;    
K = size(p95_4,1);
for j = 1:K
    id = cal_age >= p95_4(j,2) & cal_age <= p95_4(j,1);
    age = cal_age(id);
    AGE = [age(1); age; age(end); age(1)];
    pdf = prob(id);
    PDF = [0; pdf; 0; 0];
    fill(AGE,PDF,[0.929 0.694 0.125]);
end
hold on
%% plot the 68.2% hpd regions 
p68_2 = POOLED_AGE.P68_2_regions;    
L = size(p68_2,1);
for k = 1:L
    id = cal_age >= p68_2(k,2) & cal_age <= p68_2(k,1);
    age = cal_age(id);
    AGE = [age(1); age; age(end); age(1)];
    pdf = prob(id);
    PDF = [0; pdf; 0; 0];
    fill(AGE,PDF,[0.635 0.078 0.184]); 
end
set(gca,'yticklabel',[]);
set(gca,'ytick',[]);
grid on;
if strcmpi(year_scale,'BCE/CE') == 1
   xlabel('Pooled mean calendar age (BCE/CE)');  
elseif strcmpi(year_scale,'BP') == 1
   set(gca,'XDir','reverse');
   xlabel('Pooled mean calendar age (BP)');
end
set(gca, 'TickDir', 'out');
set(gca,'XMinorTick','on');
legend('Posterior probability','95.4 % HPD regions','68.2 % HPD regions','location', 'northwest'); 
grid on;