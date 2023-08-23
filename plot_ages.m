function [] = plot_ages(pdfs,CAGES,year_scale,A,B)
%% function for plotting the results of model parameters
%INPUT
%pdfs: estimated pdf of parameters
%CAGES: structure wrapping up ALPHA, BETA, and THETA
%year_scale: scale of year (BCE/CE or BP)
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%% plot pdfs
cal_age = pdfs(:,1); 
prob = pdfs(:,2:end);
N = size(prob,2);
maxprob = max(max(prob))*N/(N-1);
for i = 1:N
    prob(:,i) = prob(:,i) + maxprob*(i-1);  % blow up and shift
    ind1 = (prob(:,i) > maxprob*(i-1)+0.000001);
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
    PROB = [maxprob*(i-1); prob(ind2,i); maxprob*(i-1); maxprob*(i-1)];
    fill(CAL_AGE,PROB,[0.301 0.745 0.933]);  
    hold on
    %% plot the 95.4% hpd regions
    p95_4 = CAGES(i).P95_4_regions;    
    K = size(p95_4,1);
    for j = 1:K
        id = cal_age >= p95_4(j,2) & cal_age <= p95_4(j,1);
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf = prob(id,i);
        PDF = [maxprob*(i-1); pdf; maxprob*(i-1); maxprob*(i-1)];
        fill(AGE,PDF,[0.929 0.694 0.125]);
    end
    hold on
    %% plot the 68.2% hpd regions 
    p68_2 = CAGES(i).P68_2_regions;    
    L = size(p68_2,1);
    for k = 1:L
        id = cal_age >= p68_2(k,2) & cal_age <= p68_2(k,1);
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf = prob(id,i);
        PDF = [maxprob*(i-1); pdf; maxprob*(i-1); maxprob*(i-1)];
        fill(AGE,PDF,[0.635 0.078 0.184]); 
    end
end
ylim([-maxprob N*maxprob]);
set(gca,'yticklabel',[]);
set(gca,'ytick',[]);
grid on;
if strcmpi(year_scale,'BCE/CE') == 1
   A = A-(B-A)/4; 
   xlim([A B]);
   xlabel('Calendar age (BCE/CE)');  
elseif strcmpi(year_scale,'BP') == 1
   A = A+(A-B)/4; 
   xlim([B A]);
   set(gca,'XDir','reverse');
   xlabel('Calendar age (BP)');
end
set(gca, 'TickDir', 'out');
set(gca,'XMinorTick','on');
legend('Posterior probability','95.4 % HPD regions','68.2 % HPD regions','location', 'southeast'); 
grid on;
%% plot sampleID and radiocarbon ages 
for i = 1:N
    if strcmpi(year_scale,'BP') == 1 
       text(A-10,maxprob*0.25+maxprob*(i-1),CAGES(i).SampleID,'FontSize',8); 
    elseif strcmpi(year_scale,'BCE/CE') == 1 
       text(A+10,maxprob*0.25+maxprob*(i-1),CAGES(i).SampleID,'FontSize',8); 
    end
end
end