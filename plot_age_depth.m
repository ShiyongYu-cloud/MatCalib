function [] = plot_age_depth(pdfs,CAGES,P,alpha_depth,beta_depth,depth_unit,year_scale,A,B)
%% function for plotting the modeled ages against depth
%INPUT
%pdfs: estimated pdf of parameters
%CAGES: structure wrapping up ALPHA, BETA, and THETA
%P: structure containing laboratory data of all periods
%alpha_depth: vector containing depths of early boundaries of all periods
%beta_depth: vector containing depths of late boundaries of all periods
%year_scale: scale of year (BCE/CE or BP)
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%% extract depth info and assemble them in a sequence
M = length(P);
theta_depth = struct('depth',[]);
for i = 1:M
    theta_depth(i).depth = P(i).depth;
end 
S_depth = P2S_depth(theta_depth,alpha_depth,beta_depth);
%% ploting ages against depths
cal_age = pdfs(:,1); 
prob = pdfs(:,2:end);
[M,N] = size(prob);
prob_u = zeros(M,N);
dy = abs(diff(S_depth));
for i = 1:N-1
    if dy(i) == 0
        dy(i) = dy(i+1);
    end
end
DY = min(dy)*ones(1,N);
for i = 1:N
    factor = 0.8*DY(i)/max(prob(:,i));
    prob_u(:,i) = S_depth(i) + prob(:,i)*factor;  % blow up and shift
    ind1 = (prob_u(:,i) > S_depth(i)+0.000001);
    if strcmpi(year_scale,'BCE/CE') == 1 
        early = min(cal_age(ind1));
        late = max(cal_age(ind1)); 
        ind2 = (cal_age >= early) & (cal_age <= late);
    elseif strcmpi(year_scale,'BP') == 1
        early = max(cal_age(ind1));
        late = min(cal_age(ind1)); 
        ind2 = (cal_age < early) & (cal_age > late);
    end        
    calage = cal_age(ind2);
    CAL_AGE = [calage(1); calage; calage(end); calage(1)];
    PROB = [S_depth(i); prob_u(ind2,i); S_depth(i); S_depth(i)];
    fill(CAL_AGE,PROB,[0.301 0.745 0.933]);  % plot pdfs
    hold on
    p95_4 = CAGES(i).P95_4_regions;    % plot the 95.4% pdf
    K = size(p95_4,1);
    for j = 1:K
        id = cal_age >= p95_4(j,2) & cal_age <= p95_4(j,1);
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf = prob_u(id,i);
        PDF = [S_depth(i); pdf; S_depth(i); S_depth(i)];
        fill(AGE,PDF,[0.929 0.694 0.125]);
    end
    hold on
    p68_2 = CAGES(i).P68_2_regions;    % plot the 68.2% pdf 
    L = size(p68_2,1);
    for k = 1:L
        id = cal_age >= p68_2(k,2) & cal_age <= p68_2(k,1);
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf = prob_u(id,i);
        PDF = [S_depth(i); pdf; S_depth(i); S_depth(i)];
        fill(AGE,PDF,[0.635 0.078 0.184]); 
    end
end
ylim([S_depth(1)-DY(1) S_depth(end)+DY(end)])
ylabel(strcat('Depth (', depth_unit,')'));
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
set(gca,'YMinorTick','on');
legend('Posterior probability','95.4 % HPD regions','68.2 % HPD regions','location', 'southeast'); 

% plot sampleID 
for i = 1:N
    if strcmpi(year_scale,'BP') == 1 
       if i > 1 && (S_depth(i) == S_depth(i-1))
           text(A-10,S_depth(i)+DY(i)/4,strcat('/',CAGES(i).SampleID),'FontSize',8); 
       else
           text(A-10,S_depth(i)+DY(i)/4,CAGES(i).SampleID,'FontSize',8);
       end    
    elseif strcmpi(year_scale,'BCE/CE') == 1 
       if i > 1 && (S_depth(i) == S_depth(i-1))
           text(A+10,S_depth(i)+DY(i)/4,strcat('/',CAGES(i).SampleID)); 
       else
           text(A+10,S_depth(i)+DY(i)/4,CAGES(i).SampleID);           
       end    
    end
end
return;
%
function S_depth = P2S_depth(theta_depth,alpha_depth,beta_depth)
%% function for converting a structure to a vectori of depths
%INPUT
%theta_depth: structure containing the depths of calendar ages of all periods 
%alpha_depth: vector containing the depths of ealy boundary of all periods
%beta_depth: vector containing the depth of late boundary of alll periods
%OUTPUT
%S_depth: vector containing depths of all model paramters
%%
M = length(theta_depth); %number of periods
N_theta = zeros(1,M);   
for i = 1:M
    N_theta(i) = length(theta_depth(i).depth);
end
S_depth = zeros(1,2*M+sum(N_theta)); %set up the dimention of the age sequence
ind_alpha = cumsum([1 N_theta(1:end-1)+2]); %index of alpha in the age sequence
ind_b_theta = cumsum([1 N_theta(1:end-1)+2]) + 1; %index of beginning of theta in the sequence 
ind_e_theta = cumsum([1 N_theta(1:end-1)+2]) + N_theta; %index of end of theta in the sequence
ind_beta = cumsum([1 N_theta(1:end-1)+2]) + N_theta + 1;  %index of beta in the age sequence
for i = 1:M
    S_depth(ind_alpha(i)) = alpha_depth(i);
    S_depth(ind_b_theta(i):ind_e_theta(i)) = theta_depth(i).depth;
    S_depth(ind_beta(i)) = beta_depth(i);
end
return;