function [] = savemcmc(filename,P,samples)
%% function for saving samples from the MCMC run to a file
%INPUT
%filename: name of file
%R: information of radiocarbon ages assembled in a structure
%samples: mcmcsamples 
%% creat a header
M = length(P);      %number of periods
N = zeros(1,M);     %number of ages in each period
for i = 1:M
    N(i) = length(P(i).age); 
end
header = cell(1,2*M+sum(N)); 
ind_alpha = cumsum([1 N(1:end-1)+2]); %index of alpha in the age sequence
ind_b_theta = cumsum([1 N(1:end-1)+2]) + 1; %index of beginning of theta in the sequence 
%ind_e_theta = cumsum([1 N(1:end-1)+2]) + N; %index of end of theta in the sequence
ind_beta = cumsum([1 N(1:end-1)+2]) + N + 1;  %index of beta in the age sequence
for i = 1:M
    header{1,ind_alpha(i)} = strcat('Phase_',num2str(i),' early boundary');
    header{1,ind_beta(i)} = strcat('Phase_',num2str(i),' late boundary');
    for j = 1:N(i)
        header{1,ind_b_theta(i)+j-1} = P(i).lab_code{j};
    end
end
%% write header and data to the file
fmt = repmat('%25d ',1,size(samples,2)); 
fmt = [fmt,'\n'];
fid = fopen(filename,'wt');
fprintf(fid,'%25s ',header{1:end});
fprintf(fid,'\n');
for k = 1:size(samples,1)
    fprintf(fid,fmt,samples(k,:));
end    
fclose(fid);
end