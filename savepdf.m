function [] = savepdf(filename,P,year_scale,pdfs)
%% function for saving the empirical pdfs to a file
%INPUT
%filename: file name to which for data to be saved
%P: information of radiocarbon ages assembled in a structure
%year_scale: scale of year (BCE/CE or BP)
%pdfs: probabilities of the calendar ages
%% creat a header
M = length(P);      %number of periods
N = zeros(1,M);     %number of ages in each period
for i = 1:M
    N(i) = length(P(i).age); 
end
header = cell(1,2*M+sum(N)+1); 
ind_alpha = cumsum([1 N(1:end-1)+2]); %index of alpha in the age sequence
ind_b_theta = cumsum([1 N(1:end-1)+2]) + 1; %index of beginning of theta in the sequence 
%ind_e_theta = cumsum([1 N(1:end-1)+2]) + N; %index of end of theta in the sequence
ind_beta = cumsum([1 N(1:end-1)+2]) + N + 1;  %index of beta in the age sequence
header{1} = strcat('Age (',year_scale,')');
for i = 1:M
    %header{1,1+ind_alpha(i)} = strcat('Alpha_',num2str(i)); %alphas
    header{1,1+ind_alpha(i)} = strcat('Period_',num2str(i),' early boundary');
    %header{1,1+ind_beta(i)} = strcat('Beta_',num2str(i)); %betas
    header{1,1+ind_beta(i)} = strcat('Period_',num2str(i),' late boundary');
    for j = 1:N(i)
        %header{1,1+ind_b_theta(i)+j-1}=strcat('Theta_',[num2str(i),',',num2str(j)]); %thetas
        header{1,1+ind_b_theta(i)+j-1} = P(i).lab_code{j};
    end
end
%% write header and data to a file
fmt = repmat('%25g ',1,2*M+sum(N)); 
fmt = ['%25d ' fmt];
fid = fopen(filename,'wt');
fprintf(fid,'%25s ',header{1:end});
fprintf(fid,'\n');
for k = 1:size(pdfs,1)
    fprintf(fid,fmt,pdfs(k,:));
    fprintf(fid,'\n');
end    
fclose(fid);
end