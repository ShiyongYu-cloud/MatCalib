function [] = save_pooled_mean(filename,pooled_pdfs)
%% function for saving the empirical pdfs to a file
%INPUT
%filename: file name to which for data to be saved
%pooled_pdfs: probabilities of pooled mean age
%% creat a header
header = cell(1,3); 
header{1} = 'Mean_age (years)';
header{2} = 'Probability';
header{3} = 'HPD_68.2%';
header{4} = 'HPD_95.4%';
%% write the header and data to a file
fmt = repmat('%25g ',1,3); 
fmt = ['%25d ' fmt];
fid = fopen(filename,'wt');
fprintf(fid,'%25s ',header{1:end});
fprintf(fid,'\n');
for k = 1:size(pooled_pdfs,1)
    fprintf(fid,fmt,pooled_pdfs(k,:));
    fprintf(fid,'\n');
end    
fclose(fid);
end