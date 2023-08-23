function [] = save_age_diff(filename,diff_pdfs)
%% function for saving the empirical pdfs to a file
%INPUT
%filename: file name to which for data to be saved
%diff_pdfs: probabilities of age difference
%% creat a header
header = cell(1,3); 
header{1} = 'Age_diff (years)';
header{2} = 'Probability';
header{3} = 'HPD_68.2%';
header{4} = 'HPD_95.4%';
%% write the header and data to a file
fmt = repmat('%25g ',1,3); 
fmt = ['%25d ' fmt];
fid = fopen(filename,'wt');
fprintf(fid,'%25s ',header{1:end});
fprintf(fid,'\n');
for k = 1:size(diff_pdfs,1)
    fprintf(fid,fmt,diff_pdfs(k,:));
    fprintf(fid,'\n');
end    
fclose(fid);
end