function [fname] = getfnames(datapath,cam)
% grab all file names from folder that start with a number or letter

fnames = dir([datapath,'/',cam,'*']);
count = 0;

for i = 1:length(fnames)
    if ~startsWith(fnames(i).name,'.')
        if ~startsWith(fnames(i).name,'System ')
            count = count+1;
            fname(count,:) = fnames(i).name;
        end
    end
end