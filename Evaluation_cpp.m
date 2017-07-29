clear;
prFiles = dir('/Users/buzkent/Downloads/Results/Precision/');
srFiles = dir('/Users/buzkent/Downloads/Results/Success/');
prFiles(1:2) = [];
srFiles(1:2) = [];

pr = zeros(51,1);
sr = zeros(101,1);
for i = 1:size(prFiles,1)

    prFile = dlmread(['/Users/buzkent/Downloads/Results/Precision/' prFiles(i).name]);
    srFile = dlmread(['/Users/buzkent/Downloads/Results/Success/' srFiles(i).name]);
    
    for j = 1:size(prFile,1)
    
        pr(j) = pr(j) + prFile(j);
        
    end
    
    for j = 1:size(srFile,1)
    
        sr(j) = sr(j) + srFile(j);
        
    end
    
end
pr = pr ./ size(prFiles,1);
sr = sr ./ size(srFiles,1);

runTime = 1/mean(dlmread('/Users/buzkent/Downloads/Results/RunTime/RunTimes.txt'));
