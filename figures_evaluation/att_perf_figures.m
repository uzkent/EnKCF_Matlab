% This function computes the precision and success curves w.r.t to the
% target attributes
% Home Directory for the Object Specific Results
clear;
close all;
home_dir_pr = '/Users/buzkent/Downloads/Results_EnKCF_cpp/Precision/';
home_dir_sr = '/Users/buzkent/Downloads/Results_EnKCF_cpp/Success/';
att_dir = '/Users/buzkent/Downloads/UAV123/anno/UAV123/att/';
att_number = 12;
list_atts{1} = 'SV'; list_atts{2} = 'ARC'; list_atts{3} = 'LR'; list_atts{4} = 'FM';
list_atts{5} = 'FOC'; list_atts{6} = 'POC'; list_atts{7} = 'OV'; list_atts{8} = 'BC';
list_atts{9} = 'IV'; list_atts{10} = 'VC'; list_atts{11} = 'CM'; list_atts{12} = 'SOB';

%% Preicison and Success Figures
seqs = dir(home_dir_pr);
for i = 1:att_number
    pr{i} = zeros(51,1);
    sr{i} = zeros(101,1);
    counter = 0;
    % Iterate over sequences
    for j = 3:size(seqs,1)

        % Read the Attribute File
        att_object = dlmread([att_dir seqs(j).name]);
        
        % If the attribute of interest is 1, count it
        if att_object(i) == 0
           pr{i} = pr{i} + dlmread([home_dir_pr seqs(j).name]);
           sr{i} = sr{i} + dlmread([home_dir_sr seqs(j).name]);
           counter = counter + 1;
        end
    end
    pr{i} = pr{i} ./ counter;
    sr{i} = sr{i} ./ counter;
    counter
end

%% Have the Curve Colors for the Trackers
plotDrawStyleAll={   struct('color',[1,0,0],'lineStyle','-','lineWidth',1.5),...
    struct('color',[0,1,0],'lineStyle','-','lineWidth',1.5),...
    struct('color',[0,0,1],'lineStyle','-','lineWidth',1.5),...
    struct('color',[0,0,0],'lineStyle','-','lineWidth',1.5),...%    struct('color',[1,1,0],'lineStyle','-'),...%yellow
    struct('color',[1,0,1],'lineStyle','-','lineWidth',1.5),...%pink
    struct('color',[0,1,1],'lineStyle','-','lineWidth',1.5),...
    struct('color',[0.5,0.5,0.5],'lineStyle','-','lineWidth',1.5),...%gray-25%
    struct('color',[136,0,21]/255,'lineStyle','-','lineWidth',1.5),...%dark red
    struct('color',[255,127,39]/255,'lineStyle','-','lineWidth',1.5),...%orange
    struct('color',[0,162,232]/255,'lineStyle','-','lineWidth',1.5),...%Turquoise
    struct('color',[163,73,164]/255,'lineStyle','-','lineWidth',1.5),...%purple    %%%%%%%%%%%%%%%%%%%%
    struct('color',[0.7,0.7,0.7],'lineStyle','-','lineWidth',1.5),...%gray-25%
    struct('color',[136,145,21]/255,'lineStyle','-','lineWidth',1.5),...%dark red
    struct('color',[20,12,75]/255,'lineStyle','-','lineWidth',1.5),...%orange
    struct('color',[0,82,23]/255,'lineStyle','-','lineWidth',1.5),...%Turquoise
    struct('color',[123,93,130]/255,'lineStyle','-','lineWidth',1.5),...%purple
    };
% Rank the attributes w.r.t AUC and PR at 20 px.
for i = 1:size(list_atts,2)
    AUC(i) = round(mean(sr{i}),3);
    Prec(i) = round(pr{i}(21),3);    
end

% Sort them from larger to smaller values
[Val, AUC_Ind] = sort(AUC,'descend');
[Val, Prec_Ind] = sort(Prec,'descend');
% Draw the curve for each attribute - PR
for i = 1:att_number
    figure(1);
    plot(0:50, pr{Prec_Ind(i)}, plotDrawStyleAll{Prec_Ind(i)});
    hold on;
    legendInfoPrec{i} = [list_atts{Prec_Ind(i)} '[' num2str(pr{Prec_Ind(i)}(21)) ']'];    
end
% Add Legends and Make the figure high resolution.
h2 = legend(legendInfoPrec);
xlabel('Location Error Threshold','fontsize',10);
ylabel('Precision','fontsize',10);
set(gcf,'Units','inches','Position',[0 0 3.5 3.5]);
set(h2,'FontSize',8);
set(gca,'XTick',(0:5:50));

% Draw the curve for each attribute - SR
for i = 1:att_number
    figure(2);
    plot(0:.01:1, sr{AUC_Ind(i)}, plotDrawStyleAll{AUC_Ind(i)});
    hold on;
    legendInfoPrec{i} = [list_atts{AUC_Ind(i)} '[' num2str(mean(mean(sr{AUC_Ind(i)}))) ']'];    
end
% Add Legends and Make the figure high resolution.
h1 = legend(legendInfoPrec);
set(gcf,'Units','inches','Position',[0 0 4.0 3.5]);
set(h1,'FontSize',8);
xlabel('Overlap Threshold','fontsize',10);
ylabel('Success Rate','fontsize',10);
set(gca,'XTick',(0:0.1:1));