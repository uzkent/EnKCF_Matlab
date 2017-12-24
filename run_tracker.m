%
%  High-Speed Tracking with Kernelized Correlation Filters
%
%  Joao F. Henriques, 2014
%  http://www.isr.uc.pt/~henriques/
%
%  Main interface for Kernelized/Dual Correlation Filters (KCF/DCF).
%  This function takes care of setting up parameters, loading video
%  information and computing precisions. For the actual tracking code,
%  check out the TRACKER function.
%
%  RUN_TRACKER
%    Without any parameters, will ask you to choose a video, track using
%    the Gaussian KCF on HOG, and show the results in an interactive
%    figure. Press 'Esc' to stop the tracker early. You can navigate the
%    video using the scrollbar at the bottom.
%
%  RUN_TRACKER VIDEO
%    Allows you to select a VIDEO by its name. 'all' will run all videos
%    and show average statistics. 'choose' will select one interactively.
%
%  RUN_TRACKER VIDEO KERNEL
%    Choose a KERNEL. 'gaussian'/'polynomial' to run KCF, 'linear' for DCF.
%
%  RUN_TRACKER VIDEO KERNEL FEATURE
%    Choose a FEATURE type, either 'hog' or 'gray' (raw pixels).
%
%  RUN_TRACKER(VIDEO, KERNEL, FEATURE, SHOW_VISUALIZATION, SHOW_PLOTS)
%    Decide whether to show the scrollable figure, and the precision plot.
%
%  Useful combinations:
%  >> run_tracker choose gaussian hog  %Kernelized Correlation Filter (KCF)
%  >> run_tracker choose linear hog    %Dual Correlation Filter (DCF)
%  >> run_tracker choose gaussian gray %Single-channel KCF (ECCV'12 paper)
%  >> run_tracker choose linear gray   %MOSSE filter (single channel)
%
%
%   revised by: Yang Li, August, 2014
%   http://ihpdep.github.io

function [precision, success] = run_tracker(video, kernel_type, ~, show_visualization, ~)

	%default settings
	if nargin < 1, video = 'choose'; end
	if nargin < 2, kernel_type = 'gaussian'; end
	if nargin < 3, feature_type = 'hogcolor'; end
	if nargin < 4, show_visualization = ~strcmp(video, 'all'); end
	if nargin < 5, show_plots = ~strcmp(video, 'all'); end

	%parameters according to the paper. at this point we can override
	%parameters based on the chosen kernel or feature type
	kernel.type = kernel_type;
	
    %Features for Different Correlation Filters
    features.lroi_hogcolor = true;
    features.lroi_hog_orientations = 9;
	features.sroi_hog = true;
	features.sroi_hog_orientations = 9;
    features.scale_hogcolor = true;
    features.scale_hog_orientations = 9;
    features.deep = true;
	cell_size = 4;
        
    %model_trans = '/Volumes/Burak_HardDrive/Moving_Platform_CNN_Training/VGG16/deploy_trans.prototxt';
    %model_scale = '/Volumes/Burak_HardDrive/Moving_Platform_CNN_Training/VGG16/deploy_scale.prototxt';
    %weights = '/Volumes/Burak_HardDrive/Moving_Platform_CNN_Training/VGG16/VGG_ILSVRC_16_layers.caffemodel';
    %cnn_model.trans = caffe.Net(model_trans, weights, 'test'); % create net and load weights
    %cnn_model.scale = caffe.Net(model_scale, weights, 'test'); % create net and load weights
	
    padding.sroi = 1.5;  %Padding for the Small Area Translation Filter
	lambda = 1e-4;       %regularization
	output_sigma_factor.sroi = 0.12;  %spatial bandwidth (proportional to target)
	interp_factor.sroi = 0.020;  %linear interpolation factor for adaptation
	kernel.sroi_sigma = 0.6;  %gaussian kernel bandwidth

    padding.lroi = 1.50; %Padding for the Large Area Translation Filter
    output_sigma_factor.lroi = 0.1; %spatial bandwith for the Large ROI TF
    interp_factor.lroi = 0.02;  %linear interpolation factor for adaptation    
	kernel.lroi_sigma = 0.5;  %gaussian kernel bandwidth

    padding.scale = 0.00; %Padding for the Scale filter
    output_sigma_factor.scale = 0.08; %spatial bandwith for the Large ROI TF
    interp_factor.scale = 0.020;  %linear interpolation factor for adaptation    
	kernel.scale_sigma = 0.7;  %gaussian kernel bandwidth
    
	assert(any(strcmp(kernel_type, {'linear', 'polynomial', 'gaussian'})), 'Unknown kernel.')
		
    %we were given the name of a single video to process.

    %get image file names, initial state, and ground truth for evaluation
    base_path = '/Users/buzkent/Downloads/UAV123/data_seq/UAV123/';
    fid   = fopen('/Users/buzkent/Downloads/UAV123/startFrames_UAV123.txt');
    seqInfo = textscan(fid, '%d%d%s%s\n');
    for i = 1:size(seqInfo{1},1)
        video = seqInfo{4}{i};
        seq = seqInfo{3}{i};
        firstFrame = seqInfo{1}(i);
        lastFrame = seqInfo{2}(i);
        try
            [img_files, pos, target_sz, ground_truth, video_path] = load_video_info(base_path, video, seq, firstFrame, lastFrame);

            %call tracker function with all the relevant parameters
            [positions,rect_results, fps] = tracker(video_path, img_files, pos, target_sz, ...
                padding, kernel, lambda, output_sigma_factor, interp_factor, ...
                cell_size, features, show_visualization,[]);
            i
            %calculate and show precision plot, as well as frames-per-second
            [precision(i,:),success(i,:)] = precision_plot(positions, rect_results, ground_truth, video, 0);
            mean(precision)
            mean(success)
            fps_all(i) = mean(fps);
            mean(fps_all)
        catch err
           continue;
        end
    end
    save results2.mat precision success;
end
