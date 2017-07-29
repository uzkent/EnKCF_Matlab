
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

function [precision, fps] = run_tracker(video, kernel_type, feature_type, show_visualization, show_plots)

	%path to the videos (you'll be able to choose one with the GUI).
	base_path ='.\data';

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
        
    model = '/Users/buzkent/Documents/caffe/models/bvlc_reference_caffenet/deploy_conv2.prototxt';
    weights = '/Users/buzkent/Documents/caffe/models/bvlc_reference_caffenet/bvlc_reference_caffenet.caffemodel';
    cnn_model = caffe.Net(model, weights, 'test'); % create net and load weights
	
    padding.sroi = 2.5;  %Padding for the Small Area Translation Filter
	lambda = 1e-4;  %regularization
	output_sigma_factor.sroi = 0.1;  %spatial bandwidth (proportional to target)
	interp_factor.sroi = 0.020;  %linear interpolation factor for adaptation
	kernel.sroi_sigma = 0.6;  %gaussian kernel bandwidth

    padding.lroi = 3.0; %Padding for the Large Area Translation Filter
    output_sigma_factor.lroi = 0.1; %spatial bandwith for the Large ROI TF
    interp_factor.lroi = 0.020;  %linear interpolation factor for adaptation    
	kernel.lroi_sigma = 0.6;  %gaussian kernel bandwidth

    padding.scale = 1.0; %Padding for the Scale filter
    output_sigma_factor.scale = 0.1; %spatial bandwith for the Large ROI TF
    interp_factor.scale = 0.020;  %linear interpolation factor for adaptation    
	kernel.scale_sigma = 0.6;  %gaussian kernel bandwidth
    
	assert(any(strcmp(kernel_type, {'linear', 'polynomial', 'gaussian'})), 'Unknown kernel.')
		
    %we were given the name of a single video to process.

    %get image file names, initial state, and ground truth for evaluation
    base_path = '/Users/buzkent/Downloads/UAV123/data_seq/UAV123/';
    video = 'bike1';
    [img_files, pos, target_sz, ground_truth, video_path] = load_video_info(base_path, video);

    %call tracker function with all the relevant parameters
    [positions,~, time] = tracker(video_path, img_files, pos, target_sz, ...
        padding, kernel, lambda, output_sigma_factor, interp_factor, ...
        cell_size, features, show_visualization,cnn_model);

    %calculate and show precision plot, as well as frames-per-second
    precisions = precision_plot(positions, ground_truth, video, show_plots);
    fps = numel(img_files) / time;

    fprintf('%12s - Precision (20px):% 1.3f, FPS:% 4.2f\n', video, precisions(20), fps)

    if nargout > 0,
        %return precisions at a 20 pixels threshold
        precision = precisions(20);
    end

end
