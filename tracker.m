function [positions, rect_results, fps] = tracker(video_path, img_files, pos, target_sz, ...
	padding, kernel, lambda, output_sigma_factor, interp_factor, cell_size, ...
	features, ~, ~)
%TRACKER Kernelized/Dual Correlation Filter (KCF/DCF) tracking.
%   This function implements the pipeline for tracking with the KCF (by
%   choosing a non-linear kernel) and DCF (by choosing a linear kernel).
%
%   It is meant to be called by the interface function RUN_TRACKER, which
%   sets up the parameters and loads the video information.
%
%   Parameters:
%     VIDEO_PATH is the location of the image files (must end with a slash
%      '/' or '\').
%     IMG_FILES is a cell array of image file names.
%     POS and TARGET_SZ are the initial position and size of the target
%      (both in format [rows, columns]).
%     PADDING is the additional tracked region, for context, relative to 
%      the target size.
%     KERNEL is a struct describing the kernel. The field TYPE must be one
%      of 'gaussian', 'polynomial' or 'linear'. The optional fields SIGMA,
%      POLY_A and POLY_B are the parameters for the Gaussian and Polynomial
%      kernels.
%     OUTPUT_SIGMA_FACTOR is the spatial bandwidth of the regression
%      target, relative to the target size.
%     INTERP_FACTOR is the adaptation rate of the tracker.
%     CELL_SIZE is the number of pixels per cell (must be 1 if using raw
%      pixels).
%     FEATURES is a struct describing the used features (see GET_FEATURES).
%     SHOW_VISUALIZATION will show an interactive video if set to true.
%
%   Outputs:
%    POSITIONS is an Nx2 matrix of target positions over time (in the
%     format [rows, columns]).
%    TIME is the tracker execution time, without video loading/rendering.
%
%   Joao F. Henriques, 2014
%
%   revised by: Yang Li, August, 2014
%   http://ihpdep.github.io


addpath('./utility');
temp = load('/Users/buzkent/GITHUB/HKCF_Tracker/w2crs');
w2c = temp.w2crs;
	%if the target is large, lower the resolution, we don't need that much
	%detail
	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
	if resize_image,
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
    end
    
    % window size, taking padding into account
	window_sz.lroi = floor(target_sz * (1 + padding.lroi));
    
    %Large ROI Translation Filter
    %create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma.lroi = sqrt(prod(target_sz)) * output_sigma_factor.lroi ...
        / cell_size;
	yf.lroi = fft2(gaussian_shaped_labels(output_sigma.lroi, ...
        floor(window_sz.lroi / cell_size)));

	%store pre-computed cosine window
	cos_window.lroi = hann_window(size(yf.lroi,1))' * hann_window(size(yf.lroi,2));
   
	time = 0;  %to calculate FPS
	positions = zeros(numel(img_files), 2);  %to calculate precision
	rect_results = zeros(numel(img_files), 4);  %to calculate 
	for frame = 1:numel(img_files)
		%load image
		im = imread([video_path img_files{frame}]);
    	if resize_image
            im = im_resize(im,0.5); % Downsample image
        end
        %Resize the Image - Can be Disabled
        tic();
		if frame > 1
            %% TRANSLATION ESTIMATION
            %obtain a subwindow for detection at the position from last
            %frame, and convert to Fourier domain (its size is unchanged)
            %patch = get_subwindow(im, pos, window_sz);
            tmp_sz = floor(target_sz * (1 + padding.lroi));
            patch = get_subwindow(im, pos, tmp_sz);
            zf = fft2(get_features(patch, features, cell_size, ...
                cos_window,w2c, 0, []));

            %calculate response of the classifier at all shifts
            kzf = gaussian_correlation(zf, model_xf.lroi, kernel.lroi_sigma);

            response_lroi(:,:) = real(ifft2(model_alphaf.lroi .* kzf));  %equation for fast detection

            %target location is at the maximum response. we must take into
            %account the fact that, if the target doesn't move, the peak
            %will appear at the top-left corner, not at the center (this is
            %discussed in the paper). the responses wrap around cyclically.
            [vert_delta, horiz_delta] = find(response_lroi == ...
                max(response_lroi(:)), 1);

            % From feature space to image space translation
            if vert_delta > size(zf,1) / 2  %wrap around to negative half-space of vertical axis
                vert_delta = vert_delta - size(zf,1);
            end
            if horiz_delta > size(zf,2) / 2  %same for horizontal axis
                horiz_delta = horiz_delta - size(zf,2);
            end
            % Translation Update                
			pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];
           
        end

        %% TRAINING        
        % Update Large ROI Trans.
        tmp_sz = floor(target_sz * (1 + padding.lroi));
        patch = get_subwindow(im,pos, tmp_sz);
        xf.lroi = fft2(get_features(patch, features, cell_size,...
            cos_window,w2c,0,[]));

        %Kernel Ridge Regression, calculate alphas (in Fourier domain)
        kf = gaussian_correlation(xf.lroi, xf.lroi, kernel.lroi_sigma);

        alphaf.lroi = yf.lroi ./ (kf + lambda);   %equation for fast training
   
		if frame == 1  %first frame, train with a single image
			model_alphaf.lroi = alphaf.lroi;
			model_xf.lroi = xf.lroi;     
        else
            model_alphaf.lroi = (1 - interp_factor.lroi) * model_alphaf.lroi ...
                + interp_factor.lroi * alphaf.lroi;
            model_xf.lroi = (1 - interp_factor.lroi) * model_xf.lroi + ...
                interp_factor.lroi * xf.lroi;
        end

		%save position and timing
		positions(frame,:) = pos;
        time = toc();
        fps(frame) = 1/time;
		box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        rect_results(frame,:) = box;
        
%         %Display the results
        figure(1); imshow(im);
        rectangle('Position',box);
        drawnow;
    end

	if resize_image % Scale the new bounding box
		positions = positions * 2;
        rect_results = rect_results * 2;
	end
end

