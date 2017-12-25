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
temp = load('w2crs');
w2c = temp.w2crs;
	%if the target is large, lower the resolution, we don't need that much
	%detail
	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
	if resize_image,
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
    end

    % mod factor for running the ensemble of filters
    mod_divisor = 5;
    % Scale Confidence for the Candidates other than 1
    scale_confidence = 0.95;
	%window size, taking padding into account
    if target_sz(1) > target_sz(2)
        window_sz.sroi(1) = 96;
        window_sz.lroi(1) = 96;
        window_sz.scale(1) = 64;
        window_sz.sroi(2) = round(96/target_sz(1) * target_sz(2));
        window_sz.lroi(2) = round(96/target_sz(1) * target_sz(2));
        window_sz.scale(2) = round(64/target_sz(1) * target_sz(2));
    else
        window_sz.sroi(2) = 96;
        window_sz.lroi(2) = 96;
        window_sz.scale(2) = 64;
        window_sz.sroi(1) = round(96/target_sz(2) * target_sz(1));
        window_sz.lroi(1) = round(96/target_sz(2) * target_sz(1));
        window_sz.scale(1) = round(64/target_sz(2) * target_sz(1));
    end
    
    % window size, taking padding into account
	window_sz.sroi = floor(target_sz * (1 + padding.sroi));
	window_sz.lroi = floor(target_sz * (1 + padding.lroi));
	window_sz.scale = floor(target_sz * (1 + padding.scale));
    
    %Small ROI Translation Filter
	%create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma.sroi = sqrt(prod(target_sz)) * output_sigma_factor.sroi / cell_size;
	yf.sroi = fft2(gaussian_shaped_labels(output_sigma.sroi, ...
        floor(window_sz.sroi / cell_size)));

	%store pre-computed cosine window
	cos_window.sroi = hann(size(yf.sroi,1)) * hann(size(yf.sroi,2))';
    
    %Large ROI Translation Filter
    %create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma.lroi = sqrt(prod(target_sz)) * output_sigma_factor.lroi ...
        / cell_size;
	yf.lroi = fft2(gaussian_shaped_labels(output_sigma.lroi, ...
        floor(window_sz.lroi / cell_size)));

	%store pre-computed cosine window
	cos_window.lroi = hann(size(yf.lroi,1)) * hann(size(yf.lroi,2))';
    
    %Scale Filter
    %create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma.scale = sqrt(prod(target_sz)) * output_sigma_factor.scale ...
        / cell_size;
	yf.scale = fft2(gaussian_shaped_labels(output_sigma.scale, ...
        floor(window_sz.scale / cell_size)));

	%store pre-computed cosine window
	cos_window.scale = hann(size(yf.scale,1)) * hann(size(yf.scale,2))';
	
	search_size = [1  1.05 1.00/1.05]; % Scale Space

	%note: variables ending with 'f' are in the Fourier domain.

	time = 0;  %to calculate FPS
	positions = zeros(numel(img_files), 2);  %to calculate precision
	rect_results = zeros(numel(img_files), 4);  %to calculate 
	for frame = 1:numel(img_files)
		%load image
		im = imread([video_path img_files{frame}]);
    	if resize_image
            im = imresize(im,0.5); % Downsample image
        end
        %Resize the Image - Can be Disabled
        tic();
        szid = 1;
		if frame > 1
            %% TRANSLATION ESTIMATION
            if (mod(frame, mod_divisor) > 0 && mod(frame, mod_divisor) < 3) % Large ROI Trans.
                %obtain a subwindow for detection at the position from last
                %frame, and convert to Fourier domain (its size is unchanged)
                %patch = get_subwindow(im, pos, window_sz);
                tmp_sz = floor(target_sz * (1 + padding.lroi));
                patch = get_subwindow(im, pos, tmp_sz);
                patch = imresize(patch, window_sz.lroi);
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
                tmp_sz = floor(target_sz * (1 + padding.lroi));
                scale_ratio = tmp_sz(2) / window_sz.lroi(2);
                pos = pos + scale_ratio*cell_size * [vert_delta - 1, horiz_delta - 1];
            end
            if (mod(frame, mod_divisor) == 3 || mod(frame, mod_divisor) == 4) % Small ROI Trans.
                %obtain a subwindow for detection at the position from last
                %frame, and convert to Fourier domain (its size is unchanged)
                tmp_sz = floor(target_sz * (1 + padding.sroi));
                patch = get_subwindow(im,pos,tmp_sz);
                patch = imresize(patch,window_sz.sroi);
                zf = fft2(get_features(patch, features, cell_size, cos_window,w2c,2,[]));

                %calculate response of the classifier at all shifts
                kzf = gaussian_correlation(zf, model_xf.sroi, kernel.sroi_sigma);

                response_sroi = real(ifft2(model_alphaf.sroi .* kzf));  %equation for fast detection
                %target location is at the maximum response. we must take into
                %account the fact that, if the target doesn't move, the peak
                %will appear at the top-left corner, not at the center (this is
                %discussed in the paper). the responses wrap around cyclically.
                [vert_delta, horiz_delta] = find(response_sroi == max(response_sroi(:)), 1);

                if vert_delta > size(zf,1) / 2  %wrap around to negative half-space of vertical axis
                    vert_delta = vert_delta - size(zf,1);
                end
                if horiz_delta > size(zf,2) / 2  %same for horizontal axis
                    horiz_delta = horiz_delta - size(zf,2);
                end
                % Translation Update                
                tmp_sz = floor(target_sz * (1 + padding.sroi));
                scale_ratio = tmp_sz(2)/window_sz.sroi(2);
                pos = pos + scale_ratio * cell_size * [vert_delta - 1, horiz_delta - 1];           
            end
            
            %% SCALE ESTIMATION
            %obtain a subwindow for detection at the position from last
			%frame, and convert to Fourier domain (its size is unchanged)
			%patch = get_subwindow(im, pos, window_sz);
            response = [];
            if (mod(frame, mod_divisor) == 0) % Scale Filter
                for i=1:size(search_size,2)
                    tmp_sz = floor(target_sz * (1 + padding.scale) * search_size(i));
                    patch = get_subwindow(im, pos, tmp_sz);
                    patch = imresize(patch, window_sz.scale);
                    zf = fft2(get_features(patch, features, cell_size,... 
                    cos_window,w2c,1,[]));

                    %calculate response of the classifier at all shifts
                    kzf = gaussian_correlation(zf, model_xf.scale, kernel.scale_sigma);

                    response(:,:,i) = real(ifft2(model_alphaf.scale .* kzf));  %equation for fast detection
                end
                % New Scale Index to Update Scale
                response(:,:,2:3) = scale_confidence * response(:,:,2:3);
                [~,tmp,~] = find(response == max(response(:)), 1);
                szid = floor((tmp-1)/(size(cos_window.scale,2)))+1;
            end
        end

        %% TRAINING        
        % Update Scale First
        target_sz = target_sz * search_size(szid);
        % Update Large ROI Trans.
        if ((mod(frame, mod_divisor) > 0 && mod(frame, mod_divisor) < 3) || (frame == 1))
            tmp_sz = floor(target_sz * (1 + padding.lroi));
            patch = get_subwindow(im,pos, tmp_sz);
            patch = imresize(patch, window_sz.lroi);            
            xf.lroi = fft2(get_features(patch, features, cell_size,...
                cos_window,w2c,0,[]));

            %Kernel Ridge Regression, calculate alphas (in Fourier domain)
            kf = gaussian_correlation(xf.lroi, xf.lroi, kernel.lroi_sigma);

            alphaf.lroi = yf.lroi ./ (kf + lambda);   %equation for fast training
            flag = 0;
        end
        
        %Update SROI Trans.
        if (mod(frame, mod_divisor) == 3 || mod(frame, mod_divisor) == 4 || frame == 1)
            tmp_sz = floor(target_sz * (1 + padding.sroi));
            patch = get_subwindow(im, pos, tmp_sz);
            patch = imresize(patch, window_sz.sroi);            
            xf.sroi = fft2(get_features(patch, features, cell_size,...
                cos_window, w2c, 2, []));

            %Kernel Ridge Regression, calculate alphas (in Fourier domain)
            kf = gaussian_correlation(xf.sroi, xf.sroi, kernel.sroi_sigma);

            alphaf.sroi = yf.sroi ./ (kf + lambda);   %equation for fast training
            flag = 1;
        end
        %Update Scale Filter
        if (mod(frame, mod_divisor) == 0 || frame == 1)
            tmp_sz = floor(target_sz * (1 + padding.scale));
            patch = get_subwindow(im, pos, tmp_sz);
            patch = imresize(patch, window_sz.scale);
            xf.scale = fft2(get_features(patch, features, cell_size,...
                cos_window, w2c, 1, []));

            %Kernel Ridge Regression, calculate alphas (in Fourier domain)
            kf = gaussian_correlation(xf.scale, xf.scale, kernel.scale_sigma);

            alphaf.scale = yf.scale ./ (kf + lambda);   %equation for fast training 
            flag = 2;
        end
        
		if frame == 1  %first frame, train with a single image
			model_alphaf.lroi = alphaf.lroi;
			model_xf.lroi = xf.lroi;
            model_alphaf.sroi = alphaf.sroi;
			model_xf.sroi = xf.sroi;
            model_alphaf.scale = alphaf.scale;
			model_xf.scale = xf.scale;
		else
			%subsequent frames, interpolate model - LROI
            if flag == 0
                model_alphaf.lroi = (1 - interp_factor.lroi) * model_alphaf.lroi ...
                    + interp_factor.lroi * alphaf.lroi;
                model_xf.lroi = (1 - interp_factor.lroi) * model_xf.lroi + ...
                    interp_factor.lroi * xf.lroi;
            end
            %subsequent frames, interpolate model - SROI
			if flag == 1
                model_alphaf.sroi = (1 - interp_factor.sroi) * model_alphaf.sroi + ...
                    interp_factor.sroi * alphaf.sroi;
                model_xf.sroi = (1 - interp_factor.sroi) * model_xf.sroi + ...
                    interp_factor.sroi * xf.sroi;
            end
            %subsequent frames, interpolate model - Scale
			if flag == 2
                model_alphaf.scale = (1 - interp_factor.scale) * model_alphaf.scale + ...
                    interp_factor.scale * alphaf.scale;
                model_xf.scale = (1 - interp_factor.scale) * model_xf.scale + ...
                    interp_factor.scale * xf.scale;
            end
        end

		%save position and timing
		positions(frame,:) = pos;
        fps_all(frame) = toc();
        fps = mean(fps_all);
		box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        rect_results(frame,:) = box;
        
        %Display the results
        figure(1); imshow(im);
        rectangle('Position',box);
        drawnow;
    end

	if resize_image % Scale the new bounding box
		positions = positions * 2;
        rect_results = rect_results * 2;
	end
end

