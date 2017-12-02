function [precisions, success] = precision_plot(positions, rect_results, ground_truth, title, show)
%PRECISION_PLOT
%   Calculates precision for a series of distance thresholds (percentage of
%   frames where the distance to the ground truth is within the threshold).
%   The results are shown in a new figure if SHOW is true.
%
%   Accepts positions and ground truth as Nx2 matrices (for N frames), and
%   a title string.
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/

	
	max_threshold = 50;  %used for graphs in the paper
	max_threshold_success = 100;
	precisions = zeros(max_threshold, 1); %Initiate Precisions
	success = zeros(max_threshold_success, 1); %Initiate Precisions    
	
	if size(positions,1) ~= size(ground_truth,1),
% 		fprintf('%12s - Number of ground truth frames does not match number of tracked frames.\n', title)
		
		%just ignore any extra frames, in either results or ground truth
		n = min(size(positions,1), size(ground_truth,1));
		positions(n+1:end,:) = [];
		ground_truth(n+1:end,:) = [];
	end
	
	%calculate distances to ground truth over all frames
	distances = sqrt((positions(:,1) - ground_truth(:,1)).^2 + ...
				 	 (positions(:,2) - ground_truth(:,2)).^2);
	distances(isnan(distances)) = [];

	%compute precisions
	for p = 1:max_threshold,
		precisions(p) = nnz(distances <= p) / numel(distances);
    end
    %Compute Success Overlap
    rect_results(:,1:2) = positions(:,1:2);
    for i = 1:size(ground_truth,1)
        if isnan(ground_truth(i,1))% Handle NaN Ground Truth
            continue
        end
        intersectionArea = rectint(rect_results(i,:),ground_truth(i,:));
        unionArea = (rect_results(i,3)*rect_results(i,4))+(rect_results(i,3)*rect_results(i,4))-intersectionArea;
        iou(i) = 100*intersectionArea/unionArea;
    end
	
    for s = 1:max_threshold_success,
		success(s) = nnz(iou >= s) / numel(iou);
    end
	
end

