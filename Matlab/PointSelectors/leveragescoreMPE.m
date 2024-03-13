function [P_opt, MPE_error] = leveragescoreMPE(stream, U0, nr_points, StartPoints, beta)
% select magic points randomly according to their leverage scores
%
% This follows Section 5 of
% Saibaba, A.K.: Randomized discrete empirical interpolation method for
% nonlinear model reduction. SIAM Journal on Scientific Computing 42(3),
% 1582–1608 (2020)
%
% MATLAB' statistics and machine learning toolbox is required
%
%
% Inputs
% stream       = random stream for reproducibility
% U0           = basis matrix of dims Nxr
% nr_points    = point total to be selected from the N rows
% StartPoiunts = list of r points already selected
% beta         = parameter: convex combi of uniform and LS distribution
%                beta = 1:LS sampling, beta = 0, uniform sampling
% leverage score lj is norm of row j U0(j,:);

[N,r] = size(U0);

LS = sum(U0.*U0,2);     % sum up over columns
probs_LS = zeros(1,N);  % initialize as array, so that the next line
                        % also works with beta = 0.
probs_LS = (beta/r)*LS + (1-beta)/N;

% the statistics and machine learning toolbox provides
% a method to sample according to prescribed probabilities

k     = length(StartPoints);  % points already drawn
P_opt = StartPoints;
while length(P_opt) < nr_points
    % draw nr_points-k samples from LS-distribution
    R = randsample(stream, 1:N, nr_points-k, true, probs_LS);
    P_opt = unique([R, P_opt]);
    k = length(P_opt)
end

% compute oblique projection error
sigmas  = svd(U0(P_opt,:), 0);
MPE_error = 1.0/sigmas(length(sigmas));

end