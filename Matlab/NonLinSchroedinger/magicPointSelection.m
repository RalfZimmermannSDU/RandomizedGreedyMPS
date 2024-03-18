function [Pout, DEIM_op, MPE_error] = magicPointSelection(stream, ...
                                                          Upod, ...
                                                          Qdeim, ...
                                                          P0,...
                                                          nr_points,...
                                                          criterion,...
                                                          beta,...
                                                          MPE_mode)
%--------------------------------------------------------------------------
% wrapper function for magic point selection,
% perform either
% * random  greedy oversampling          (MPE_mode==1)
% * random leverage score oversampling   (MPE_mode==2)
% * or stay with original DEIM selection (MPE_mode, else case)
%--------------------------------------------------------------------------

if MPE_mode==1
    disp('random greedy oversampling')
    [Pout, ~] = randfastMPE(stream, Qdeim, nr_points, P0, criterion);
elseif MPE_mode==2 
    disp('leverage score oversampling')
    % leverage score point selection
    % with beta = 0, this is uniform random point selection
    [Pout, ~] = leveragescoreMPE(stream, Qdeim, nr_points, P0, beta);
else % use initial DEIM only
    disp('classic DEIM, no oversampling')
    Pout = P0;
end % of magic point selection
        
% oblique DEIM operators
% for b-term
PTQ      = Qdeim(Pout,:);
[Phi, sigmas, Psi]  = svd(PTQ, 0);
MPE_error = 1.0/sigmas(end,end);
DEIM_op = (Upod'*Qdeim)*Psi*((diag(1.0./diag(sigmas)))*Phi');
end 