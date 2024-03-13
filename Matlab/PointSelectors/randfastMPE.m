function [P_opt, MPE_error] = randfastMPE(stream, ...
                                            U0,...
                                            nr_points,...
                                            StartPoints,...
                                            mode)
% This is Algorithm 3 from 
%
% Randomized greedy magic point selection schemes for nonlinear model reduction
% Ralf Zimmermann and Kai Cheng, ACOM
%
% Input: 
% stream     = fixed (matlab  built-in) stream of random numbers
%     U0     = (POD) basis matrix
% nr_points  = number of magic points to be selected
%StartPoints = initial set of magic points
% mode       = switch to change between the various screening criteria
%              sc0, sc1, sc2 of the associated paper.

[n,p] = size(U0);
P_init = 1:n;       %initialize list of MPE candidate points
P_opt = zeros(1,nr_points);   %initialize list of MPE points
P_opt(1:p) = StartPoints;
n_rand_preselect = 2*nr_points;    %nr points to be screened
[~, Sigma_up, Phi_up] = svd(U0(P_opt(1:p),:), 0); 
for counter = p:nr_points-1
    % candidate set of random points:
    P_rand = randi(stream, size(U0,1), 1, n_rand_preselect);
    P_init = unique(setdiff(P_rand, P_opt(1:counter)));
    if mode == 1                             % this is (sc1) from the paper
        ScreenMat = U0(P_init,:)*Phi_up(:, p-1:p);        
        v_onlinescreen = abs(ScreenMat(:,1)).^2 - ...
                         abs(ScreenMat(:,2)).^2;
    elseif mode == 2                         % this is (sc2) from the paper
        ScreenMat = U0(P_init,:)*Phi_up(:, p-1:p);
        v_onlinescreen = abs(ScreenMat(:,1)) +...
                         1.0./abs(ScreenMat(:,2));
    else                                     % this is (sc0) from the paper
        ScreenMat = U0(P_init,:)*Phi_up(:,p); 
        v_onlinescreen = -abs(ScreenMat);
    end
    [~, ind_sort] = min(v_onlinescreen);                     % pick minimum
    P_opt(counter+1) = P_init(ind_sort);               % keep optimal index
    u_up = U0(P_init(ind_sort),:)*Phi_up;                      % update Phi
    [~, Sigma_up, Phi_tilde] = svd([Sigma_up;u_up], 0); 
    Phi_up = Phi_up*Phi_tilde;
    MPE_error = 1.0/Sigma_up(end,end);                    % for the records
end % end "for counter = p:nr_points-1"
return;
end
