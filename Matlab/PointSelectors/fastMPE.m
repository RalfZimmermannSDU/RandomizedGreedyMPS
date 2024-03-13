function [P_opt, MPE_error] = fastMPE(U0, nr_points, StartPoints, mod_switch)
%--------------------------------------------------------------------------
% Accelerated MPE point selection.
%This function implements the method introduced in
% Zimmermann, R., Willcox, K.: An accelerated greedy missing point
% estimation procedure. SIAM Journal on Scientific Computing 38(5),
% 2827–2850 (2016)
% see Alg. B.1 "Efficient greedy step based on (3.7)"
%
% Input: 
%     U0     = (POD) basis matrix
% nr_points  = number of magic points to be selected
%StartPoints = initial set of magic points
% mod_switch = various choices for targeting different 
%              eigenvalues in the greedy procedure
%              mod_switch == 1
%                    address penultimate eigenvalue every third iteration
%
%              mod_switch > 1
%                    address eigenvalue according to growth potential (4.2)
%--------------------------------------------------------------------------


[n,p] = size(U0);

% list of values for the error bound
MPE_error = [];

%initialize list of MPE candidate points
P_init = 1:n;

%initialize list of MPE points
P_opt = [StartPoints];

P_init = setdiff(P_init, P_opt);

%--------------------
% start greedy alg:
%--------------------
tau = 0.05;  % gap threshold between two consecutive eig.vals
counter = 1;

% keep track of which modified eigvals are addressed
target_highest_dk = 0;     %check if targeting alg. falls back to dk = dmin
counter_targetdk = zeros(size(U0,2), 1); % count how often each dk is addressed

while size(P_opt, 2) < nr_points
    [Psi, sigmas, Phi] = svd(U0(P_opt,:),0);
            
    %************************************************
    %screening
    %************************************************
    ScreenMat = U0(P_init,:)*Phi;
    cols_ScreenMat = size(ScreenMat,2);
    v_onlinescreen = zeros(size(ScreenMat,1),1);
    d = diag(sigmas).*diag(sigmas); % compute eig.vals from sing.vals
    dim = cols_ScreenMat;

    if mod_switch == 1
        % switch cyclically between the last three components
        % this is Alg. 4.1 from the paper
        if mod(counter-1, 3) < 2
            shift_this = dim;
        else
            shift_this = dim-1;
        end
    elseif mod_switch > 1
        % switch between the components based on the growth potential
        % this is Alg. 4.2 from the paper
        shift_this = dim;
        for c = 1:dim-1
            if (d(dim-c) - d(dim-c+1))/d(dim-c) > tau
                shift_this = dim-c+1;
                counter_targetdk(c) = counter_targetdk(c) +1;
                break;
            end
            if c == dim-1
                target_highest_dk = target_highest_dk +1;
            end
        end
    else
        shift_this = dim;
    end
    counter = counter + 1;
    
    %----------------------------------------------------------------------
    % solve the eigval approximations based on the quadratic equation
    %----------------------------------------------------------------------
    d_last = d(shift_this);
    d_pen  = d(shift_this -1);
    c_up_vec = d - d_pen; % vector minus constant in each component
    c_up_vec(1:shift_this-2) = 1.0./c_up_vec(1:shift_this-2);
    c_up_vec(shift_this:end) = 1.0./c_up_vec(shift_this:end);
    %omit components 'dim-shift_this -1' and 'dim-shift_this'
    c_up_vec(shift_this -1) = 0.0;
    c_up_vec(shift_this) = 0.0;
    gap1 = d_last + d_pen;
    gap2 = d_last*d_pen;

    %precompute squared absolute values of the rows of ScreenMat
    % conj() function necessary for complex case
    z_sq_mat = ScreenMat.*conj(ScreenMat);
    c_ups    = ones(size(z_sq_mat,1),1) + z_sq_mat*c_up_vec;
    alphas1  = gap1 + (z_sq_mat(:,shift_this-1)        + z_sq_mat(:,shift_this))./c_ups;
    alphas2  = gap2 + (z_sq_mat(:,shift_this-1)*d_last + z_sq_mat(:,shift_this)*d_pen)./c_ups;

    v_onlinescreen = d_pen - 0.5*alphas1  + sign(c_ups).*sqrt(0.25*alphas1.*alphas1 - alphas2);

    % sort in ascending order. This means that the indices of the rows
    % with the best screening value come first
    [v_screen_sort, ind_sort] = min(v_onlinescreen);

    %add optimal index to current MPE selection
    P_opt = [P_opt, P_init(ind_sort)];
    


    %remove best index from list so that it is not chosen more than once
    P_init = setdiff(P_init, P_init(ind_sort(1)));
end
% compute value of error bound associated with the current point
% selection
sigmas_opt = svd(U0(P_opt,:), 0);
sig_mpe = 1.0/sigmas_opt(length(sigmas_opt));
MPE_error = [MPE_error, sig_mpe];


% uncomment for some onscreen output of interest
%counter_targetdk
%target_highest_dk
%sum(counter_targetdk) + target_highest_dk
%sqrt(d)

return;
end

