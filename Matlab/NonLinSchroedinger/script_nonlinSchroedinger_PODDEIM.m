%
% This script implements a symplectic one-step method
% for the 
% "nonlinear 1D Schroedinger" equation as featured in
% Section 3.3 of 
%
% Randomized greedy magic point selection schemes for nonlinear model reduction
% Ralf Zimmermann and Kai Cheng, ACOM
%
%
%
%
%clear; close all;
%
% add path to time stepping schemes
addpath('Symp_time_stepping_schemes')
addpath('../PointSelectors/')
%
%% USER PARAMETERS
model         = 1;   % 1:FOM,            2:POD,               3:DEIM
tstep_scheme  = 2;   % 1:Stormer-Verlet, 2:symplectic Euler
MPE_mode      = 1;   % 1:onlinefast MPE, 2:leverage score sampling, 
                     % else: no oversampling
beta          = 0.0  % weight parameter for leverage score sampling
r             = 40;
rDEIM         = 40;
nr_points     = 2*rDEIM% choose point selection scheme
nr_randruns   = 1; % set to 100 for reproducing the paper results
loadFOMstring = 'Snapshots_Schroed/SchroedNew_FOM_N400_tsteps2000_Tend20.mat';
%
stream = RandStream('mt19937ar')
%
%
%% model paramters, cf. Maboudi Afkham, Hesthaven
l = 0.11;
L = 2*pi/l;
c = 1.0;
epsilon = 1.0932;
dt= 0.01;
T = 20.0;
n_tsteps = floor(T/dt);

if model == 1         % do a full-order run
    N         = input("Enter FOM dimension (for reproducing the paper results(N=400)): ")
    state_dim = N;  % extra variable that works for FOM and ROM
    dx= L/N;        % spatial resolution
    % set up FD operator
    e = ones(N,1);
    Dxx = spdiags([e -2*e e], -1:1, N, N);
    Dxx(N,1) = 1;
    Dxx(1,N) = 1;
    Dxx = Dxx/dx^2;
    % spatial discretization
    x0 = 0.0;
    x_space = linspace(x0-L/2,x0+L/2,N+1)';
    % remove last point for use of periodic boundary conditions
    x_space(end) = [];
    % initial condition
    u0 = (sqrt(2.0)./cosh(x_space-x0)).*exp(i*(c/2)*(x_space-x0));
    q0 = imag(u0);
    p0 = real(u0);
else 
    % i.e. in POD or POD-DEIM mode: 
    % It is expected that snapshots are precomputed!
    state_dim = r;
    % load full order snapshots
    FOMstruct = load(loadFOMstring);
    FOM       = FOMstruct.FOM
    N         = FOM.dimN
    n_tsteps  = FOM.n_tsteps
    % FOM.snapshots contains the q- and p-snapshots stacked as
    %     |q1,q2,...,qm|
    % Y = |p1,p2,...,pm| 
    %
    Yq        = FOM.snapshots(1:N,:);
    Yp        = FOM.snapshots(N+1:2*N,:);
    Dxx       = FOM.operator;
    x_space   = FOM.x_space;
    a_nonlin  = FOM.a_nonlin(:,1:n_tsteps);
    b_nonlin  = FOM.b_nonlin(:,1:n_tsteps);
    % create separate POD bases for Yq and Yp
    % set up POD states: Uq*qr(t) ~ q(t),  Up*pr(t) ~ p(t)
    [Uq,Sq,Vq] = svd(Yq);
    Uq         = Uq(:,1:r);
    [Up,Sp,Vp] = svd(Yp);
    Up         = Up(:,1:r);
    dSq = diag(Sq);
    dSp = diag(Sp);
    figure;
    semilogy(1:N, diag(Sq), 1:N, diag(Sq), 'k-', 1:10:N, dSq(1:10:N), 'k*',...
             1:N, diag(Sp), 1:N, diag(Sp), 'b--', 5:10:N, dSq(5:10:N), 'bs');
    legend('Singular values of Yq','Singular values of Yp')

    % the projected ODE is
    % qr(t) =  (Uq'*Dxx*Up)*pr + eps*Uq'*b(Uq*qr(t),Up*pr(t)) = \nabla_pH
    % pr(t) = -(Up'*Dxx*Uq)*qr + eps*Up'*a(Uq*qr(t),Up*pr(t)) =-\nabla_qH

    % The linear operator is projected as
    Dxxr = Up'*(Dxx*Uq); 
    % Dxx is sym. => Dxxr' = (Up'*Dxx*Uq)
    % Dxxr  features in grad_qH
    % Dxxr' features in grad_pH
    if model > 2
        % create DEIM bases for a_nonlin, b_nonlin
        [Qa,Ea,Ra] = svd(a_nonlin);
        Qa         = Qa(:,1:rDEIM);
        [Qb,Eb,Rb] = svd(b_nonlin);
        Qb         = Qb(:,1:rDEIM);
        %figure;
        %semilogy(1:N, diag(Eb), 1:N, diag(Ea));
        % perform DEIM, output is point lists Pa, Pb
        Pa0 = deim(Qa);
        Pb0 = deim(Qb);
    end
end


% right-hand side operators defined as functions grad_qH, grad_pH
% function handles with fixed operator and parameter epsilon
if model == 1  % FOM
    disp("in FOM mode")
    grad_qH_fun = @(q,p) grad_qH(q,p, Dxx, epsilon);
    grad_pH_fun = @(q,p) grad_pH(q,p, Dxx, epsilon);
elseif model == 2 % POD-ROM
    disp("in POD mode")
    grad_qH_fun = @(q,p) grad_qH_POD(q,p, Dxxr,  epsilon, Uq, Up);
    grad_pH_fun = @(q,p) grad_pH_POD(q,p, Dxxr', epsilon, Uq, Up);
end

%% storage of  trajectory  and initial states
q_traj = zeros(state_dim,n_tsteps+1);
p_traj = zeros(state_dim,n_tsteps+1);

% set initial states
if model == 1 % FOM mode
    q_traj(:,1) = q0;
    p_traj(:,1) = p0;
else % i.e. POD or POD-DEIM
    %project initial state
    q_traj(:,1) = Uq'*Yq(:,1);
    p_traj(:,1) = Up'*Yp(:,1);
end

% for recording results
MPE_err_b    = zeros(1,nr_randruns);
MPE_err_a    = zeros(1,nr_randruns);
traj_error_q = zeros(1,nr_randruns);
traj_error_p = zeros(1,nr_randruns);
for k=1:nr_randruns
    if model>2 % i.e., POD-DEIM
        %% enhance the magic points
        if MPE_mode==1
            disp('random greedy oversampling')
            [Pa, ~] = randfastMPE(stream, Qa, nr_points, Pa0, 2);
            [Pb, ~] = randfastMPE(stream, Qb, nr_points, Pb0, 2);
        elseif MPE_mode==2 
            disp('leverage score oversampling')
            % leverage score point selection
            % with beta = 0, this is uniform random point selection
            [Pa, ~] = leveragescoreMPE(stream, Qa, nr_points, Pa0, beta);
            [Pb, ~] = leveragescoreMPE(stream, Qb, nr_points, Pb0, beta);
        else % use initial DEIM only
            disp('classic DEIM, no oversampling')
            Pa = Pa0;
            Pb = Pb0;
        end % of magic point selection
        
        % oblique DEIM operators
        % for b-term
        PTQb      = Qb(Pb,:);
        [Phi_b, sigmas_b, Psi_b]  = svd(Qb(Pb,:), 0);
        MPE_err_b(k) = 1.0/sigmas_b(end,end);
        b_DEIM_op = (Uq'*Qb)*Psi_b*((diag(1.0./diag(sigmas_b)))*Phi_b');
        % for a-term
        PTQa      = Qa(Pa,:);
        [Phi_a, sigmas_a, Psi_a]  = svd(Qa(Pa,:), 0);
        MPE_err_a(k) = 1.0/sigmas_a(end,end);
        a_DEIM_op = (Up'*Qa)*Psi_a*((diag(1.0./diag(sigmas_a)))*Phi_a');
        % right hand side operatores in reduced POD-DEIM ODE
        disp("in DEIM mode")
        grad_qH_fun = @(q,p) grad_qH_POD_DEIM(q,p, Dxxr,  epsilon, Uq, Up, a_DEIM_op, Pa);
        grad_pH_fun = @(q,p) grad_pH_POD_DEIM(q,p, Dxxr', epsilon, Uq, Up, b_DEIM_op, Pb);    
    end% end DEIM MODE


%% the actual time step loop
    tic;
    for j=1:n_tsteps
        %disp(['t = ', num2str(j*dt)]);
        if tstep_scheme == 1
            % time stepping scheme
            [qrn_p1, prn_p1] = StormerVerlet_step(dt, q_traj(:,j), p_traj(:,j),...
                                      grad_qH_fun, grad_pH_fun);
        else
            [qrn_p1, prn_p1] = SympEuler_step(dt, q_traj(:,j), p_traj(:,j),...
                                        grad_qH_fun, grad_pH_fun);
        end
        q_traj(:,j+1)  = qrn_p1;
        p_traj(:,j+1)  = prn_p1;
    end
    looptime = toc;

    disp(['Time step loop took ', num2str(looptime), 's.'])
    if model > 1  %i.e POD or POD DEIM
        % compute reconstruction errors
        traj_error_q(k) = norm(Uq*q_traj - Yq, 'fro')/norm(Yq, 'fro');
        traj_error_p(k) = norm(Up*p_traj - Yp, 'fro')/norm(Yp, 'fro');
    end
end % k=1:nr_rand_runs


if model == 1 % in FOM mode, write data to disc
    %  all snapshots, q stacked above p
    FOM.snapshots = [q_traj;p_traj];
    FOM.operator  = Dxx;
    FOM.dimN      = N;
    FOM.epsilon   = epsilon;
    FOM.n_tsteps  = n_tsteps;
    FOM.x_space   = x_space;
    % collect DEIM snapshots
    FOM.a_nonlin  = (q_traj.^2 + p_traj.^2).*q_traj;
    FOM.b_nonlin  = (q_traj.^2 + p_traj.^2).*p_traj;
    eval(['save Snapshots_Schroed/', 'SchroedNew_FOM_N', num2str(N), '_tsteps', num2str(n_tsteps), '_Tend', num2str(T), '.mat FOM']);
end


% in POD or DEIM mode, the trajectories have to be lifted
% to the FOM dimension before plotting
if model > 1 % i.e. POD or POD-DEIM
    p_traj = Up*p_traj;
    q_traj = Uq*q_traj;
end
% plot some graphs
pt0 = p_traj(:,1); 
qt0 = q_traj(:,1);

pt10 = p_traj(:, floor(n_tsteps/2) + 1); 
qt10 = q_traj(:, floor(n_tsteps/2) + 1);

pt20 = p_traj(:,n_tsteps+1); 
qt20 = q_traj(:,n_tsteps+1);


if model>1 % i.e. POD or POD-DEIM
    % error statisics
    mean_MPE_a = mean(MPE_err_a)
    var_MPE_a = var(MPE_err_a)
    mean_MPE_b = mean(MPE_err_b)
    var_MPE_b = var(MPE_err_b)

    mean_error_traj_q = mean(traj_error_q)
    mean_error_traj_p = mean(traj_error_p)
    var_error_traj_q = var(traj_error_q)
    var_error_traj_p = var(traj_error_p)


    % for illustration purposes: Compute quantum probabilities
    % and compare ROM to FOM

    Schroedprob_traj    = sqrt(p_traj.*p_traj + q_traj.*q_traj);
    SchroedprobFOM_traj = sqrt(Yp.*Yp + Yq.*Yq);

    lw = 1.0; %linewidth
    figure
    plot(x_space, SchroedprobFOM_traj(:, floor(n_tsteps/2) + 1), 'k-',  ...
     x_space, SchroedprobFOM_traj(:, n_tsteps + 1), 'k-',...
     x_space, Schroedprob_traj(:, floor(n_tsteps/2) + 1), 'r-.',  ...
     x_space, Schroedprob_traj(:, n_tsteps + 1), 'r-.', ...
     'LineWidth', lw)
    xlabel("x")
    ylabel("|u|")
    legend('FOM t = 10','FOM t = 20','ROM t = 10', 'ROM t = 20', 'Location','northwest')


    % plot sequence of snapshots
    figure;
    subplot(1,2,1)
    for k = 1:50:n_tsteps+1;
        plot3( (k/n_tsteps)*ones(N,1), x_space, SchroedprobFOM_traj(:,k), 'k-'); 
        hold on
    end
    xlim( [0 1])
    ylim([-30 30])
    zlim([0 1.6])
    box on
    grid on
    hold on
    subplot(1,2,2)
    for k = 1:50:n_tsteps+1;
        plot3( (k/n_tsteps)*ones(N,1), x_space, Schroedprob_traj(:,k), 'r-'); 
        hold on
    end
    xlim( [0 1])
    ylim([-30 30])
    zlim([0 1.6])
    box on
    grid on
end





