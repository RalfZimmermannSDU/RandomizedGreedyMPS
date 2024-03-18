%
% This script implements a symplectic one-step method
% for the 
% "nonlinear 1D Schroedinger" equation as featured in
% Section 3.3 of 
%
% "Randomized greedy magic point selection schemes for 
%  nonlinear model reduction"
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
beta          = 1.0; % weight parameter for leverage score sampling
r             = 40;
rDEIM         = 40;
nr_points     = 2*rDEIM% choose point selection scheme
nr_randruns   = 100; % set to 100 for reproducing the paper results
loadFOMstring = 'Snapshots_Schroed/SchroedNew_FOM_N400_tsteps2000_Tend20.mat';
%
stream = RandStream('mt19937ar');
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

if model == 1         
    % FOM
    N         = input("Enter FOM dimension (for reproducing the paper results(N=400)): ")
    state_dim = N;  % extra variable that works for FOM and ROM
    % set up FOM operator and initial states
    [Dxx, q0, p0, x_space] = prepare_FOM(N, L, c);
else 
    % i.e. in POD or POD-DEIM mode: 
    % It is expected that snapshots are precomputed and stored under
    % "loadFOMstring" !!!
    state_dim = r;
    % load full order snapshots
    FOMstruct = load(loadFOMstring);
    FOM       = FOMstruct.FOM
    N         = FOM.dimN
    n_tsteps  = FOM.n_tsteps
    %% create ROM struct
    [ROM] = prepare_PODROM(FOM, r);
    if model == 3
        % create DEIM bases and DEIM magic points 
        % for a_nonlin, b_nonlin
        ROM = prepare_PODDEIMROM(FOM, ROM, rDEIM);
    end
end


% right-hand side lin_ops defined as functions grad_qH, grad_pH
% function handles with fixed lin_op and parameter epsilon
% these are independent from any random point selection
if model == 1  % FOM
    % FOM is not random
    % do only one full-order run, set back to 1, if set otherwise
    nr_randruns = 1;
    disp("in FOM mode")
    grad_qH_fun = @(q,p) grad_qH(q,p, Dxx, epsilon);
    grad_pH_fun = @(q,p) grad_pH(q,p, Dxx, epsilon);
elseif model == 2 % POD-ROM
    % POD ROM is not random
    % do only one POD-ROM run, set back to 1, if set otherwise
    nr_randruns = 1;
    disp("in POD mode")
    grad_qH_fun = @(q,p) grad_qH_POD(q,p, ROM.lin_op,  epsilon, ROM.Uq, ROM.Up);
    grad_pH_fun = @(q,p) grad_pH_POD(q,p, ROM.lin_op', epsilon, ROM.Uq, ROM.Up);
end

%% storage of  trajectory  and initial states
q_traj = zeros(state_dim,n_tsteps+1);   % state_dim is N for FOM, r for ROM
p_traj = zeros(state_dim,n_tsteps+1);
% for recording results
MPE_err_b    = zeros(1,nr_randruns);
MPE_err_a    = zeros(1,nr_randruns);
traj_error_q = zeros(1,nr_randruns);
traj_error_p = zeros(1,nr_randruns);

% set initial states
if model == 1 % FOM mode
    q_traj(:,1) = q0;
    p_traj(:,1) = p0;
else % i.e. POD or POD-DEIM
    %project initial state
    q_traj(:,1) = ROM.Uq'*ROM.Yq(:,1);
    p_traj(:,1) = ROM.Up'*ROM.Yp(:,1);
end

for k=1:nr_randruns
    if model>2 % i.e., POD-DEIM
        %% enhance the magic points
        % the DEIM approximation depends on the random selection
        % --> has to be recomputed in every run
        [Pa, a_DEIM_op, MPE_error_a] = magicPointSelection(stream, ...
                                                         ROM.Up, ...
                                                         ROM.Qa, ...
                                                         ROM.Pa0,...
                                                         nr_points,...
                                                         2,...
                                                         beta,...
                                                         MPE_mode);
        MPE_err_a(k) = MPE_error_a;
        [Pb, b_DEIM_op, MPE_error_b] = magicPointSelection(stream, ...
                                                         ROM.Uq, ...
                                                         ROM.Qb, ...
                                                         ROM.Pb0,...
                                                         nr_points,...
                                                         2,...
                                                         beta,...
                                                         MPE_mode);
        MPE_err_b(k) = MPE_error_b;
        % right hand side operatores in reduced POD-DEIM ODE
        disp("in DEIM mode")
        grad_qH_fun = @(q,p) grad_qH_POD_DEIM(q,p, ROM.lin_op,...
                                              epsilon, ROM.Uq, ROM.Up,...
                                              a_DEIM_op, Pa);
        grad_pH_fun = @(q,p) grad_pH_POD_DEIM(q,p, ROM.lin_op',...
                                              epsilon, ROM.Uq, ROM.Up,...
                                              b_DEIM_op, Pb);    
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
        traj_error_q(k) = norm(ROM.Uq*q_traj - ROM.Yq, 'fro')/norm(ROM.Yq, 'fro');
        traj_error_p(k) = norm(ROM.Up*p_traj - ROM.Yp, 'fro')/norm(ROM.Yp, 'fro');
    end
end % k=1:nr_rand_runs


if model == 1 % in FOM mode, write data to disc
    %  all snapshots, q stacked above p
    FOM.snapshots = [q_traj;p_traj];
    FOM.lin_op    = Dxx;
    FOM.dimN      = N;
    FOM.epsilon   = epsilon;
    FOM.n_tsteps  = n_tsteps;
    FOM.x_space   = x_space;
    % collect DEIM snapshots
    FOM.a_nonlin  = (q_traj.^2 + p_traj.^2).*q_traj;
    FOM.b_nonlin  = (q_traj.^2 + p_traj.^2).*p_traj;
    eval(['save Snapshots_Schroed/', 'SchroedNew_FOM_N', num2str(N), '_tsteps', num2str(n_tsteps), '_Tend', num2str(T), '.mat FOM']);
end


%% from down here: display and illustration of results

% in POD or DEIM mode, the trajectories have to be lifted
% to the FOM dimension before plotting
if model > 1 % i.e. POD or POD-DEIM
    p_traj = ROM.Up*p_traj;
    q_traj = ROM.Uq*q_traj;
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
    disp(['       mean_MPE_a = ', num2str(mean(MPE_err_a))]);
    disp(['        var_MPE_a = ', num2str(var(MPE_err_a))]);
    disp(['       mean_MPE_b = ', num2str(mean(MPE_err_b))]);
    disp(['        var_MPE_b = ', num2str(var(MPE_err_b))]);

    disp(['mean_error_traj_q = ', num2str(mean(traj_error_q))]);
    disp([' var_error_traj_q = ', num2str(var(traj_error_q))]);  
    disp(['mean_error_traj_p = ', num2str(mean(traj_error_p))]);
    disp([' var_error_traj_p = ', num2str(var(traj_error_p))]);


    % for illustration purposes: Compute quantum probabilities
    % and compare ROM to FOM

    Schroedprob_traj    = sqrt(p_traj.*p_traj + q_traj.*q_traj);
    SchroedprobFOM_traj = sqrt(ROM.Yp.*ROM.Yp + ROM.Yq.*ROM.Yq);

    lw = 1.0; %linewidth
    figure
    plot(FOM.x_space, SchroedprobFOM_traj(:, floor(n_tsteps/2) + 1), 'k-',  ...
         FOM.x_space, SchroedprobFOM_traj(:, n_tsteps + 1), 'k-',...
         FOM.x_space, Schroedprob_traj(:, floor(n_tsteps/2) + 1), 'r-.',  ...
         FOM.x_space, Schroedprob_traj(:, n_tsteps + 1), 'r-.', 'LineWidth', lw)
    xlabel("x")
    ylabel("|u|")
    legend('FOM t = 10','FOM t = 20','ROM t = 10', 'ROM t = 20', 'Location','northwest')


    % plot sequence of snapshots
    figure;
    subplot(1,2,1)
    for k = 1:50:n_tsteps+1;
        plot3( (k/n_tsteps)*ones(N,1), FOM.x_space, SchroedprobFOM_traj(:,k), 'k-'); 
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
        plot3( (k/n_tsteps)*ones(N,1), FOM.x_space, Schroedprob_traj(:,k), 'r-'); 
        hold on
    end
    xlim( [0 1])
    ylim([-30 30])
    zlim([0 1.6])
    box on
    grid on
end





