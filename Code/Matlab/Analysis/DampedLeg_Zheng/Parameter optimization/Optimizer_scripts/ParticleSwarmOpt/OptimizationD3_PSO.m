%% Parameter Swarm Optimization for D3 log leg - underdamped

% This script is intended to optimize the dog leg to a scaled version of
% the rat leg, changing frequency and damping as outline in supplementary
% material from Sutton et al to achieve a less damped dog leg

clear; clc; close('all');

% Add paths needed for loading data and using functions
addpath('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Zheng\Parameter optimization\Optimizer functions and data')
addpath('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Zheng\Parameter optimization\IC_check')
addpath('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Zheng\Parameter optimization\dyn_prop_calc')

% Create folder to save results to and add file path
% folder_title = 'PSOD3_results';
% mkdir(strcat('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Krnacik\Haonan\Parameter optimization\Results\', folder_title))
% addpath(strcat('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Krnacik\Haonan\Parameter optimization\Results\', folder_title)) 
% clc;


%% Load and format data to be used

% Choose mechnical system to optimize
sysName = 'MechPropDog';
sysProp = load('-mat', sysName);

% Load calculated dynamic properties of each joint
jdataProp = load('-mat', 'jointdyn3');

% Load the data file for all joint data
load('-mat', 'jdata');

% Load start and end indices
load('-mat', 'start_indices')
load('-mat', 'end_indices')

% Define muscles and trials, respectively, to optimize to (chosen by me)
muscles = [1 2 3 4 5 6 7];
trials  = [5 1 1 1 1 1 1];

% Define which trials for each joint to place higher value on. Rows
% correspond to hip, knee, and ankle, and columns refer to each trial, 1
% indicating a higher value should be placed
trial_value = 5 *  [1 1 1 0 0 1 0; ...
                    1 0 1 1 0 1 0; ...
                    1 0 1 0 0 1 0]';
                
trial_value(trial_value==0) = 1;

% Define muscle names
muscle_names = {'IP', 'GS', 'ST', 'ST2', 'VL', 'BFp', 'BFa'};
data_title = cell(1,length(muscles));

% Create data titles
for n_i = 1:length(muscles)
   data_title(n_i) = strcat('Muscle', muscle_names(n_i), 'OPT');
end

% Set start time step
time_step = 0.001;                                                  % [ s ]


%% Set optimization specifications (IC, BC, options)

% Load results from hand optimized starting conditions
load('-mat', 'HandOPTD3_IC');
IC_data     = [b_i(1) k_i(1) b_i(2) k_i(2) b_i(3) k_i(3)];

% Set options for the particle swarm optimization
n_particles     = 20;
n_iterations    = 70;

% Set search boundary conditions
UB = [ 3 35 0.3 4 0.2 1 ];
LB = [ 0.01 1 0.001 0.1 0.001 0.01 ];


%% Load EOM

% Define equations of motion from saved EOM file
load('EOM.mat')
fprintf('\nEOM loaded.\n')

% Compile EOM into single variable
dui = [du1, du2, du3, du4, du5, du6];



%% OPTIMIZATION 

% Compile non-optimized variables into single structure
data.dui                    = dui;                  % Symbolic EOM
data.jdata                  = jdata;                % All joint angle data
data.thetabias_sym          = thetabias_sym;        % symbolic spring offset
data.muscles                = muscles;              % Muscles to be used
data.trials                 = trials;               % Trials to be used
data.start_indices          = start_indices;        % Start indices for data
data.end_indices            = end_indices;          % End indices for data
data.time_step              = time_step;            % Time step for interpolation
data.UB                     = UB;                   % Upper boundary conditions for param
data.LB                     = LB;                   % Lower boundary conditions for param
data.sysProp                = sysProp;              % Mechanical system properties
data.jdataProp              = jdataProp;            % Dynamic properties for system
data.trial_value            = trial_value;          % indicator for which trials to place higher cost



% Run fminsearch on cost function
fprintf('Beginning particle swarm optimization...\n')
[ PSO_results ] = pso_fsn_start(data, IC_data, n_particles, n_iterations);
fprintf('\nSaving results and creating graphs.\n')
 

%% DISPLAY AND SAVE RESULTS 

fig = figure(1);
% Save convergence graph
% saveas(fig, strcat('C:\Users\krnac\OneDrive\Desktop\School\Dynamic leg\Krnacik\Parameter optimization\Results\', folder_title, '\PSOconvergence.fig'))

% Create variable to save cost calculation parameters for the hip, knee,
% and ankle
JAdiff_sim          = zeros(length(muscles), 1);
risetime_sim        = zeros(length(muscles), 1);
omega_sim           = zeros(length(muscles), 1);
zeta_sim            = zeros(length(muscles), 1);

% Produce figure
fig = figure( 'Color', 'w');
    
% Run simulation on each trial with solved values
for m = 1:length(muscles)
    
    % Choose muscle and trial, and starting index value (check figure from
    % "plotjdat" in RawDataPlottingProcessing folder. All starting values
    % have been manually chosen and saved to "start_indices" data file.
    muscle      = muscles(m);
    trial       = trials(m);
    start_index = start_indices(muscle, trial);
    end_index   = end_indices(muscle, trial);

    % From chosen muscle and trial extract joint angles and time data. Note that
    % first cell of jdata corresponds to theta values, the second corresponds
    % to time values, and the third corresponds to "force" values.
    time    = jdata{2}{muscle, trial}(start_index:end_index);                      % [s]
    time    = time - time(1);
    thetas  = jdata{1}{muscle, trial}(start_index:end_index, :) * (2 * pi)/360;  % [rad]

    % Interpolate the data using spline method                                    
    thetas  = interp1(time, thetas, 0:time_step:time(end), 'spline');
    time    = 0:time_step:time(end);
    time    = time * 2;
    
    % Calculate starting joint angle velocites
    w0s = [ ( thetas(2,1) - thetas(1,1) ) ( thetas(2,2) - thetas(1,2) ) ( thetas(2,3) - thetas(1,3) ) ] / time_step;
    
    % Set parameters to run ODE
    param = [ w0s  PSO_results.gbest(1) PSO_results.gbest(3) PSO_results.gbest(5) PSO_results.gbest(2) PSO_results.gbest(4) PSO_results.gbest(6) 0 0 0];

    % Run ODE
    [ time_ode, thetas_ode, thetabias_val ] = dynamics_func(param, dui, time, thetas, thetabias_sym, sysName);

    
    
%%%%%%%%%%%% DYNAMIC PROPERTY CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate the risetime for simulation. 
    risetime_sim(m) = risetime_ek(thetas_ode(:,1), time);
    
    % Find peaks and location of sim response and actual response for knee
    [pks_x, loc_x] = findpeaks(thetas_ode(:,3)); [pks_xneg, loc_xneg] = (findpeaks(-thetas_ode(:,3)));   

    % Consolidate into single array
    pks_x = abs([pks_x' pks_xneg']); loc_x = [loc_x' loc_xneg'];

    % Define peaks based on offset from steady state value
    SteadyStateX = thetas_ode(end,3);
    pks_x = abs(pks_x - SteadyStateX);

    % sort positive and negative peaks into singular array
    [pks_x, I_x] = sort(pks_x, 'descend'); loc_x = loc_x(I_x);

    if length(pks_x) > 2

        % Calculate log dec and damping ratio
        logdec_x = log(pks_x(1) / pks_x(3));
        zeta_x = ( logdec_x / (2*pi) ) / sqrt( 1 + (logdec_x / (2*pi))^2);

        % Calculate frequency
        omega_x = 2* pi * (0.002 * (loc_x(3) - loc_x(1)))^-1;

    else
        % if system is over/criticlly damped, assign arbitrarily low
        % omega_x and high zeta
        omega_x = 0.1;
        zeta_x = 300;
    end

    % Save data into array
    omega_sim(m)    = omega_x;
    zeta_sim(m)     = zeta_x;
    
    % Calculate average joint angle difference for the ankle
	JAdiff_sim(m) = sum( thetas_ode(:,5) - thetas(:,3) ) / length(time);
    
    
    
    % Create plot of simulated vs optimized response
    subplot(2,4,m); hold on
    
    % Plot hip angles
    plot(time, rad2deg(thetas(:,1)), '--k'); plot(time_ode, rad2deg(thetas_ode(:,1)), '-k')
    
    % Plot knee angles
    plot(time, rad2deg(thetas(:,2)), '--r'); plot(time_ode, rad2deg(thetas_ode(:,3)), '-r')
    
    % Plot ankle angles
    plot(time, rad2deg(thetas(:,3)), '--b'); plot(time_ode, rad2deg(thetas_ode(:,5)), '-b')
    
    % Graph set-up
    title(sprintf('%s D3 PSO ', data_title{m}))
    xlabel('Time (s)')
    ylabel('Joint angles (deg)')
    xlim([0 1])
    ylim([60 180])
   

end

legend('Scaled rat hip', 'D3 hip', 'Scaled rat knee', 'D3 knee', 'Scaled rat ankle', 'D3 ankle')
    
% save results - optimized parameters and figure
% save(strcat('C:\Users\krnac\OneDrive\Desktop\School\Dynamic leg\Krnacik\Parameter optimization\Results\', folder_title, '\AllTrialsPSOResults.mat'), 'PSO_results', 'thetabias_val', 'risetime_sim', 'omega_sim', 'zeta_sim', 'JAdiff_sim')
% saveas(fig, strcat('C:\Users\krnac\OneDrive\Desktop\School\Dynamic leg\Krnacik\Parameter optimization\Results\', folder_title, '\AllTrialsPSOresults.fig'))
% fprintf('Data saved.\n')



