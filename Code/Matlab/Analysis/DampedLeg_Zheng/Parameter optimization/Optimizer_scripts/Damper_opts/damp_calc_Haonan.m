%% Parameter optimization for big rat leg

% Optimizing repsonse b value of each robot joint only based on gathered
% response data.

clear; clc; close('all');

% Add paths needed for loading data and using functions
addpath('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Zheng\Parameter optimization\Optimizer functions and data')
addpath('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Zheng\Parameter optimization\Results')

% Create folder to save results to and add file path
% folder_title = 'Robot_damping_calcs_Haonan';
% mkdir(strcat('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Krnacik\Haonan\Parameter optimization\Results\', folder_title))
% addpath(strcat('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Krnacik\Haonan\Parameter optimization\Results\', folder_title)) 

%% Load and format data to be used

% Choose mechanical system to optimize
sysName = 'MechPropBigRat';
sysProp = load('-mat', sysName);

% Load the data file for all joint data
DampData = load('-mat', 'DampDataHaonanK119');

% Set start time step
time_step = 0.001;                                                  % [ s ]

%% Set optimization specifications (IC, BC, options)

% Set initial conditions based on expected damper values
IC_data = [0.7 0.105 0.063];

% Set options for the search function
options = optimset('PlotFcns','optimplotfval', 'MaxIter', 500);

% Set search boundary conditions based on a scalar multiple of inital guess
scalar_b = 5;

UB = scalar_b .* IC_data;
LB = IC_data ./ scalar_b;

% Set high damping constants to be used to hold joints steady
set_var = [1000 1000 1000];

%% Load EOM

% Define equations of motion from saved EOM file
load('EOM_no_springs.mat')
fprintf('\nEOM loaded.\n')

% Compile EOM into single variable
dui = [du1, du2, du3, du4, du5, du6];

%% OPTIMIZATION 

% Compile non-optimized variables into single structure
data.dui = dui;
data.DampData = DampData;
data.set_var = set_var;
data.time_step = time_step;
data.sysProp = sysProp;

% Run fminsearch on cost function for all three joints
b_sols = [ 0 0 0 ];

for n = 1:3
    
fprintf('Running optimization for hip...\n')
[b_sol,fval1,exitflag,output] = fminsearchbnd(@(U) cost_calc_damp(U, data, n),IC_data(n),LB(n),UB(n),options);
fprintf('Optimization for hip finished.\n\n')

b_sols(n) = b_sol;

end

%% DISPLAY RESULTS 

% Run ODE with achieved parameters and save results, only for one trial per
% joint 
fprintf('\nParameter optimization achieved. Producing results.\n')

% Produce figure
fig = figure( 'Color', 'w');
    
% Create array of b values to use
b = [b_sols(1) set_var(2) set_var(3); set_var(1) b_sols(2) set_var(3); set_var(1) set_var(2) b_sols(3)];

% Define time
time = DampData.time;

% Define universal constants.
g_value = 9.81;                             % [m/s^2]

rt_ode = [0 0 0];
rt_data = [0 0 0];
joints = ["Hip" "Knee" "Ankle"];

% State symbolic variables used in EOM solver
syms M1 M2 M3;
syms theta1(t) dtheta1(t) ddtheta1(t) theta2(t) dtheta2(t) ddtheta2(t) theta3(t) dtheta3(t) ddtheta3(t);
syms L1 L2 L3;
syms R1 R2 R3;
syms b1 b2 b3 b4 b5 b6;
syms I1 I2 I3;
syms g;
syms a1 a2 a3 w1 w2 w3;
syms u1 u2 u3 u4 u5 u6 du1 du2 du3 du4 du5 du6;
syms tau1 tau2 tau3


m1_value = data.sysProp.m1_value; m2_value = data.sysProp.m2_value; m3_value = data.sysProp.m3_value;
R1_value = data.sysProp.R1_value; R2_value = data.sysProp.R2_value; R3_value = data.sysProp.R3_value;
L1_value = data.sysProp.L1_value; L2_value = data.sysProp.L2_value; L3_value = data.sysProp.L3_value;

% Define x limits for each trial
lim = [0 2; 0 2; 0 5];

for n = 1:3
    
    b1_value = b(n,1);
    b2_value = b(n,2);
    b3_value = b(n,3);
    
    if n == 1
        thetas = DampData.hip_dat; 
        time = DampData.time.hip;
    elseif n == 2
        thetas = DampData.knee_dat; 
        time = DampData.time.knee;
    else
        thetas = DampData.ankle_dat;
        time = DampData.time.ankle;
    end

    % Assign equations of motion
    du1 = dui(1); 
    du2 = dui(2); 
    du3 = dui(3); 
    du4 = dui(4); 
    du5 = dui(5); 
    du6 = dui(6);

    % Substitute numerical values into the dynamical system flow.
    du1_value = subs( du1, [ L1 L2 L3 R1 R2 R3 M1 M2 M3 b1 b2 b3 g ], [ L1_value L2_value L3_value R1_value R2_value R3_value m1_value m2_value m3_value b1_value b2_value b3_value g_value ] );
    du2_value = subs( du2, [ L1 L2 L3 R1 R2 R3 M1 M2 M3 b1 b2 b3 g ], [ L1_value L2_value L3_value R1_value R2_value R3_value m1_value m2_value m3_value b1_value b2_value b3_value g_value ] );
    du3_value = subs( du3, [ L1 L2 L3 R1 R2 R3 M1 M2 M3 b1 b2 b3 g ], [ L1_value L2_value L3_value R1_value R2_value R3_value m1_value m2_value m3_value b1_value b2_value b3_value g_value ] );
    du4_value = subs( du4, [ L1 L2 L3 R1 R2 R3 M1 M2 M3 b1 b2 b3 g ], [ L1_value L2_value L3_value R1_value R2_value R3_value m1_value m2_value m3_value b1_value b2_value b3_value g_value ] );
    du5_value = subs( du5, [ L1 L2 L3 R1 R2 R3 M1 M2 M3 b1 b2 b3 g ], [ L1_value L2_value L3_value R1_value R2_value R3_value m1_value m2_value m3_value b1_value b2_value b3_value g_value ] );
    du6_value = subs( du6, [ L1 L2 L3 R1 R2 R3 M1 M2 M3 b1 b2 b3 g ], [ L1_value L2_value L3_value R1_value R2_value R3_value m1_value m2_value m3_value b1_value b2_value b3_value g_value ] );

    % Create anonymous functions from the dynamical system flow components.
    fdu1_temp = matlabFunction( du1_value ); fdu1 = @( t, u, tau ) fdu1_temp( u(2) );
    fdu2_temp = matlabFunction( du2_value ); fdu2 = @( t, u, tau ) fdu2_temp( tau(1), tau(2), tau(3), u(1), u(2), u(3), u(4), u(5), u(6) );
    fdu3_temp = matlabFunction( du3_value ); fdu3 = @( t, u, tau ) fdu3_temp( u(4) );
    fdu4_temp = matlabFunction( du4_value ); fdu4 = @( t, u, tau ) fdu4_temp( tau(1), tau(2), tau(3), u(1), u(2), u(3), u(4), u(5), u(6) );
    fdu5_temp = matlabFunction( du5_value ); fdu5 = @( t, u, tau ) fdu5_temp( u(6) );
    fdu6_temp = matlabFunction( du6_value ); fdu6 = @( t, u, tau ) fdu6_temp( tau(1), tau(2), tau(3), u(1), u(2), u(3), u(4), u(5), u(6) );
    fdu = @( t, u, tau ) [ fdu1( t, u, tau ); fdu2( t, u, tau ); fdu3( t, u, tau ); fdu4( t, u, tau ); fdu5( t, u, tau ); fdu6( t, u, tau ) ];

    % Define the applied torques.
    taus = zeros( 3, 1 );

    % Define the initial joint angles.
    theta0s = [thetas(1,1), thetas(1,2), thetas(1,3)];
    
    % Calculate initial joint angular velocity
    w0s = zeros(3, 1);
   
    % Assemble the state variable initial condition.
    u0s = [ theta0s( 1 ); w0s( 1 ); theta0s( 2 ); w0s( 2 ); theta0s( 3 ); w0s( 3 ) ];

    % Simulate the triple pendulum dynamics.
    [ t, x ] = ode15s( @( t, u ) fdu( t, u, taus ), time, u0s );

    % Calculate cost by risetime
    rt_ode(n) = risetime_ek(x( :,(n*2-1) ),t);
    rt_data(n) = risetime_ek(thetas( :,n ),time);

    subplot(2,3,n)
    plot(time, rad2deg(thetas( :,n )), ':b', 'LineWidth', 1); hold on
    plot(t, rad2deg(x( :,(n*2-1) )), '-r', 'LineWidth', 1)
    title(sprintf('Damped %s response', joints(n)))
    xlabel('Time(s)'); ylabel('Joint angle (deg)')
    xlim(lim(n,:))
    ylim([0 185])
    
    % Create scaled plot
    x_scale = x( :,(n*2-1) ); x_scale = x_scale - x_scale(1); x_scale = x_scale/x_scale(end);
    thetas_scale = thetas(:,n); thetas_scale = thetas_scale - thetas_scale(1); thetas_scale = thetas_scale/thetas_scale(end);
    
    subplot(2, 3, n+3)
    plot(time, thetas_scale, ':b', 'LineWidth', 1); hold on
    plot(t, x_scale, '-r', 'LineWidth', 1)
    title(sprintf('Damped scaled %s response', joints(n)))
    xlabel('Time(s)'); ylabel('Scaled response')
    xlim(lim(n,:))
    ylim([0 1.2])
end
legend('Actual response', 'ODE response')
%% Save results

% save results - optimized parameters and figure
% save(strcat('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Krnacik\Haonan\Parameter optimization\Results\', folder_title, '\Results_HaonanK118_adjusted.mat'), 'b_sols', 'rt_ode', 'rt_data')
% saveas(fig, strcat('C:\GitHub\Quadruped_Robot\Code\Matlab\Analysis\DampedLeg_Krnacik\Haonan\Parameter optimization\Results\', folder_title, '\Results_HaonanK118_adjusted.fig'))
% fprintf('Data saved.\n')



