%% Masuda Network Simulation

% This script simulates the Masuda network.

% Clear Everything.
clear, close('all'), clc


%% Define Neuron Properties.

% Define the number of neurons.
num_neurons = 8;                                                                                        % [#] Number of Neurons in the Network.

% Define Hip Extensor Motor Neuron Properties (Neuron 1).
Cm1 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm1 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er1 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R1 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am1 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm1 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm1 = R1;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah1 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh1 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh1 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna1 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh1_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna1 = 0;                                                                                               % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Define Hip Flexor Motor Neuron Properties (Neuron 2).
Cm2 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm2 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er2 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R2 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am2 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm2 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm2 = R2;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah2 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh2 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh2 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna2 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh2_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna2 = 0;                                                                                               % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Define Knee Extensor Motor Neuron Properties (Neuron 3).
Cm3 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm3 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er3 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R3 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am3 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm3 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm3 = R3;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah3 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh3 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh3 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna3 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh3_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna3 = 0;                                                                                               % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Define Knee Flexor Motor Neuron Properties (Neuron 4).
Cm4 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm4 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er4 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R4 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am4 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm4 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm4 = R4;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah4 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh4 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh4 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna4 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh4_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna4 = 0;                                                                                               % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Define Hip Feedback Neuron Properties (Neuron 5).
Cm5 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm5 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er5 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R5 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am5 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm5 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm5 = R5;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah5 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh5 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh5 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna5 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh5_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna5 = TwoNeuronCPGSubnetworkNaConductance(R5, Gm5, Am5, Sm5, dEm5, Ah5, Sh5, dEh5, dEna5);             % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Define Knee Feedback Neuron Properties (Neuron 6).
Cm6 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm6 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er6 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R6 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am6 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm6 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm6 = R6;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah6 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh6 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh6 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna6 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh6_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna6 = TwoNeuronCPGSubnetworkNaConductance(R6, Gm6, Am6, Sm6, dEm6, Ah6, Sh6, dEh6, dEna6);             % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Define Hip Bistable Interneuron Properties (Neuron 7).
Cm7 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm7 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er7 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R7 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am7 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm7 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm7 = R7;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah7 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh7 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh7 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna7 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh7_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna7 = TwoNeuronCPGSubnetworkNaConductance(R7, Gm7, Am7, Sm7, dEm7, Ah7, Sh7, dEh7, dEna7);             % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Define Knee Bistable Intereuron Properties (Neuron 8).
Cm8 = 5e-9;                                                                                             % [C] Membrane Capacitance.
Gm8 = 1e-6;                                                                                             % [S] Membrane Conductance.
Er8 = -60e-3;                                                                                           % [V] Membrane Resting (Equilibrium) Potential.
R8 = 20e-3;                                                                                             % [V] Biphasic Equilibrium Voltage Range.
Am8 = 1;                                                                                                % [-] Sodium Channel Activation Parameter A.
Sm8 = -50;                                                                                              % [-] Sodium Channel Activation Parametter S.
dEm8 = R8;                                                                                              % [V] Sodium Channel Activation Reversal Potential w.r.t. Equilibrium Potential.
Ah8 = 0.5;                                                                                              % [-] Sodium Channel Deactivation Parameter A.
Sh8 = 50;                                                                                               % [-] Sodium Channel Deactivation Parameter S.
dEh8 = 0;                                                                                               % [V] Sodium Channel Deactivation Reversal Potential  w.r.t. Equilibrium Potential.
dEna8 = 110e-3;                                                                                         % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.
tauh8_max = 0.300;                                                                                      % [s] Maximum Sodium Channel Deactivation Time Constant.
Gna8 = TwoNeuronCPGSubnetworkNaConductance(R8, Gm8, Am8, Sm8, dEm8, Ah8, Sh8, dEh8, dEna8);             % [S] Sodium Channel Conductance.  (A zero value means that sodium channel currents will not be applied to this neuron.)

% Store the neuron properties into arrays.
Cms = [Cm1; Cm2; Cm3; Cm4; Cm5; Cm6; Cm7; Cm8];
Gms = [Gm1; Gm2; Gm3; Gm4; Gm5; Gm6; Gm7; Gm8];
Ers = [Er1; Er2; Er3; Er4; Er5; Er6; Er7; Er8];
Rs = [R1; R2; R3; R4; R5; R6; R7; R8]; Rs = repmat(Rs', [num_neurons, 1]);
Ams = [Am1; Am2; Am3; Am4; Am5; Am6; Am7; Am8];
Sms = [Sm1; Sm2; Sm3; Sm4; Sm5; Sm6; Sm7; Sm8];
dEms = [dEm1; dEm2; dEm3; dEm4; dEm5; dEm6; dEm7; dEm8];
Ahs = [Ah1; Ah2; Ah3; Ah4; Ah5; Ah6; Ah7; Ah8];
Shs = [Sh1; Sh2; Sh3; Sh4; Sh5; Sh6; Sh7; Sh8];
dEhs = [dEh1; dEh2; dEh3; dEh4; dEh5; dEh6; dEh7; dEh8];
dEnas = [dEna1; dEna2; dEna3; dEna4; dEna5; dEna6; dEna7; dEna8];
tauh_maxs = [tauh1_max; tauh2_max; tauh3_max; tauh4_max; tauh5_max; tauh6_max; tauh7_max; tauh8_max];
Gnas = [Gna1; Gna2; Gna3; Gna4; Gna5; Gna6; Gna7; Gna8];


%% Define Applied Current Magnitudes.

% Compute the necessary applied current magnitudes.
Iapp1 = 0;                     % [A] Applied Current.
Iapp2 = Gm2*R2;                % [A] Applied Current.
Iapp3 = 0;                     % [A] Applied Current.
Iapp4 = Gm4*R4;                % [A] Applied Current.
Iapp5 = Gm5*R5;                % [A] Applied Current.
Iapp6 = Gm6*R6;                % [A] Applied Current.
Iapp7 = Gm7*R7;                % [A] Applied Current.
Iapp8 = Gm8*R8;                % [A] Applied Current.     


%% Define Synapse Properties.

% Define the Bistable CPG subnetwork bifurcation parameters.
delta57 = -0.01e-3;          % [V] Voltage Difference Between Inhibited Neuron's Equilibrium Potential & the Presynaptic Threshold.
delta68 = -0.01e-3;          % [V] Voltage Difference Between Inhibited Neuron's Equilibrium Potential & the Presynaptic Threshold.

% Define synapse reversal potentials.
dEsyn61 = 2*R1;              % [V] Synapse Reversal Potential.
dEsyn64 = -2*R4;             % [V] Synapse Reversal Potential.
dEsyn52 = -2*R2;             % [V] Synapse Reversal Potential.
dEsyn53 = 2*R3;              % [V] Synapse Reversal Potential.
dEsyn57 = -40e-3;            % [V] Synapse Reversal Potential.
dEsyn75 = -40e-3;            % [V] Synapse Reversal Potential.
dEsyn68 = -40e-3;            % [V] Synapse Reversal Potential.
dEsyn86 = -40e-3;            % [V] Synapse Reversal Potential.

% Compute the synapse conductances.
gsyn61_max = Gm1*R1/(dEsyn61 - R1);
gsyn64_max = -Iapp4/dEsyn64;
gsyn52_max = -Iapp2/dEsyn52;
gsyn53_max = Gm3*R3/(dEsyn53 - R3);
gsyn57_max = TwoNeuronCPGSubnetworkSynConductance(delta57, dEsyn57, Am5, Sm5, dEm5, Ah5, Sh5, dEh5, dEna5, Gna5);
gsyn75_max = gsyn57_max;
gsyn68_max = TwoNeuronCPGSubnetworkSynConductance(delta68, dEsyn68, Am6, Sm6, dEm6, Ah6, Sh6, dEh6, dEna6, Gna6);
gsyn86_max = gsyn68_max;

% Store the synapse properties into arrays.
dEsyns = zeros(num_neurons, num_neurons);
dEsyns(1, 6) = dEsyn61;
dEsyns(4, 6) = dEsyn64;
dEsyns(2, 5) = dEsyn52;
dEsyns(3, 5) = dEsyn53;
dEsyns(7, 5) = dEsyn57;
dEsyns(5, 7) = dEsyn75;
dEsyns(8, 6) = dEsyn68;
dEsyns(6, 8) = dEsyn86;

gsyn_maxs = zeros(num_neurons, num_neurons);
gsyn_maxs(1, 6) = gsyn61_max;
gsyn_maxs(4, 6) = gsyn64_max;
gsyn_maxs(2, 5) = gsyn52_max;
gsyn_maxs(3, 5) = gsyn53_max;
gsyn_maxs(7, 5) = gsyn57_max;
gsyn_maxs(5, 7) = gsyn75_max;
gsyn_maxs(8, 6) = gsyn68_max;
gsyn_maxs(6, 8) = gsyn86_max;


%% Define Simulation Properties.

% Set the simulation time.
tf = 1;         % [s] Simulation Duration.
dt = 1e-3;      % [s] Simulation Time Step.

% Compute the simulation time vector.
ts = 0:dt:tf;

% Compute the number of time steps.
num_timesteps = length(ts);

% Set the network initial conditions.
Us0 = zeros(num_neurons, 1);
hs0 = zeros(num_neurons, 1);

% Define the number of cycles.
num_cycles = 5;

% Define the applied currents over time.
Iapp1s = Iapp1*ones(1, num_timesteps);
Iapp2s = Iapp2*ones(1, num_timesteps);
Iapp3s = Iapp3*ones(1, num_timesteps);
Iapp4s = Iapp4*ones(1, num_timesteps);

Iapp5s = zeros(1, num_timesteps); Iapp5s(1, 1) = Iapp5;
Iapp6s = zeros(1, num_timesteps); Iapp6s(1, 1) = Iapp6;
Iapp7s = zeros(1, num_timesteps);
Iapp8s = zeros(1, num_timesteps);

% Iapp5s = zeros(1, num_timesteps);
% Iapp6s = zeros(1, num_timesteps);
% Iapp7s = zeros(1, num_timesteps); Iapp7s(1, 1) = Iapp7;
% Iapp8s = zeros(1, num_timesteps); Iapp8s(1, 1) = Iapp8;

% Store the applied currents into arrays.
Iapps = [Iapp1s; Iapp2s; Iapp3s; Iapp4s; Iapp5s; Iapp6s; Iapp7s; Iapp8s];


%% Simulate the Network

% Simulate the network.
[ts, Us, hs, dUs, dhs, Gsyns, Ileaks, Isyns, Inas, Itotals, minfs, hinfs, tauhs] = SimulateNetwork(Us0, hs0, Gms, Cms, Rs, gsyn_maxs, dEsyns, Ams, Sms, dEms, Ahs, Shs, dEhs, tauh_maxs, Gnas, dEnas, Iapps, tf, dt);


%% Plot the Motor Neuron Subnetwork States vs Time.

% Plot the network membranve voltage vs time.
fig_U = figure('Color', 'w', 'Name', 'Motor Neuron Subnetwork States vs Time');

subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage, $U$ [V]', 'Interpreter', 'Latex'), title('Membrane Voltage vs Time')
plot(ts, Us(1:4, :), '-', 'Linewidth', 3)
legend({'(1) Hip Ext MN', '(2) Hip Flx MN', '(3) Knee Ext MN', '(4) Knee Flx MN'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage Derivative, $\dot{U}$ [V/s]', 'Interpreter', 'Latex'), title('Membrane Voltage Derivative vs Time')
plot(ts(1:end-1), dUs(1:4, 1:end-1), '-', 'Linewidth', 3)
legend({'(1) Hip Ext MN', '(2) Hip Flx MN', '(3) Knee Ext MN', '(4) Knee Flx MN'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')


%% Plot the Bistable Subnetwork States Over Time.

% Plot the bistable hip subnetwork states over time.
fig_HipBistableOverTime = figure('Color', 'w', 'Name', 'Bistable Hip Subnetwork States vs Time');

subplot(2, 2, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage, $U$ [V]', 'Interpreter', 'Latex'), title('Membrane Voltage vs Time')
plot(ts, Us(5, :), '-', 'Linewidth', 3)
plot(ts, Us(7, :), '-', 'Linewidth', 3)
legend({'(5) Hip Feedback', '(7) Hip Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

subplot(2, 2, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Sodium Channel Deactivation Parameter, $h$ [-]', 'Interpreter', 'Latex'), title('Sodium Channel Deactivation Parameter vs Time')
plot(ts, hs(5, :), '-', 'Linewidth', 3)
plot(ts, hs(7, :), '-', 'Linewidth', 3)
legend({'(5) Hip Feedback', '(7) Hip Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

subplot(2, 2, 3), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage Derivative, $\dot{U}$ [V/s]', 'Interpreter', 'Latex'), title('Membrane Voltage Derivative vs Time')
plot(ts(1:end-1), dUs(5, 1:end-1), '-', 'Linewidth', 3)
plot(ts(1:end-1), dUs(7, 1:end-1), '-', 'Linewidth', 3)
legend({'(5) Hip Feedback', '(7) Hip Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

subplot(2, 2, 4), hold on, grid on, xlabel('Time [s]'), ylabel('Sodium Channel Deactivation Parameter Derivative, $\dot{h}$ [-]', 'Interpreter', 'Latex'), title('Sodium Channel Deactivation Parameter Derivative vs Time')
plot(ts(1:end-1), dhs(5, 1:end-1), '-', 'Linewidth', 3)
plot(ts(1:end-1), dhs(7, 1:end-1), '-', 'Linewidth', 3)
legend({'(5) Hip Feedback', '(7) Hip Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

% Plot the bistable knee subnetwork states over time.
fig_KneeBistableOverTime = figure('Color', 'w', 'Name', 'Bistable Knee Subnetwork States vs Time');

subplot(2, 2, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage, $U$ [V]', 'Interpreter', 'Latex'), title('Membrane Voltage vs Time')
plot(ts, Us(6, :), '-', 'Linewidth', 3)
plot(ts, Us(8, :), '-', 'Linewidth', 3)
legend({'(6) Knee Feedback', '(8) Knee Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

subplot(2, 2, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Sodium Channel Deactivation Parameter, $h$ [-]', 'Interpreter', 'Latex'), title('Sodium Channel Deactivation Parameter vs Time')
plot(ts, hs(6, :), '-', 'Linewidth', 3)
plot(ts, hs(8, :), '-', 'Linewidth', 3)
legend({'(6) Knee Feedback', '(8) Knee Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

subplot(2, 2, 3), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage Derivative, $\dot{U}$ [V/s]', 'Interpreter', 'Latex'), title('Membrane Voltage Derivative vs Time')
plot(ts(1:end-1), dUs(6, 1:end-1), '-', 'Linewidth', 3)
plot(ts(1:end-1), dUs(8, 1:end-1), '-', 'Linewidth', 3)
legend({'(6) Knee Feedback', '(8) Knee Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

subplot(2, 2, 4), hold on, grid on, xlabel('Time [s]'), ylabel('Sodium Channel Deactivation Parameter Derivative, $\dot{h}$ [-]', 'Interpreter', 'Latex'), title('Sodium Channel Deactivation Parameter Derivative vs Time')
plot(ts(1:end-1), dhs(6, 1:end-1), '-', 'Linewidth', 3)
plot(ts(1:end-1), dhs(8, 1:end-1), '-', 'Linewidth', 3)
legend({'(6) Knee Feedback', '(8) Knee Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')


%% Plot the Bistable Subnetwork State Space Trajectories.

% Plot the bistable hip subnetwork state space trajectories.
fig_HipBistableStateTrajectory = figure('Color', 'w', 'Name', 'Bistable Hip Subnetwork State Space Trajectories');
hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage, $U$ [V]', 'Interpreter', 'Latex'), title('Membrane Voltage vs Time')
plot(Us(5, :), hs(5, :), '-', 'Linewidth', 3)
plot(Us(7, :), hs(5, :), '-', 'Linewidth', 3)
legend({'(5) Hip Feedback', '(7) Hip Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')

% Plot the bistable knee subnetwork state space trajectories.
fig_KneeBistableStateTrajectory = figure('Color', 'w', 'Name', 'Bistable Knee Subnetwork State Space Trajectories');
hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage, $U$ [V]', 'Interpreter', 'Latex'), title('Membrane Voltage vs Time')
plot(Us(6, :), hs(6, :), '-', 'Linewidth', 3)
plot(Us(8, :), hs(8, :), '-', 'Linewidth', 3)
legend({'(6) Hip Feedback', '(8) Hip Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')



%% Plot the Network Currents

fig_I = figure('Color', 'w', 'Name', 'Network Currents vs Time');
subplot(5, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Leak Current, $I_{leak}$ [A]', 'Interpreter', 'Latex'), title('Leak Current vs Time'), plot(ts(1:end-1), Ileaks(:, 1:end-1), '-', 'Linewidth', 3)
subplot(5, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Synaptic Current, $I_{syn}$ [A]', 'Interpreter', 'Latex'), title('Synaptic Current vs Time'), plot(ts(1:end-1), Isyns(:, 1:end-1), '-', 'Linewidth', 3)
subplot(5, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Sodium Channel Current, $I_{na}$ [A]', 'Interpreter', 'Latex'), title('Sodium Channel Current vs Time'), plot(ts(1:end-1), Inas(:, 1:end-1), '-', 'Linewidth', 3)
subplot(5, 1, 4), hold on, grid on, xlabel('Time [s]'), ylabel('Applied Current, $I_{app}$ [A]', 'Interpreter', 'Latex'), title('Applied Current vs Time'), plot(ts(1:end-1), Iapps(:, 1:end-1), '-', 'Linewidth', 3)
subplot(5, 1, 5), hold on, grid on, xlabel('Time [s]'), ylabel('Total Current, $I_{total}$ [A]', 'Interpreter', 'Latex'), title('Total Current vs Time'), plot(ts(1:end-1), Itotals(:, 1:end-1), '-', 'Linewidth', 3)
legend({'Hip Ext MN', 'Hip Flx MN', 'Knee Ext MN', 'Knee Flx MN', 'Hip Feedback', 'Knee Feedback', 'Hip Bistable', 'Knee Bistable'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal')



