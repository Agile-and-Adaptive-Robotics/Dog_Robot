%% CPG Model

% Clear Everything.
clear, close('all'), clc

%% Simulate the System.

% Define the simulation properties.
% V0 = [1; 1; 1; 1];
% h0 = [0; 0; 0; 0];
V0 = [-50e-3; -50e-3; -50e-3; -50e-3];
h0 = [0.67; 0.67; 0.67; 0.67];
x0 = [V0; h0];
% tspan = [0 10];
tspan = [0 0.1];

% Simulate the system.
[ts, ys] = ode45(@HalfCenter, tspan, x0);

%% Plot the System Response.

% Define the number of neurons in the simulation.
num_neurons = 4;

% Create figures to store the simulation data.
fig_neuron_voltages_vs_time = figure('color', 'w', 'name', 'CPG Neurons Membrane Voltage vs Time');
fig_neuron_Na_deactivation_vs_time = figure('color', 'w', 'name', 'CPG Neurons Na Ch. Deactivation vs Time');
fig_neuron_state_trajectories = figure('color', 'w', 'name', 'CPG Neuron State Trajectories');

% Plot the states of the CPG neurons.
for k = 1:num_neurons           % Iterate through each of the neurons...

    % Plot the CPG neuron membrane voltage over time.
    figure(fig_neuron_voltages_vs_time), subplot(2, 2, k), hold on, grid on, xlabel('Time [s]', 'Interpreter', 'latex'), ylabel('Neuron Membrane Voltage, $V$ [V]', 'Interpreter', 'latex'), title(sprintf('Neuron %0.0f: Membrane Voltage, $V$ [V] vs Time', k), 'Interpreter', 'latex'), plot(ts, ys(:, k), '-', 'Linewidth', 3)
    figure(fig_neuron_Na_deactivation_vs_time), subplot(2, 2, k), hold on, grid on, xlabel('Time [s]', 'Interpreter', 'latex'), ylabel('Neuron Na Ch. Deactivation, $h$ [-]', 'Interpreter', 'latex'), title(sprintf('Neuron %0.0f: Na Ch. Deactivation, $h$ [-] vs Time', k), 'Interpreter', 'latex'), plot(ts, ys(:, k + 4), '-', 'Linewidth', 3)


    % Plot the CPG neuron state trajectory over time.
    figure(fig_neuron_state_trajectories), subplot(2, 2, k), hold on, grid on, xlabel('Membrane Voltage, $V$ [V]', 'Interpreter', 'latex'), ylabel('Na Ch. Deactivation, $h$, [-]', 'Interpreter', 'latex'), title(sprintf('Neuron %0.0f: State Trajectory', k), 'Interpreter', 'latex'), plot(ys(:, k), ys(:, k + 4), '-', 'Linewidth', 3)
    
end

% legend('Membrane Voltage, V', 'Na Ch. Deactivation, h', 'Location', 'South', 'Orientation', 'Horizontal')

%% Define the System Dynamics.

function dxdt = HalfCenter(t, x)

% Retrieve the components of the input vector.
Vs = x(1:4)'; hs = x(5:end)';

% Define the input current.
Iapps = [(10*(10^(-9))) 0 0 0]; 
% Iapps = [0 0 0 0]; 

% Define membrane properties.
Cms = (2.5e-9)*ones(size(Vs));
Gms = (500e-9)*ones(size(Vs));
Ers = (-50e-3)*ones(size(Vs));

% Define synapse properties.
Elos = (-65e-3)*ones(size(Vs));
Ehis = (-20e-3)*ones(size(Vs));
gmaxs = (0.5e-6)*ones(size(Vs));
Ess = (-10e-3)*ones(size(Vs));

% Define sodium channel properties.
Ams = 1*ones(size(Vs));
Sms = 0.046*ones(size(Vs));
Ems = (-40e-3)*ones(size(Vs));
Ahs = 0.5*ones(size(Vs));
Shs = -0.046*ones(size(Vs));
Ehs = (-60e-3)*ones(size(Vs));
Gnas = (1e-6)*ones(size(Vs));
Enas = (-50e-3)*ones(size(Vs));
tauhmaxs = 0.005*ones(size(Vs));

% Compute the synapse conductance.
Gss = gmaxs.*min( max( (circshift(Vs, 1) - Elos)./(Ehis - Elos), 0), 1);

% Compute the steady state sodium channel activation and deactivation parameters.
minfs = 1./(1 + Ams.*exp(Sms.*(Vs - Ems))); hinfs = 1./(1 + Ahs.*exp(Shs.*(Vs - Ehs)));

% Compute the sodium channel deactivation time constant.
tauhs = tauhmaxs.*hinfs.*sqrt(Ahs.*exp(Shs.*(Vs - Ehs)));

% Compute the leak current.
Ileaks = Gms.*(Ers - Vs);

% Compute the synaptic current.
Isyns = Gss.*(Ess - Vs);

% Compute the sodium current.
Inas = Gnas.*minfs.*hs.*(Enas - Vs);

% Compute the total current.
Itotals = Ileaks + Isyns + Inas + Iapps;

% Compute the membrane voltage derivative.
dVdts = Itotals./Cms;

% Compute the sodium channel deactivation derivative.
dhdts = (hinfs - hs)./tauhs;

% Compute the state derivative.
dxdt = [dVdts'; dhdts'];

end





