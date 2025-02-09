%% CPG Model

% Clear Everything.
clear, close('all'), clc

%% Define the Network Properties.

% Define the number of neurons in the simulation.
num_neurons = 2;

% Define the voltage range for the CPG to oscillator over.
R = 20e-3;                                          % [V] Biphasic Equilibrium Voltage Range.

% Define membrane properties.
Cm = 5e-9;                   % [C] Membrane Capacitance.
Gm = 1e-6;                   % [S] Membrane Conductance.
Er = -60e-3;                 % [V] Membrane Resting (Equilibrium) Potential.

% Define synapse properties.
Elo = Er;                    % [V] Presynaptic Threshold.
Ehi = Elo + R;               % [V] Presynaptic Saturation Level.
Es = -100e-3;                % [V] Synaptic Equilibrium Potential.
Es_tilde = Es - Er;         % [V] Synaptic Equilibrium Potential With Respect to
delta = 0.01e-3;                 % [V] Voltage Difference Between Inhibited Neuron's Equilibrium Potential & the Presynaptic Threshold.
% delta = 7.0e-3;                 % [V] Voltage Difference Between Inhibited Neuron's Equilibrium Potential & the Presynaptic Threshold.

% Define sodium channel activation properties.
Am = 1;
Sm = -50;
Em = Ehi;
Em_tilde = Em - Er;

% Define sodium channel deactivation properties.
Ah = 0.5;
Sh = 50;
Eh = Elo;
Eh_tilde = Eh - Er;

% Define the steady state sodium channel activation & deactivation parameters intermediate calculations.
minfterm_func = @(U) Am.*exp(-Sm.*(Em_tilde - U));
hinfterm_func = @(U) Ah.*exp(-Sh.*(Eh_tilde - U));

% Define the steady state sodium channel activation & deactivation parameters.
minf_func = @(U) 1./(1 + Am.*exp(-Sm.*(Em_tilde - U)));
hinf_func = @(U) 1./(1 + Ah.*exp(-Sh.*(Eh_tilde - U)));

% Define the sodium channel reversal potential.
Ena = 50e-3;                % [V] Sodium Channel Reversal Potential.
Ena_tilde = Ena - Er;       % [V] Sodium Channel Reversal Potential With Respect to the Resting Potential.

% Compute the sodium channel conductance.
Gna = (Gm*R)/(minf_func(R)*hinf_func(R)*(Ena_tilde - R));       % [S] Sodium Channel Conductance.

% Define the maximum sodium channel time constant.
tauhmax = 0.3;             % [s] Maximum Sodium Channel Time Constant.

% Define the sodium channel time constant intermediate calculation.
tauhterm1_func = @(U) Ah.*exp(-Sh.*(Eh_tilde - U));
tauhterm2_func = @(U) sqrt(Ah.*exp(-Sh.*(Eh_tilde - U)));

% Define the sodium channel time constant.
tauh_func = @(U) tauhmax*hinf_func(U).*sqrt(Ah.*exp(-Sh.*(Eh_tilde - U)));       % [s] Sodium Channel Time Constant.

% Compute the maximum synaptic conductance.
gsynmax = (-delta*(10^(-6)) - delta*Gna*minf_func(delta)*hinf_func(delta) + Gna*minf_func(delta)*hinf_func(delta)*Ena_tilde)/(delta - Es_tilde);

% Define a function to compute the synaptic conductance.
gsyn_func = @(Upre) gsynmax*min( max( Upre/R, 0 ), 1 );

% Define a function to compute the leak current.
Ileak_func = @(U) -Gm*U;

% Define a function to compute the synaptic current.
Isyn_func = @(Upre, Upost) gsyn_func(Upre).*( Es_tilde - Upost );

% Define a function to compute the sodium channel current.
Ina_func = @(U, h) Gna*minf_func(U).*h.*( Ena_tilde - U );

% Define a function to compute the applied current.
Iapp_func = @(t) zeros(num_neurons, 1);


Us = (-4*R):0.0001:(6*R);
gsyns = gsyn_func(Us);
minfterms = minfterm_func(Us);
hinfterms = hinfterm_func(Us);
minfs = minf_func(Us);
hinfs = hinf_func(Us);
tauhterms1 = tauhterm1_func(Us);
tauhterms2 = tauhterm2_func(Us);
tauhs = tauh_func(Us);
Ileaks = Ileak_func(Us);
Isyns = Isyn_func(Us, Us);

[USpre, USpost] = meshgrid(Us, Us);

Isyns_surf = Isyn_func(USpre, USpost);

figure, hold on, grid on, xlabel('U [V]'), ylabel('$m_{\infty}$ [-] Intermediate Calculation', 'Interpreter', 'latex'), title('$m_{\infty}$ vs U (Intermediate Calculation)', 'Interpreter', 'latex'), plot(Us, minfterms)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$h_{\infty}$ [-] Intermediate Calculation', 'Interpreter', 'latex'), title('$h_{\infty}$ vs U (Intermediate Calculation)', 'Interpreter', 'latex'), plot(Us, hinfterms)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$m_{\infty}$ [-]', 'Interpreter', 'latex'), title('$m_{\infty}$ vs U', 'Interpreter', 'latex'), plot(Us, minfs)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$h_{\infty}$ [-]', 'Interpreter', 'latex'), title('$h_{\infty}$ vs U', 'Interpreter', 'latex'), plot(Us, hinfs)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$\tau_{h}$ [s] Intermediate Calculation 1', 'Interpreter', 'latex'), title('$\tau_{h}$ [s] Intermediate Calculation 1 vs U', 'Interpreter', 'latex'), plot(Us, tauhterms1)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$\tau_{h}$ [s] Intermediate Calculation 2', 'Interpreter', 'latex'), title('$\tau_{h}$ [s] Intermediate Calculation 2 vs U', 'Interpreter', 'latex'), plot(Us, tauhterms2)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$\tau_{h}$ [s] Complete Calculation', 'Interpreter', 'latex'), title('$\tau_{h}$ [s] Full Calculation vs U', 'Interpreter', 'latex'), plot(Us, tauhs)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$I_{leak}$ [A]', 'Interpreter', 'latex'), title('$I_{leak}$ [s] vs U', 'Interpreter', 'latex'), plot(Us, Ileaks)
figure, hold on, grid on, xlabel('U [V]'), ylabel('$I_{syn}$ [A]', 'Interpreter', 'latex'), title('$I_{syn}$ [s] vs U', 'Interpreter', 'latex'), plot(Us, Isyns)
figure, hold on, grid on, rotate3d on, xlabel('$U_{pre}$ [V]'), ylabel('$U_{post}$ [V]', 'Interpreter', 'latex'), zlabel('$I_{syn}$ [A]', 'Interpreter', 'latex'), title('$I_{syn}$ vs $U_{pre}$ \& $U_{post}$', 'Interpreter', 'latex'), surf(USpre, USpost, Isyns_surf, 'Edgecolor', 'none')
plot3(Us, Us, Isyns, '-r', 'Linewidth', 3)

%% Simulate the System.


% Define simulation properties.
t_sim_duration = 5;         % [s] Simulation Duration.
dt_sim = 1e-3;              % [s] Simulation Time Step.
% dt_sim = 1e-2;              % [s] Simulation Time Step.

% Create the simulation time vector.
ts = 0:dt_sim:t_sim_duration;       % [s] Simulation Time Vector.

% Compute the number of simulation time steps.
num_steps = length(ts) - 1;

% Create variables to store the simulation data.
[Us, hs, dUs, dhs, Ileaks, Isyns, Inas, Iapps, Itotals, Gsyns, minfs, hinfs, tauhs] = deal( zeros(num_neurons, length(ts)) );

% Initialize the simulation.
% hs(:, 1) = hinf_func(Us(:, 1));
Iapps(1, 1) = 1e-9;
% Iapps(1, :) = 1e-9*ones(1, size(Iapps, 2));

% Run the Simulation.
for k = 1:num_steps             % Iterate through each of the simulation time steps...
    
    % Compute the steady state sodium channel activation and deactivation parameters.
    minfs(:, k) = minf_func(Us(:, k));
    hinfs(:, k) = hinf_func(Us(:, k));

    % Compute the sodium channel time constant.
    tauhs(:, k) = tauh_func(Us(:, k));
    
    % Compute the synaptic conductances.
    Gsyns(:, k) = gsyn_func(flipud(Us(:, k)));

    % Compute the leak currents.
    Ileaks(:, k) = Ileak_func(Us(:, k));
    
    % Compute the synaptic currents.
    Isyns(:, k) = Isyn_func(flipud(Us(:, k)), Us(:, k));
    
    % Compute the sodium channel currents.
    Inas(:, k) = Ina_func(Us(:, k), hs(:, k));
    
%     % Compute the applied currents.
%     Iapps(:, k) = Iapp_func(ts(k));
    
    % Compute the total current.
    Itotals(:, k) = Ileaks(:, k) + Isyns(:, k) + Inas(:, k) + Iapps(:, k);
%     Itotals(:, k) = Ileaks(:, k) + Iapps(:, k);
%     Itotals(:, k) = Ileaks(:, k) + Isyns(:, k) + Iapps(:, k);

    % Compute the membrane voltage derivative.
    dUs(:, k) = Itotals(:, k)/Cm;
    
    % Compute the sodium channel deactivation parameter derivative.
    dhs(:, k) = (hinfs(:, k) - hs(:, k))./tauhs(:, k);
    
    % Estimate the next states.
    Us(:, k + 1) = Us(:, k) + dt_sim*dUs(:, k);
    hs(:, k + 1) = hs(:, k) + dt_sim*dhs(:, k);

end


%% Plot the Simulation Results.

% Create a figure to store the system state space trajectory.
fig_state_trajectory = figure('color', 'w', 'name', 'Network State Space Trajectory'); hold on, grid on, xlabel('Membrane Voltage w.r.t Resting Potential, U [V]'), ylabel('Na Ch. Deactivation, h [-]'), title('Network State Space Trajectory')
plot(Us(1, :), hs(1, :), '-', 'Linewidth', 3)
plot(Us(2, :), hs(2, :), '-', 'Linewidth', 3)

% Create a figure to store the system states over time.
fig_states_vs_time = figure('color', 'w', 'name', 'Network States Over Time');

subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage w.r.t. Resting Voltage, $U$ [V]', 'Interpreter', 'latex'), title('Membrane Voltage w.r.t. Resting Voltage, $U$ [V] vs Time [s]', 'Interpreter', 'latex')
plot(ts, Us, '-', 'Linewidth', 3)

subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Sodium Channel Deactivation Parameter $h$, [-]', 'Interpreter', 'latex'), title('Sodium Channel Deactivation Parameter $h$, [-] vs Time [s]', 'Interpreter', 'latex')
plot(ts, hs, '-', 'Linewidth', 3)

% Create a figure to store the system state derivatives over time.
fig_state_derivatives_vs_time = figure('color', 'w', 'name', 'Network State Derivatives Over Time');

subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Membrane Voltage w.r.t Resting Voltage Derivative, $\dot{U}$ [V/s]', 'Interpreter', 'latex'), title('Membrane Voltage w.r.t Resting Voltage Derivative, $\dot{U}$ [V/s] vs Time [s]', 'Interpreter', 'latex')
plot(ts, dUs, '-', 'Linewidth', 3)

subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Na Ch. Deactivation Parameter Derivatie, $\dot{h}$ [-/s]', 'Interpreter', 'latex'), title('Na Ch. Deactivation Parameter Derivatie, $\dot{h}$ [-/s] vs Time [s]', 'Interpreter', 'latex')
plot(ts, dhs, '-', 'Linewidth', 3)

% Create a figure to store the currents over time.
fig_currents_vs_time = figure('color', 'w', 'name', 'Network Currents Over Time');

subplot(5, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Leak Currents, $I_{leak}$ [A]', 'Interpreter', 'latex'), title('Leak Currents, $I_{leak}$ [A] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), Ileaks(:, 1:end-1), '-', 'Linewidth', 3) 

subplot(5, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Synaptic Currents, $I_{syn}$ [A]', 'Interpreter', 'latex'), title('Synaptic Currents, $I_{syn}$ [A] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), Isyns(:, 1:end-1), '-', 'Linewidth', 3) 

subplot(5, 1, 3), hold on, grid on, xlabel('Time [s]'), ylabel('Na Ch. Currents, $I_{Na}$ [A]', 'Interpreter', 'latex'), title('Na Ch. Currents, $I_{Na}$ [A] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), Inas(:, 1:end-1), '-', 'Linewidth', 3) 

subplot(5, 1, 4), hold on, grid on, xlabel('Time [s]'), ylabel('Applied Currents, $I_{app}$ [A]', 'Interpreter', 'latex'), title('Applied Currents, $I_{app}$ [A] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), Iapps(:, 1:end-1), '-', 'Linewidth', 3) 

subplot(5, 1, 5), hold on, grid on, xlabel('Time [s]'), ylabel('Total Currents, $I_{total}$ [A]', 'Interpreter', 'latex'), title('Total Currents, $I_{total}$ [A] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), Itotals(:, 1:end-1), '-', 'Linewidth', 3) 


% Create a figure to store the sodium channel activation & deactivation parameters over time.
fig_mhinf_vs_time = figure('color', 'w', 'name', 'Steady State Sodium Channel Activation & Deactivation Over Time');

subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Steady State Sodium Channel Activation, $m_{\infty}$ [-]', 'Interpreter', 'latex'), title('Steady State Sodium Channel Activation, $m_{\infty}$ [-] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), minfs(:, 1:end-1), '-', 'Linewidth', 3)

subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Steady State Sodium Channel Deactivation, $h_{\infty}$ [-]', 'Interpreter', 'latex'), title('Steady State Sodium Channel Deactivation, $h_{\infty}$ [-] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), hinfs(:, 1:end-1), '-', 'Linewidth', 3)


% Create a figure to store the sodium channel deactivation time constant over time.
fig_tauh_vs_time = figure('color', 'w', 'name', 'Sodium Channel Deactivation Time Constants Over Time'); hold on, grid on, xlabel('Time [s]'), ylabel('Na Ch. Deactivation Time Constant, $\tau_{h}$ [s]', 'Interpreter', 'latex'), title('Na Ch. Deactivation Time Constant, $\tau_{h}$ [s] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), tauhs(:, 1:end-1), '-', 'Linewidth', 3)

% Create a figure to store the synaptic conductances.
fig_gsyns_vs_time = figure('color', 'w', 'name', 'Synaptic Conductances Over Time'); hold on, grid on, xlabel('Time [s]'), ylabel('Synaptic Conductance, $G_{s,i}$ [S]', 'Interpreter', 'latex'), title('Synaptic Conductance, $G_{s,i}$ [S] vs Time [s]', 'Interpreter', 'latex')
plot(ts(1:end-1), Gsyns(:, 1:end-1), '-', 'Linewidth', 3)

% Make the network states over time figure active.
figure(fig_states_vs_time)

