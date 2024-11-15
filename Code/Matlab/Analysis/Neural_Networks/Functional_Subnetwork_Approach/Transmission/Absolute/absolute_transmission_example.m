%% Absolute Transmission Subnetwork Example.

% Clear Everything.
clear, close( 'all' ), clc


%% Define Simulation Parameters.

% Define the save and load directories.
save_directory = '.\Save';                       	% [str] Save Directory.
load_directory = '.\Load';                         	% [str] Load Directory.

% Define the level of verbosity.
verbose_flag = true;                             	% [T/F] Printing Flag.

% Define the undetected option.
undetected_option = 'error';                        % [str] Undetected Option.

% Define the network integration step size.
network_dt = 1e-3;                                  % [s] Simulation Timestep.
% network_dt = 1e-4;                            	% [s] Simulation Timestep.

% Define the network simulation duration.
network_tf = 0.5;                                 	% [s] Simulation Duration.
% network_tf = 1;                                 	% [s] Simulation Duration.
% network_tf = 3;                                 	% [s] Simulation Duration.

% Define the integration method.
integration_method = 'RK4';                         % [str] Integration Method (Either FE for Forward Euler or RK4 for Fourth Order Runge-Kutta).


%% Define Absolute Transmission Subnetwork Parameters.

% Define the transmission subnetwork design parameters.
c = 1.0;                                            % [-] Absolute Transmission Subnetwork Gain.
R1 = 20e-3;                                         % [V] Maximum Membrane Voltage (Neuron 1).
Gm1 = 1e-6;                                         % [S] Membrane Conductance (Neuron 1).
Gm2 = 1e-6;                                       	% [S] Membrane Conductance (Neuron 2).
Cm1 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 1).
Cm2 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 2).

% Store the transmission subnetwork design parameters in a cell.
transmission_parameters = { c, R1, Gm1, Gm2, Cm1, Cm2 };

% Define the encoding scheme.
encoding_scheme = 'absolute';


%% Define the Absolute Transmission Subnetwork Input Current Parameters.

% Define the current identification properties.
input_current_ID = 1;                               % [#] Input Current ID.
input_current_name = 'Applied Current 1';           % [str] Input Current Name.
input_current_to_neuron_ID = 1;                     % [#] Neuron ID to Which Input Current is Applied.

% Compute the number of simulation timesteps.
n_timesteps = floor( network_tf/network_dt ) + 1;   % [#] Number of Simulation Timesteps.

% Construct the simulation times associated with the input currents.
ts = ( 0:network_dt:network_tf )';                 	% [s] Simulation Times.

% Define the current magnitudes.
Ias1_mag = R1*Gm1;                                  % [A] Applied Current Magnitude.

% Define the magnitudes of the applied current input.
Ias1 = Ias1_mag*ones( n_timesteps, 1 );             % [A] Applied Currents.


%% Create Absolute Transmission Subnetwork.

% Create an instance of the netwo5rk class.
network = network_class( network_dt, network_tf );

% Create a transmission subnetwork.
[ c, Gnas, R2, dEs21, gs21, Ia2, neurons, synapses, neuron_manager, synapse_manager, network ] = network.create_transmission_subnetwork( transmission_parameters, encoding_scheme, network.neuron_manager, network.synapse_manager, network.applied_current_manager, true, true, false, undetected_option );

% Create the input applied current.
[ ~, ~, ~, network.applied_current_manager ] = network.applied_current_manager.create_applied_current( input_current_ID, input_current_name, input_current_to_neuron_ID, ts, Ias1, true, network.applied_current_manager.applied_currents, true, false, network.applied_current_manager.array_utilities );


%% Print Absolute Transmission Subnetwork Information.

% Print transmission subnetwork information.
network.print( network.neuron_manager, network.synapse_manager, network.applied_current_manager, verbose_flag );


%% Compute Absolute Transmission Numerical Stability Analysis Parameters.

% Define the property retrieval settings.
as_matrix_flag = true;

% Retrieve properties from the existing network.
Cms = network.neuron_manager.get_neuron_property( 'all', 'Cm', as_matrix_flag, network.neuron_manager.neurons, undetected_option );         % [F] Membrane Capacitance.
Gms = network.neuron_manager.get_neuron_property( 'all', 'Gm', as_matrix_flag, network.neuron_manager.neurons, undetected_option );         % [S] Membrane Conductance.
Rs = network.neuron_manager.get_neuron_property( 'all', 'R', as_matrix_flag, network.neuron_manager.neurons, undetected_option );           % [V] Maximum Membrane Voltage.
gs = network.get_gs( 'all', network.neuron_manager, network.synapse_manager );                                                              % [S] Synaptic Conductance.
dEs = network.get_dEs( 'all', network.neuron_manager, network.synapse_manager );                                                            % [V] Synaptic Reversal Potential.
Us = zeros( 1, network.neuron_manager.num_neurons );                                                                                        % [V] Membrane Voltage.

% Define the stability analysis timestep seed.
dt0 = 1e-6;                                                                                                                                 % [s] Stability Analysis Time Step Seed.

% Compute the maximum RK4 step size and condition number.
[ As, dts, condition_numbers ] = network.RK4_stability_analysis( Cms, Gms, Rs, gs, dEs, Us, dt0, network.neuron_manager, network.synapse_manager, undetected_option, network.network_utilities );

% Print out the stability information.
network.numerical_method_utilities.print_numerical_stability_info( As, dts, network_dt, condition_numbers );


%% Simulate the Absolute Transmission Subnetwork.

% Set additional simulation properties.
filter_disabled_flag = true;                % [T/F] Filter Disabled Flag.
set_flag = true;                            % [T/F] Set Flag.
process_option = 'None';                    % [str] Process Option.
undetected_option = 'Ignore';               % [str] Undetected Option.

% Start the timer.
tic

% Simulate the network.
[ ts, Us, hs, dUs, dhs, Gs, I_leaks, I_syns, I_nas, I_apps, I_totals, m_infs, h_infs, tauhs, neurons, synapses, neuron_manager, synapse_manager, network ] = network.compute_simulation( network_dt, network_tf, integration_method, network.neuron_manager, network.synapse_manager, network.applied_current_manager, network.applied_voltage_manager, filter_disabled_flag, set_flag, process_option, undetected_option, network.network_utilities );

% End the timer.
toc


%% Plot the Absolute Transmission Subnetwork Results.

% Compute the decoded output.
Us_decoded = Us*( 10^3 );

% Retrieve the neuron IDs.
neuron_IDs = network.neuron_manager.get_all_neuron_IDs( network.neuron_manager.neurons );

% Plot the network currents over time.
fig_network_currents = network.network_utilities.plot_network_currents( ts, I_leaks, I_syns, I_nas, I_apps, I_totals, neuron_IDs );

% Plot the network states over time.
fig_network_states = network.network_utilities.plot_network_states( ts, Us, hs, neuron_IDs );

% Plot the absolute network decoding over time.
fig_network_decoding = figure( 'Color', 'w', 'Name', 'Absolute Transmission Decoding vs Time' ); hold on, grid on, xlabel( 'Time, t [s]' ), ylabel( 'Network Decoding [-]' ), title( 'Absolute Transmission Decoding vs Time' )
plot( ts, Us_decoded( 1, : ), '-', 'Linewidth', 3 )
plot( ts, Us_decoded( 2, : ), '-', 'Linewidth', 3 )
legend( 'Input', 'Output' )
saveas( fig_network_decoding, [ save_directory, '\', 'absolute_transmission_decoding_example' ] )

% Plot the absolute network dynamic decoding example.
fig_network_decoding = figure( 'Color', 'w', 'Name', 'Absolute Transmission Dynamic Decoding Example' ); hold on, grid on, xlabel( 'Input [-]' ), ylabel( 'Output [-]' ), title( 'Absolute Transmission Dynamic Decoding Example' )
plot( Us_decoded( 1, : ), Us_decoded( 2, : ), '-', 'Linewidth', 3 )
saveas( fig_network_decoding, [ save_directory, '\', 'absolute_transmission_dynamic_decoding_example' ] )

% Animate the network states over time.
fig_network_animation = network.network_utilities.animate_network_states( Us, hs, neuron_IDs );

