%% Absolute Inversion Subnetwork Example.

% Clear Everything.
clear, close( 'all' ), clc


%% Define Simulation Parameters.

% Define the save and load directories.
save_directory = '.\Save';                        	% [str] Save Directory.
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

% Compute the number of simulation timesteps.
n_timesteps = floor( network_tf/network_dt ) + 1;   % [#] Number of Simulation Timesteps.

% Construct the simulation times associated with the input currents.
ts = ( 0:network_dt:network_tf )';                 	% [s] Simulation Times.

% Define the integration method.
integration_method = 'RK4';                         % [str] Integration Method (Either FE for Forward Euler or RK4 for Fourth Order Runge-Kutta).

% Define the encoding scheme.
encoding_scheme = 'absolute';


%% Define the Desired Inversion Subnetwork Parameters.

% Create an instance of the network utilities class.
network_utilities = network_utilities_class(  );

% Define the inversion subnetwork parameters.
c1 = 0.40e-9;                                       % [W] Design Constant 1.
c3 = 20e-9;                                         % [A] Design Constant 3.
delta = 1e-3;                                       % [V] Membrane Voltage Offset.

% Define the desired mapping operation.
f_desired = @( x, c2 ) network_utilities.compute_desired_inversion_sso( x, c1, c2, c3 );


%% Define the Encoding & Decoding Operations.

% Define the domain of the input and output signals.
x_max = 20;

% Define the encoding scheme.
f_encode = @( x ) x*( 10^( -3 ) );

% Define the decoding scheme.
f_decode = @( U ) U*( 10^3 );


%% Define Additional Absolute Inversion Design Subnetwork Parameters.

% Define the inversion subnetwork design parameters.
R1 = 20e-3;                                         % [V] Maximum Membrane Voltage (Neuron 1).
Gm1 = 1e-6;                                         % [S] Membrane Conductance (Neuron 1).
Gm2 = 1e-6;                                       	% [S] Membrane Conductance (Neuron 2).
% Cm1 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 1).
% Cm2 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 2).
Cm1 = 15e-9;                                     	% [F] Membrane Capacitance (Neuron 1).
Cm2 = 15e-9;                                      % [F] Membrane Capacitance (Neuron 2).

% Store the inversion subnetwork design parameters in a cell.
inversion_parameters = { c1, c3, delta, R1, Gm1, Gm2, Cm1, Cm2 };


%% Define the Desired Input Signal.

% Define the desired decoded input signal.
% xs_desired = 0*ones( n_timesteps, 1 );
xs_desired = x_max*ones( n_timesteps, 1 );

% Encode the input signal.
Us1_desired = f_encode( xs_desired );


%% Define the Absolute Inversion Subnetwork Input Current Parameters.

% Define the current identification properties.
input_current_ID = 1;                               % [#] Input Current ID.
input_current_name = 'Applied Current 1';           % [str] Input Current Name.
input_current_to_neuron_ID = 1;                     % [#] Neuron ID to Which Input Current is Applied.

% Define the magnitudes of the applied current input.
Ias1 = Us1_desired*Gm1;                           	% [A] Applied Currents.


%% Create Absolute Inversion Subnetwork.

% Create an instance of the network class.
network = network_class( network_dt, network_tf );

% Create a inversion subnetwork.
[ cs, Gnas, R2, dEs21, gs21, Ia2, ~, ~, ~, ~, ~, ~, network ] = network.create_inversion_subnetwork( inversion_parameters, encoding_scheme, network.neuron_manager, network.synapse_manager, network.applied_current_manager, true, true, false, undetected_option );

% Update the input current ID and name.
[ ~, network.applied_current_manager ] = network.applied_current_manager.set_applied_current_property( network.applied_current_manager.applied_currents( 1 ).ID, 2, 'ID', network.applied_current_manager.applied_currents, true );
[ ~, network.applied_current_manager ] = network.applied_current_manager.set_applied_current_property( network.applied_current_manager.applied_currents( 1 ).ID, { 'Applied Current 2' }, 'name', network.applied_current_manager.applied_currents, true );

% Create the input applied current.
[ ~, ~, ~, network.applied_current_manager ] = network.applied_current_manager.create_applied_current( input_current_ID, input_current_name, input_current_to_neuron_ID, ts, Ias1, true, network.applied_current_manager.applied_currents, true, false, network.applied_current_manager.array_utilities );

% Reverse the order of the applied currents in the applied current manager for cleanliness.
temporary_applied_current = network.applied_current_manager.applied_currents( 1 );
network.applied_current_manager.applied_currents( 1 ) = network.applied_current_manager.applied_currents( 2 );
network.applied_current_manager.applied_currents( 2 ) = temporary_applied_current;


%% Print Absolute Inversion Subnetwork Parameters.

% Print inversion subnetwork information.
network.print( network.neuron_manager, network.synapse_manager, network.applied_current_manager, verbose_flag );


%% Compute Derived Encoding & Decoding Parameters.

% Retrieve the second gain variable.
c2 = cs( 2 );

% Compute the maximum output signal.
y_max = f_desired( 0, c2 );


%% Compute Absolute Inversion Numerical Stability Analysis Parameters.

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


%% Print the Numerical Stability Information.

% Print out the stability information.
network.numerical_method_utilities.print_numerical_stability_info( As, dts, network_dt, condition_numbers );


%% Simulate the Absolute Inversion Subnetwork.

% Set additional simulation properties.
filter_disabled_flag = true;                % [T/F] Filter Disabled Flag.
set_flag = true;                            % [T/F] Set Flag.
process_option = 'None';                    % [str] Process Option.
undetected_option = 'Ignore';               % [str] Undetected Option.

% Start the timer.
tic

% Simulate the network.
[ ts, Us, hs, dUs, dhs, Gs, I_leaks, I_syns, I_nas, I_apps, I_totals, m_infs, h_infs, tauhs, ~, ~, ~, ~, network ] = network.compute_simulation( network_dt, network_tf, integration_method, network.neuron_manager, network.synapse_manager, network.applied_current_manager, network.applied_voltage_manager, filter_disabled_flag, set_flag, process_option, undetected_option, network.network_utilities );

% End the timer.
toc


%% Decode the Absolute Inversion Subnetwork Output.

% Decode the network input.
xs = f_decode( Us( 1, : ) );

% Decode the network output.
ys = f_decode( Us( 2, : ) );


%% Plot the Absolute Inversion Subnetwork Results.

% Retrieve the neuron IDs.
neuron_IDs = network.neuron_manager.get_all_neuron_IDs( network.neuron_manager.neurons );

% Plot the network currents over time.
fig_network_currents = network.network_utilities.plot_network_currents( ts, I_leaks, I_syns, I_nas, I_apps, I_totals, neuron_IDs );

% Plot the network states over time.
fig_network_states = network.network_utilities.plot_network_states( ts, Us, hs, neuron_IDs );

% Plot the encoded network input and output over time.
fig_network_encoded = figure( 'Color', 'w', 'Name', 'AI: Encoded Input & Output vs Time' ); hold on, grid on, xlabel( 'Time, t [s]' ), ylabel( 'AI: Encoded Input & Output, U [V]' ), title( 'AI: Encoded Input & Output vs Time' )
plot( ts, Us( 1, : ), '-', 'Linewidth', 3 )
plot( ts, Us( 2, : ), '-', 'Linewidth', 3 )
legend( 'Encoded Input', 'Encoded Output' )
saveas( fig_network_encoded, [ save_directory, '\', 'absolute_inversion_example_encoded' ] )

% Plot the decoded network input and output over time.
fig_network_decoded = figure( 'Color', 'w', 'Name', 'AI: Decoded Input & Output vs Time' ); hold on, grid on, xlabel( 'Time, t [s]' ), ylabel( 'AI: Decoded Input & Output [-]' ), title( 'AI: Decoded Input & Output vs Time' )
plot( ts, xs, '-', 'Linewidth', 3 )
plot( ts, ys, '-', 'Linewidth', 3 )
legend( 'Decoded Input', 'Decoded Output' )
saveas( fig_network_decoded, [ save_directory, '\', 'absolute_inversion_example_decoded' ] )

% Plot the encoded network input and output.
fig_network_encoded = figure( 'Color', 'w', 'Name', 'AI: Decoded Output vs Decoded Input' ); hold on, grid on, xlabel( 'Encoded Input, U1 [V]' ), ylabel( 'Encoded Output, U2 [V]' ), title( 'AI: Encoded Output vs Encoded Input' )
plot( Us( 1, : ), Us( 2, : ), '-', 'Linewidth', 3 )
saveas( fig_network_encoded, [ save_directory, '\', 'absolute_inversion_dynamic_example_encoded' ] )

% Plot the decoded network input and output.
fig_network_decoding = figure( 'Color', 'w', 'Name', 'AI: Decoded Output vs Decoded Input' ); hold on, grid on, xlabel( 'Decoded Input [-]' ), ylabel( 'Decoded Output [-]' ), title( 'AI: Decoded Output vs Decoded Input' )
plot( xs, ys, '-', 'Linewidth', 3 )
saveas( fig_network_decoding, [ save_directory, '\', 'absolute_inversion_dynamic_example_decoded' ] )

% Animate the network states over time.
fig_network_animation = network.network_utilities.animate_network_states( Us, hs, neuron_IDs );

