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
% network_dt = 1.3e-4;                                % [s] Simulation Timestep.
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
c3 = 20e-9;                                         % [A] Design Constant 2.
delta = 1e-3;                                       % [V] Membrane Voltage Offset.

% Define the desired mapping operation.
f_desired = @( x, c2 ) network_utilities.compute_desired_inversion_sso( x, c1, c2, c3 );


%% Define the Encoding & Decoding Operations.

% Define the domain of the input and output signals.
x_max = 20;
% y_max = f_desired( x_max );

% Define the encoding scheme.
f_encode = @( x ) x*( 10^( -3 ) );

% Define the decoding scheme.
f_decode = @( U ) U*( 10^3 );


%% Define Additional Absolute Inversion Design Subnetwork Parameters.

% Define the inversion subnetwork design parameters.
R1 = 20e-3;                                         % [V] Maximum Membrane Voltage (Neuron 1).
Gm1 = 1e-6;                                         % [S] Membrane Conductance (Neuron 1).
Gm2 = 1e-6;                                       	% [S] Membrane Conductance (Neuron 2).
Cm1 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 1).
Cm2 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 2).
% Cm1 = 30e-9;                                     	% [F] Membrane Capacitance (Neuron 1).
% Cm2 = 30e-9;                                      % [F] Membrane Capacitance (Neuron 2).

% Store the inversion subnetwork design parameters in a cell.
inversion_parameters = { c1, c3, delta, R1, Gm1, Gm2, Cm1, Cm2 };


%% Define the Desired Input Signal.

% Define the desired decoded input signal.
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
[ c, Gnas, R2, dEs21, gs21, Ia2, neurons, synapses, neuron_manager, synapse_manager, network ] = network.create_inversion_subnetwork( inversion_parameters, encoding_scheme, network.neuron_manager, network.synapse_manager, network.applied_current_manager, true, true, false, undetected_option );

% Create the input applied current.
[ ~, ~, ~, network.applied_current_manager ] = network.applied_current_manager.create_applied_current( input_current_ID, input_current_name, input_current_to_neuron_ID, ts, Ias1, true, network.applied_current_manager.applied_currents, true, false, network.applied_current_manager.array_utilities );




%% Compute Derived Absolute Inversion Subnetwork Constraints.

% Compute the network design parameters.
c2 = ( c1 - delta*c3 )/( delta*R1 );                % [S] Design Constant 2.

% Compute the maximum membrane voltages.
R2 = c1/c3;                                         % [V] Maximum Membrane Voltage (Neuron 2).

% Compute the applied currents.
Ia2 = R2*Gm2;                                       % [A] Applied Current (Neuron 2).

% Compute the synaptic conductances.
gs21 = ( delta*Gm2 - Ia2 )/( dEs21 - delta );       % [S] Synaptic Conductance (Synapse 21).


%% Print Absolute Inversion Subnetwork Parameters.

% Print out a header.
fprintf( '\n------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )
fprintf( 'ABSOLUTE INVERSION SUBNETWORK PARAMETERS:\n' )
fprintf( '------------------------------------------------------------\n' )

% Print out neuron information.
fprintf( 'Neuron Parameters:\n' )
fprintf( 'R1 \t\t= \t%0.2f \t[mV]\n', R1*( 10^3 ) )
fprintf( 'R2 \t\t= \t%0.2f \t[mV]\n', R2*( 10^3 ) )

fprintf( 'Gm1 \t= \t%0.2f \t[muS]\n', Gm1*( 10^6 ) )
fprintf( 'Gm2 \t= \t%0.2f \t[muS]\n', Gm2*( 10^6 ) )

fprintf( 'Cm1 \t= \t%0.2f \t[nF]\n', Cm1*( 10^9 ) )
fprintf( 'Cm2 \t= \t%0.2f \t[nF]\n', Cm2*( 10^9 ) )

fprintf( 'Gna1 \t= \t%0.2f \t[muS]\n', Gna1*( 10^6 ) )
fprintf( 'Gna2 \t= \t%0.2f \t[muS]\n', Gna2*( 10^6 ) )
fprintf( '\n' )

% Print out the synapse information.
fprintf( 'Synapse Parameters:\n' )
fprintf( 'dEs21 \t= \t%0.2f \t[mV]\n', dEs21*( 10^3 ) )
fprintf( 'gs21 \t= \t%0.2f \t[muS]\n', gs21*( 10^6 ) )
fprintf( '\n' )

% Print out the applied current information.
fprintf( 'Applied Curent Parameters:\n' )
fprintf( 'Ia1 \t= \t%0.2f \t[nA]\n', current_state1*Ia1*( 10^9 ) )
fprintf( 'Ia2 \t= \t%0.2f \t[nA]\n', Ia2*( 10^9 ) )
fprintf( '\n' )

% Print out the network design parameters.
fprintf( 'Network Design Parameters:\n' )
fprintf( 'c1 \t\t= \t%0.2f \t[nW]\n', c1*( 10^9 ) )
fprintf( 'c2 \t\t= \t%0.2f \t[muS]\n', c2*( 10^6 ) )
fprintf( 'c3 \t\t= \t%0.2f \t[nA]\n', c3*( 10^9 ) )
fprintf( 'delta \t= \t%0.2f \t[mV]\n', delta*( 10^3 ) )

% Print out ending information.
fprintf( '------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )


%% Create Absolute Inversion Subnetwork.

% Create an instance of the network class.
network = network_class( network_dt, network_tf );

% Create the network components.
[ network.neuron_manager, neuron_IDs ] = network.neuron_manager.create_neurons( 2 );
[ network.synapse_manager, synapse_IDs ] = network.synapse_manager.create_synapses( 1 );
[ network.applied_current_manager, applied_current_IDs ] = network.applied_current_manager.create_applied_currents( 2 );

% Set the neuron parameters.
network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs, [ Gna1, Gna2 ], 'Gna' );
network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs, [ R1, R2 ], 'R' );
network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs, [ Gm1, Gm2 ], 'Gm' );
network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs, [ Cm1, Cm2 ], 'Cm' );

% Set the synapse parameters.
network.synapse_manager = network.synapse_manager.set_synapse_property( synapse_IDs, 1, 'from_neuron_ID' );
network.synapse_manager = network.synapse_manager.set_synapse_property( synapse_IDs, 2, 'to_neuron_ID' );
network.synapse_manager = network.synapse_manager.set_synapse_property( synapse_IDs, gs21, 'g_syn_max' );
network.synapse_manager = network.synapse_manager.set_synapse_property( synapse_IDs, dEs21, 'dE_syn' );

% Set the applied current parameters.
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs, [ 1, 2 ], 'neuron_ID' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs, [ Ia1, Ia2 ], 'I_apps' );


%% Compute Absolute Inversion Numerical Stability Analysis Parameters.

% Compute the maximum RK4 step size and condition number.
[ A, dt_max, condition_number ] = network.RK4_stability_analysis( cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Cm' ) ), cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Gm' ) ), cell2mat( network.neuron_manager.get_neuron_property( 'all', 'R' ) ), network.get_gsynmaxs( 'all' ), network.get_dEsyns( 'all' ), zeros( network.neuron_manager.num_neurons, 1 ), 1e-6 );
% [ A, dt_max, condition_number ] = network.RK4_stability_analysis( cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Cm' ) ), cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Gm' ) ), cell2mat( network.neuron_manager.get_neuron_property( 'all', 'R' ) ), network.get_gsynmaxs( 'all' ), network.get_dEsyns( 'all' ), [ 0; Ia2/Gm2 ], 1e-6 );
% [ A, dt_max, condition_number ] = network.RK4_stability_analysis( cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Cm' ) ), cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Gm' ) ), cell2mat( network.neuron_manager.get_neuron_property( 'all', 'R' ) ), network.get_gsynmaxs( 'all' ), network.get_dEsyns( 'all' ), [ ( ( delta*Gm2 - Ia2 )*R1 )/( ( dEs21 - delta )*gs21 ); delta ], 1e-6 );

% Print out the stability information.
fprintf( '\nSTABILITY SUMMARY:\n' )
fprintf( 'Linearized System Matrix: A =\n\n' ), disp( A )
fprintf( 'Max RK4 Step Size: \tdt_max = %0.3e [s]\n', dt_max )
fprintf( 'Proposed Step Size: \tdt = %0.3e [s]\n', network_dt )
fprintf( 'Condition Number: \tcond( A ) = %0.3e [-]\n', condition_number )


%% Simulate the Absolute Inversion Subnetwork.

% Start the timer.
tic

% Simulate the network.
[ network, ts, Us, hs, dUs, dhs, G_syns, I_leaks, I_syns, I_nas, I_apps, I_totals, m_infs, h_infs, tauhs, neuron_IDs ] = network.compute_set_simulation(  );

% End the timer.
toc


%% Plot the Absolute Inversion Subnetwork Results.

% Compute the decoded output.
Us_decoded = Us*( 10^3 );

% Plot the network currents over time.
fig_network_currents = network.network_utilities.plot_network_currents( ts, I_leaks, I_syns, I_nas, I_apps, I_totals, neuron_IDs );

% Plot the network states over time.
fig_network_states = network.network_utilities.plot_network_states( ts, Us, hs, neuron_IDs );

% Plot the absolute network decoding over time.
fig_network_decoding = figure( 'Color', 'w', 'Name', 'Absolute Inversion Decoding vs Time' ); hold on, grid on, xlabel( 'Time, t [s]' ), ylabel( 'Network Decoding [-]' ), title( 'Absolute Inversion Decoding vs Time' )
plot( ts, Us_decoded( 1, : ), '-', 'Linewidth', 3 )
plot( ts, Us_decoded( 2, : ), '-', 'Linewidth', 3 )
legend( 'Input', 'Output' )
saveas( fig_network_decoding, [ save_directory, '\', 'absolute_inversion_decoding_example' ] )

% Plot the absolute network dynamic decoding example.
fig_network_decoding = figure( 'Color', 'w', 'Name', 'Absolute Inversion Dynamic Decoding Example' ); hold on, grid on, xlabel( 'Input, x [-]' ), ylabel( 'Output, y [-]' ), title( 'Absolute Inversion Dynamic Decoding Example' )
plot( Us_decoded( 1, : ), Us_decoded( 2, : ), '-', 'Linewidth', 3 )
saveas( fig_network_decoding, [ save_directory, '\', 'absolute_inversion_dynamic_decoding_example' ] )

% Animate the network states over time.
fig_network_animation = network.network_utilities.animate_network_states( Us, hs, neuron_IDs );

