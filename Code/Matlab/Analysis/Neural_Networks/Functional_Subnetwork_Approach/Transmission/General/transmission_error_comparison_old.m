%% Transmission Subnetwork Encoding Comparison.

% Clear Everything.
clear, close( 'all' ), clc


%% Initialize Project Options.

% Define the save and load directories.
save_directory = '.\Save';                                  % [str] Save Directory.
load_directory = '.\Load';                                  % [str] Load Directory.

% Define the network simulation time step.
network_dt = 1e-3;                                          % [s] Simulation Time Step.

% Define the network simulation duration.
network_tf = 3;                                             % [s] Simulation Duration.


%% Define Absolute Transmission Subnetwork Parameters.

% Define the maximum membrane voltages.
R1_absolute = 20e-3;                                    	% [V] Maximum Membrane Voltage (Neuron 1).

% Define the membrane conductances.
Gm1_absolute = 1e-6;                                       	% [S] Membrane Conductance (Neuron 1)
Gm2_absolute = 1e-6;                                      	% [S] Membrane Conductance (Neuron 2) 

% Define the membrane capacitance.
Cm1_absolute = 5e-9;                                     	% [F] Membrane Capacitance (Neuron 1)
Cm2_absolute = 5e-9;                                      	% [F] Membrane Capacitance (Neuron 2)

% Define the sodium channel conductance.
Gna1_absolute = 0;                                         	% [S] Sodium Channel Conductance (Neuron 1).
Gna2_absolute = 0;                                       	% [S] Sodium Channel Conductance (Neuron 2).

% Define the synaptic conductances.
dEs21_absolute = 194e-3;                                   	% [V] Synaptic Reversal Potential (Synapse 21).

% Define the applied currents.
Ia1_absolute = R1_absolute*Gm1_absolute;                 	% [A] Applied Current (Neuron 1)
Ia2_absolute = 0;                                         	% [A] Applied Current (Neuron 2).

% Define the current state.
current_state1_absolute = 1.0;                             	% [-] Current State (Neuron 1). (Specified as a ratio of the total applied current that is active.)

% Define the network design parameters.
c_absolute = 2;                                             % [-] Design Constant.


%% Compute the Derived Absolute Transmission Subnetwork Parameters.

% Compute the maximum membrane voltages.
R2_absolute = c_absolute*R1_absolute;                                          % [V] Maximum Membrane Voltage (Neuron 2).

% Compute the synaptic conductances.
gs21_absolute = ( R2_absolute*Gm2_absolute - Ia2_absolute )/( dEs21_absolute - R2_absolute );             % [S] Synaptic Conductance (Synapse 21).


%% Print Absolute Transmission Subnetwork Parameters.

% Print out a header.
fprintf( '\n------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )
fprintf( 'ABSOLUTE TRANSMISSION SUBNETWORK PARAMETERS:\n' )
fprintf( '------------------------------------------------------------\n' )

% Print out neuron information.
fprintf( 'Neuron Parameters:\n' )
fprintf( 'R1 \t\t= \t%0.2f \t[mV]\n', R1_absolute*( 10^3 ) )
fprintf( 'R2 \t\t= \t%0.2f \t[mV]\n', R2_absolute*( 10^3 ) )

fprintf( 'Gm1 \t= \t%0.2f \t[muS]\n', Gm1_absolute*( 10^6 ) )
fprintf( 'Gm2 \t= \t%0.2f \t[muS]\n', Gm2_absolute*( 10^6 ) )

fprintf( 'Cm1 \t= \t%0.2f \t[nF]\n', Cm1_absolute*( 10^9 ) )
fprintf( 'Cm2 \t= \t%0.2f \t[nF]\n', Cm2_absolute*( 10^9 ) )

fprintf( 'Gna1 \t= \t%0.2f \t[muS]\n', Gna1_absolute*( 10^6 ) )
fprintf( 'Gna2 \t= \t%0.2f \t[muS]\n', Gna2_absolute*( 10^6 ) )
fprintf( '\n' )

% Print out the synapse information.
fprintf( 'Synapse Parameters:\n' )
fprintf( 'dEs21 \t= \t%0.2f \t[mV]\n', dEs21_absolute*( 10^3 ) )
fprintf( 'gs21 \t= \t%0.2f \t[muS]\n', gs21_absolute*( 10^6 ) )
fprintf( '\n' )

% Print out the applied current information.
fprintf( 'Applied Curent Parameters:\n' )
fprintf( 'Ia1 \t= \t%0.2f \t[nA]\n', current_state1_absolute*Ia1_absolute*( 10^9 ) )
fprintf( '\n' )

% Print out the network design parameters.
fprintf( 'Network Design Parameters:\n' )
fprintf( 'c \t\t= \t%0.2f \t[-]\n', c_absolute )

% Print out ending information.
fprintf( '------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )


%% Create an Absolute Transmission Subnetwork.

% Create an instance of the network class.
network_absolute = network_class( network_dt, network_tf );

% Create the network components.
[ network_absolute.neuron_manager, neuron_IDs_absolute ] = network_absolute.neuron_manager.create_neurons( 2 );
[ network_absolute.synapse_manager, synapse_IDs_absolute ] = network_absolute.synapse_manager.create_synapses( 1 );
[ network_absolute.applied_current_manager, applied_current_IDs_absolute ] = network_absolute.applied_current_manager.create_applied_currents( 2 );

% Set the neuron parameters.
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ Gna1_absolute, Gna2_absolute ], 'Gna' );
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ R1_absolute, R2_absolute ], 'R' );
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ Gm1_absolute, Gm2_absolute ], 'Gm' );
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ Cm1_absolute, Cm2_absolute ], 'Cm' );

% Set the synapse parameters.
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, 1, 'from_neuron_ID' );
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, 2, 'to_neuron_ID' );
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, gs21_absolute, 'g_syn_max' );
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, dEs21_absolute, 'dE_syn' );

% Set the applied current parameters.
network_absolute.applied_current_manager = network_absolute.applied_current_manager.set_applied_current_property( applied_current_IDs_absolute, [ 1, 2 ], 'neuron_ID' );
network_absolute.applied_current_manager = network_absolute.applied_current_manager.set_applied_current_property( applied_current_IDs_absolute, [ current_state1_absolute*Ia1_absolute, Ia2_absolute ], 'I_apps' );


%% Define Basic Relative Transmission Subnetwork Parameters.

% Define the maximum membrane voltages.
R1_relative = 20e-3;                                        % [V] Maximum Membrane Voltage (Neuron 1).
R2_relative = 20e-3;                                        % [V] Maximum Membrane Voltage (Neuron 2).

% Define the membrane conductances.
Gm1_relative = 1e-6;                                       	% [S] Membrane Conductance (Neuron 1)
Gm2_relative = 1e-6;                                      	% [S] Membrane Conductance (Neuron 2) 

% Define the membrane capacitance.
Cm1_relative = 5e-9;                                     	% [F] Membrane Capacitance (Neuron 1)
Cm2_relative = 5e-9;                                      	% [F] Membrane Capacitance (Neuron 2)

% Define the sodium channel conductance.
Gna1_relative = 0;                                          % [S] Sodium Channel Conductance (Neuron 1).
Gna2_relative = 0;                                          % [S] Sodium Channel Conductance (Neuron 2).

% Define the synaptic conductances.
dEs21_relative = 194e-3;                                   	% [V] Synaptic Reversal Potential (Synapse 21).

% Define the applied currents.
Ia1_relative = R1_relative*Gm1_relative;                  	% [A] Applied Current (Neuron 1)
Ia2_relative = 0;                                       	% [A] Applied Current (Neuron 2).

% Define the current state.
current_state1_relative = 1.0;                              % [-] Current State (Neuron 1). (Specified as a ratio of the total applied current that is active.)

% Define the network design parameters.
c_relative = 1;                                          	% [-] Design Constant.


%% Compute Derived Relative Transmission Subnetwork Parameters.

% Compute the synaptic conductances.
gs21_relative = ( R2_relative*Gm2_relative - Ia2_relative )/( dEs21_relative - R2_relative );             % [S] Synaptic Conductance (Synapse 21).


%% Print Relative Transmission Subnetwork Parameters.

% Print out a header.
fprintf( '\n------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )
fprintf( 'RELATIVE TRANSMISSION SUBNETWORK PARAMETERS:\n' )
fprintf( '------------------------------------------------------------\n' )

% Print out neuron information.
fprintf( 'Neuron Parameters:\n' )
fprintf( 'R1 \t\t= \t%0.2f \t[mV]\n', R1_relative*( 10^3 ) )
fprintf( 'R2 \t\t= \t%0.2f \t[mV]\n', R2_relative*( 10^3 ) )

fprintf( 'Gm1 \t= \t%0.2f \t[muS]\n', Gm1_relative*( 10^6 ) )
fprintf( 'Gm2 \t= \t%0.2f \t[muS]\n', Gm2_relative*( 10^6 ) )

fprintf( 'Cm1 \t= \t%0.2f \t[nF]\n', Cm1_relative*( 10^9 ) )
fprintf( 'Cm2 \t= \t%0.2f \t[nF]\n', Cm2_relative*( 10^9 ) )

fprintf( 'Gna1 \t= \t%0.2f \t[muS]\n', Gna1_relative*( 10^6 ) )
fprintf( 'Gna2 \t= \t%0.2f \t[muS]\n', Gna2_relative*( 10^6 ) )
fprintf( '\n' )

% Print out the synapse information.
fprintf( 'Synapse Parameters:\n' )
fprintf( 'dEs21 \t= \t%0.2f \t[mV]\n', dEs21_relative*( 10^3 ) )
fprintf( 'gs21 \t= \t%0.2f \t[muS]\n', gs21_relative*( 10^6 ) )
fprintf( '\n' )

% Print out the applied current information.
fprintf( 'Applied Curent Parameters:\n' )
fprintf( 'Ia1 \t= \t%0.2f \t[nA]\n', current_state1_relative*Ia1_relative*( 10^9 ) )
fprintf( '\n' )

% Print out the network design parameters.
fprintf( 'Network Design Parameters:\n' )
fprintf( 'c \t\t= \t%0.2f \t[-]\n', c_relative )

% Print out ending information.
fprintf( '------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )


%% Create a Relative Transmission Subnetwork.

% Create an instance of the network class.
network_relative = network_class( network_dt, network_tf );

% Create the network components.
[ network_relative.neuron_manager, neuron_IDs ] = network_relative.neuron_manager.create_neurons( 2 );
[ network_relative.synapse_manager, synapse_IDs ] = network_relative.synapse_manager.create_synapses( 1 );
[ network_relative.applied_current_manager, applied_current_IDs ] = network_relative.applied_current_manager.create_applied_currents( 2 );

% Define neuron parameters.
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs, [ R1_relative, R2_relative ], 'R' );
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs, [ Gm1_relative, Gm2_relative ], 'Gm' );
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs, [ Cm1_relative, Cm2_relative ], 'Cm' );
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs, [ Gna1_relative, Gna2_relative ], 'Gna' );

% Define synapse parameters.
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs, 1, 'from_neuron_ID' );
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs, 2, 'to_neuron_ID' );
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs, gs21_relative, 'g_syn_max' );
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs, dEs21_relative, 'dE_syn' );

% Define applied current manager.
network_relative.applied_current_manager = network_relative.applied_current_manager.set_applied_current_property( applied_current_IDs, [ 1, 2 ], 'neuron_ID' );
network_relative.applied_current_manager = network_relative.applied_current_manager.set_applied_current_property( applied_current_IDs, [ current_state1_relative*Ia1_relative, Ia2_relative ], 'I_apps' );


%% Load the Absolute & Relative Transmission Subnetworks.

% Load the simulation results.
absolute_transmission_simulation_data = load( [ load_directory, '\', 'absolute_transmission_subnetwork_error' ] );
relative_transmission_simulation_data = load( [ load_directory, '\', 'relative_transmission_subnetwork_error' ] );

% Store the absolute simulation results in separate variables.
absolute_applied_currents = absolute_transmission_simulation_data.applied_currents;
Us_achieved_absolute = absolute_transmission_simulation_data.Us_achieved;

% Store the relative simulation results in separate variables.
relative_applied_currents = relative_transmission_simulation_data.applied_currents;
Us_achieved_relative = relative_transmission_simulation_data.Us_achieved;


%% Compute the Error in the Steady State Transmission Subnetwork Responses.

% Store the absolute and relative maximum membrane voltages into arrays.
Rs_absolute = [ R1_absolute, R2_absolute ];                                     % [V] Maximum Membrane Voltages (Absolute Neurons 1 & 2).
Rs_relative = [ R1_relative, R2_relative ];                                     % [V] Maximum Membrane Voltages (Relative Neurons 1 & 2).

% Compute the decoded absolute and relative maximum membrane voltages.
Rs_decoded_absolute = Rs_absolute*( 10^3 );
Rs_decoded_relative = Rs_absolute*( 10^3 );

% Compute the desired steady state output membrane voltage.
Us_desired_absolute_output = c_absolute*Us_achieved_absolute( :, 1 );
Us_desired_relative_output = c_relative*( R2_relative/R1_relative )*Us_achieved_relative( :, 1 );

% Generate desired steady state membrane voltage matrices.
Us_desired_absolute = Us_achieved_absolute; Us_desired_absolute( :, end ) = Us_desired_absolute_output;
Us_desired_relative = Us_achieved_relative; Us_desired_relative( :, end ) = Us_desired_relative_output;

% Compute the absolute desired and achieved decoded steady state values.
Us_achieved_decoded_absolute = Us_achieved_absolute*( 10^3 );
Us_desired_decoded_absolute = Us_desired_absolute*( 10^3 );

% Compute the relative desired and achieved encoded and decoded steady state values.
Us_achieved_encoded_relative = Us_achieved_relative/Rs_relative;
Us_desired_encoded_relative = Us_desired_relative/Rs_relative;
Us_achieved_decoded_relative = ( Rs_absolute./Rs_relative ).*Us_achieved_relative*( 10^3 );
Us_desired_decoded_relative = ( Rs_absolute./Rs_relative ).*Us_desired_relative*( 10^3 );

% Compute the error between the achieved and desired membrane voltage results.
error_absolute = Us_achieved_absolute( :, end ) - Us_desired_absolute( :, end );
error_relative = Us_achieved_relative( :, end ) - Us_desired_relative( :, end );

% Compute the error between the achieved and desired decoded results.
error_decoded_absolute = Us_achieved_decoded_absolute( :, end ) - Us_desired_decoded_absolute( :, end );
error_decoded_relative = Us_achieved_decoded_relative( :, end ) - Us_desired_decoded_relative( :, end );

% Compute the percent error between the achieve and desired results.
error_absolute_percent = 100*( error_absolute/R2_absolute );
error_relative_percent = 100*( error_relative/R2_relative );
error_decoded_absolute_percent = 100*( error_decoded_absolute/Rs_decoded_absolute( 2 ) );
error_decoded_relative_percent = 100*( error_decoded_relative/Rs_decoded_relative( 2 ) );

% Compute the mean squared error.
mse_absolute = ( 1/numel( error_absolute ) )*sqrt( sum( error_absolute.^2, 'all' ) );
mse_relative = ( 1/numel( error_relative ) )*sqrt( sum( error_relative.^2, 'all' ) );
mse_decoded_absolute = ( 1/numel( error_decoded_absolute ) )*sqrt( sum( error_decoded_absolute.^2, 'all' ) );
mse_decoded_relative = ( 1/numel( error_decoded_relative ) )*sqrt( sum( error_decoded_relative.^2, 'all' ) );

% Compute the mean squared error percentage.
mse_absolute_percent = 100*( mse_absolute/R2_absolute );
mse_relative_percent = 100*( mse_relative/R2_relative );
mse_decoded_absolute_percent = 100*( mse_decoded_absolute/Rs_decoded_absolute( 2 ) );
mse_decoded_relative_percent = 100*( mse_decoded_relative/Rs_decoded_relative( 2 ) );

% Compute the standard deviation of the error.
std_absolute = std( error_absolute, 0, 'all' );
std_relative = std( error_relative, 0, 'all' );
std_decoded_absolute = std( error_decoded_absolute, 0, 'all' );
std_decoded_relative = std( error_decoded_relative, 0, 'all' );

% Compute the standard deviation of the error percentage.
std_absolute_percent = 100*( std_absolute/R2_absolute );
std_relative_percent = 100*( std_relative/R2_relative );
std_decoded_absolute_percent = 100*( std_decoded_absolute/Rs_decoded_absolute( 2 ) );
std_decoded_relative_percent = 100*( std_decoded_relative/Rs_decoded_relative( 2 ) );

% Compute the maximum errors.
[ error_absolute_max, index_absolute_max ] = max( abs( error_absolute ), [  ], 'all', 'linear' );
[ error_relative_max, index_relative_max ] = max( abs( error_relative ), [  ], 'all', 'linear' );
[ error_decoded_absolute_max, index_decoded_absolute_max ] = max( abs( error_decoded_absolute ), [  ], 'all', 'linear' );
[ error_decoded_relative_max, index_decoded_relative_max ] = max( abs( error_decoded_relative ), [  ], 'all', 'linear' );

% Compute the maximum error percentages.
error_absolute_max_percent = 100*( error_absolute_max/R2_absolute );
error_relative_max_percent = 100*( error_relative_max/R2_relative );
error_decoded_absolute_max_percent = 100*( error_decoded_absolute_max/Rs_decoded_absolute( 2 ) );
error_decoded_relative_max_percent = 100*( error_decoded_relative_max/Rs_decoded_relative( 2 ) );

% Compute the minimum errors.
[ error_absolute_min, index_absolute_min ] = min( abs( error_absolute ), [  ], 'all', 'linear' );
[ error_relative_min, index_relative_min ] = min( abs( error_relative ), [  ], 'all', 'linear' );
[ error_decoded_absolute_min, index_decoded_absolute_min ] = min( abs( error_decoded_absolute ), [  ], 'all', 'linear' );
[ error_decoded_relative_min, index_decoded_relative_min ] = min( abs( error_decoded_relative ), [  ], 'all', 'linear' );

% Compute the minimum error percentages.
error_absolute_min_percent = 100*( error_absolute_min/R2_absolute );
error_relative_min_percent = 100*( error_relative_min/R2_relative );
error_decoded_absolute_min_percent = 100*( error_decoded_absolute_min/Rs_decoded_absolute( 2 ) );
error_decoded_relative_min_percent = 100*( error_decoded_relative_min/Rs_decoded_relative( 2 ) );

% Compute the range of the error.
error_absolute_range = error_absolute_max - error_absolute_min;
error_relative_range = error_relative_max - error_relative_min;
error_decoded_absolute_range = error_decoded_absolute_max - error_decoded_absolute_min;
error_decoded_relative_range = error_decoded_relative_max - error_decoded_relative_min;

% Compute the range of the error percentages.
error_absolute_range_percent = 100*( error_absolute_range/R2_absolute );
error_relative_range_percent = 100*( error_relative_range/R2_relative );
error_decoded_absolute_range_percent = 100*( error_decoded_absolute_range/Rs_decoded_absolute( 2 ) );
error_decoded_relative_range_percent = 100*( error_decoded_relative_range/Rs_decoded_relative( 2 ) );

% Compute the difference in error between the absolute and relative schemes.
error_difference = abs( error_relative ) - abs( error_absolute );
error_difference_percent = abs( error_relative_percent ) - abs( error_absolute_percent );
error_decoded_difference = abs( error_decoded_relative ) - abs( error_decoded_absolute );
error_decoded_difference_percent = abs( error_decoded_relative_percent ) - abs( error_decoded_absolute_percent );

% Compute the mean squared error difference.
error_difference_mse = abs( mse_relative ) - abs( mse_absolute );
error_difference_mse_percent = abs( mse_relative_percent ) - abs( mse_absolute_percent );
error_decoded_difference_mse = abs( mse_decoded_relative ) - abs( mse_decoded_absolute );
error_decoded_difference_mse_percent = abs( mse_decoded_relative_percent ) - abs( mse_decoded_absolute_percent );

% Compute the standard deviation difference.
error_difference_std = abs( std_relative ) - abs( std_absolute );
error_difference_std_percent = abs( std_relative_percent ) - abs( std_absolute_percent );
error_decoded_difference_std = abs( std_decoded_relative ) - abs( std_decoded_absolute );
error_decoded_difference_std_percent = abs( std_decoded_relative_percent ) - abs( std_decoded_absolute_percent );

% Compute the maximum error difference.
error_difference_max = abs( error_relative_max ) - abs( error_absolute_max );
error_difference_max_percent = abs( error_relative_max_percent ) - abs( error_absolute_max_percent );
error_decoded_difference_max = abs( error_decoded_relative_max ) - abs( error_decoded_absolute_max );
error_decoded_difference_max_percent = abs( error_decoded_relative_max_percent ) - abs( error_decoded_absolute_max_percent );


%% Print Out the Summary Information.

% Retrieve the absolute input voltage matrices.
Us1_achieved_absolute = Us_achieved_absolute( :, 1 );
Us2_achieved_absolute = Us_achieved_absolute( :, 2 );

% Retrieve the relative input voltage matrices.
Us1_achieved_relative = Us_achieved_relative( :, 1 );
Us2_achieved_relative = Us_achieved_relative( :, 2 );

% Print out the absolute transmission membrane voltage summary statistics.
fprintf( 'Absolute Transmission Summary Statistics (Membrane Voltages)\n' )
fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_absolute, mse_absolute_percent )
fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_absolute, std_absolute_percent )
fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_absolute_max, error_absolute_max_percent, Us1_achieved_absolute( index_absolute_max ), Us2_achieved_absolute( index_absolute_max ) )
fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_absolute_min, error_absolute_min_percent, Us1_achieved_absolute( index_absolute_min ), Us2_achieved_absolute( index_absolute_min ) )
fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_absolute_range, error_absolute_range_percent )

fprintf( '\n' )
fprintf( 'Relative Transmission Summary Statistics (Membrane Voltages)\n' )
fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_relative, mse_relative_percent )
fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_relative, std_relative_percent )
fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_relative_max, error_relative_max_percent, Us1_achieved_relative( index_relative_max ), Us2_achieved_relative( index_relative_max ) )
fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_relative_min, error_relative_min_percent, Us1_achieved_relative( index_relative_min ), Us2_achieved_relative( index_relative_min ) )
fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_relative_range, error_relative_range_percent )

fprintf( '\n' )
fprintf( 'Absolute vs Relative Transmission Summary Statistics (Membrane Voltages):\n' )
fprintf( 'delta MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', error_difference_mse, error_difference_mse_percent )
fprintf( 'delta STD:\t%9.3e [V] (%6.2f [%%])\n', error_difference_std, error_difference_std_percent )
fprintf( 'delta Max Error:\t%9.3e [mV] (%6.2f [%%])\n', error_difference_max, error_difference_max_percent )

% Print out the absolute transmission decoding summary statistics.
fprintf( 'Absolute Transmission Summary Statistics (Decoded)\n' )
fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_decoded_absolute, mse_decoded_absolute_percent )
fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_decoded_absolute, std_decoded_absolute_percent )
fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_decoded_absolute_max, error_decoded_absolute_max_percent, Us_achieved_decoded_absolute( index_decoded_absolute_max, 1 ), Us_achieved_decoded_absolute( index_decoded_absolute_max, 2 ) )
fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_decoded_absolute_min, error_decoded_absolute_min_percent, Us_achieved_decoded_absolute( index_decoded_absolute_min, 1 ), Us_achieved_decoded_absolute( index_decoded_absolute_min, 2 ) )
fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_decoded_absolute_range, error_decoded_absolute_range_percent )

fprintf( '\n' )
fprintf( 'Relative Transmission Summary Statistics (Membrane Voltages)\n' )
fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_decoded_relative, mse_decoded_relative_percent )
fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_decoded_relative, std_decoded_relative_percent )
fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_decoded_relative_max, error_decoded_relative_max_percent, Us_achieved_decoded_relative( index_decoded_relative_max, 1 ), Us_achieved_decoded_relative( index_decoded_relative_max, 2 ) )
fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV])\n', error_decoded_relative_min, error_decoded_relative_min_percent, Us_achieved_decoded_relative( index_decoded_relative_min, 1 ), Us_achieved_decoded_relative( index_decoded_relative_min, 2 ) )
fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_decoded_relative_range, error_decoded_relative_range_percent )

fprintf( '\n' )
fprintf( 'Absolute vs Relative Transmission Summary Statistics (Membrane Voltages):\n' )
fprintf( 'delta MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', error_decoded_difference_mse, error_decoded_difference_mse_percent )
fprintf( 'delta STD:\t%9.3e [V] (%6.2f [%%])\n', error_decoded_difference_std, error_decoded_difference_std_percent )
fprintf( 'delta Max Error:\t%9.3e [mV] (%6.2f [%%])\n', error_decoded_difference_max, error_decoded_difference_max_percent )


%% Plot the Steady State Transmission Error Surfaces.

% Create a figure that shows the achieved and desired membrane voltage outputs for the absolute transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Absolute Transmission Steady State Response (Comparison)' ); hold on, grid on, xlabel( 'Input Membrane Voltage, U1 [mV]' ), ylabel( 'Output Membrane Voltage, U2 [mV]' ), title( 'Absolute Transmission Steady State Response (Comparison)' )
plot( Us_desired_absolute( :, 1 )*( 10^3 ), Us_desired_absolute( :, end )*( 10^3 ), '-', 'Linewidth', 3 )
plot( Us_achieved_absolute( :, 1 )*( 10^3 ), Us_achieved_absolute( :, end )*( 10^3 ), '-', 'Linewidth', 3 )
legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'absolute_transmission_ss_response_comparison.png' ] )

% Create a figure that shows the decoded achieved and desired membrane voltage outputs for the absolute transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Absolute Transmission Steady State Decoding (Comparison)' ); hold on, grid on, xlabel( 'Input Decoding [-]' ), ylabel( 'Output Decoding [-]' ), title( 'Absolute Transmission Steady State Decoding (Comparison)' )
plot( Us_desired_decoded_absolute( :, 1 ), Us_desired_decoded_absolute( :, end ), '-', 'Linewidth', 3 )
plot( Us_achieved_decoded_absolute( :, 1 ), Us_achieved_decoded_absolute( :, end ), '-', 'Linewidth', 3 )
legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'absolute_transmission_ss_decoding_comparison.png' ] )

% Create a figure that shows the differences between the achieved and desired membrane voltage outputs for the relative transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Response (Comparison)' ); hold on, grid on, xlabel( 'Membrane Voltage of Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Output Neuron, U2 [mV]' ), title( 'Relative Transmission Steady State Response (Comparison)' )
plot( Us_desired_relative( :, 1 )*( 10^3 ), Us_desired_relative( :, end )*( 10^3 ), '-', 'Linewidth', 3 )
plot( Us_achieved_relative( :, 1 )*( 10^3 ), Us_achieved_relative( :, end )*( 10^3 ), '-', 'Linewidth', 3 )
legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_comparison.png' ] )

% Create a figure that shows the achieved and desired encoding outputs for the relative transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Encoding (Comparison)' ); hold on, grid on, xlabel( 'Input Encoding [-]' ), ylabel( 'Output Encoding [-]' ), title( 'Relative Transmission Steady State Encoding (Comparison)' )
plot( Us_desired_encoded_relative( :, 1 ), Us_desired_encoded_relative( :, end ), '-', 'Linewidth', 3 )
plot( Us_achieved_encoded_relative( :, 1 ), Us_achieved_encoded_relative( :, end ), '-', 'Linewidth', 3 )
legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_encoding_comparison.png' ] )

% Create a figure that shows the achieved and desired decoding outputs for the relative transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Decoding (Comparison)' ); hold on, grid on, xlabel( 'Input Decoding [-]' ), ylabel( 'Output Decoding [-]' ), title( 'Relative Transmission Steady State Decoding (Comparison)' )
plot( Us_desired_decoded_relative( :, 1 ), Us_desired_decoded_relative( :, end ), '-', 'Linewidth', 3 )
plot( Us_achieved_decoded_relative( :, 1 ), Us_achieved_decoded_relative( :, end ), '-', 'Linewidth', 3 )
legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_decoding_comparison.png' ] )

% Create a figure that shows the achieved and desired membrane voltage outputs for the relative transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Response (Comparison)' ); hold on, grid on, xlabel( 'Membrane Voltage of Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Output Neuron, U2 [mV]' ), title( 'Transmission Steady State Response (Comparison)' )
plot( Us_desired_absolute( :, 1 )*( 10^3 ), Us_desired_absolute( :, end )*( 10^3 ), 'r-', 'Linewidth', 3 )
plot( Us_achieved_absolute( :, 1 )*( 10^3 ), Us_achieved_absolute( :, end )*( 10^3 ), 'r--', 'Linewidth', 3 )
plot( Us_desired_relative( :, 1 )*( 10^3 ), Us_desired_relative( :, end )*( 10^3 ), 'b-', 'Linewidth', 3 )
plot( Us_achieved_relative( :, 1 )*( 10^3 ), Us_achieved_relative( :, end )*( 10^3 ), 'b--', 'Linewidth', 3 )
legend( { 'Absolute Desired', 'Absolute Achieved', 'Relative Desired', 'Relative Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_ss_response_comparison.png' ] )

% Create a figure that shows the achieved and desired decoding outputs for the relative transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Decoding (Comparison)' ); hold on, grid on, xlabel( 'Input Decoding [-]' ), ylabel( 'Output Decoding [-]' ), title( 'Transmission Steady State Decoding (Comparison)' )
plot( Us_desired_decoded_absolute( :, 1 ), Us_desired_decoded_absolute( :, end ), 'r-', 'Linewidth', 3 )
plot( Us_achieved_decoded_absolute( :, 1 ), Us_achieved_decoded_absolute( :, end ), 'r--', 'Linewidth', 3 )
plot( Us_desired_decoded_relative( :, 1 ), Us_desired_decoded_relative( :, end ), 'b-', 'Linewidth', 3 )
plot( Us_achieved_decoded_relative( :, 1 ), Us_achieved_decoded_relative( :, end ), 'b--', 'Linewidth', 3 )
legend( { 'Absolute Desired', 'Absolute Achieved', 'Relative Desired', 'Relative Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_ss_response_comparison.png' ] )

% Create a surface that shows the membrane voltage error for the transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Error' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Membrane Voltage Error, E [mV]' ), title( 'Transmission Steady State Error' )
plot( Us_achieved_absolute( :, 1 )*( 10^3 ), error_absolute*( 10^3 ), '-', 'Linewidth', 3 )
plot( Us_achieved_relative( :, 1 )*( 10^3 ), error_relative*( 10^3 ), '-', 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_error_comparison.png' ] )

% Create a surface that shows the decoding error for the transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Decoding Error' ); hold on, grid on, xlabel( 'Input Decoding [-]' ), ylabel( 'Output Decoding Error [-]' ), title( 'Transmission Steady State Decoding Error' )
plot( Us_achieved_decoded_absolute( :, 1 ), error_decoded_absolute, '-', 'Linewidth', 3 )
plot( Us_achieved_decoded_relative( :, 1 ), error_decoded_relative, '-', 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_decoding_error_comparison.png' ] )

% Create a surface that shows the membrane voltage error percentage of the transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Error Percentage' ); hold on, grid on, xlabel( 'Membrane Voltage of Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage Error Percentage, E [%]' ), title( 'Transmission Steady State Error Percentage' )
plot( Us_achieved_absolute( :, 1 )*( 10^3 ), error_absolute_percent, '-', 'Linewidth', 3 )
plot( Us_achieved_relative( :, 1 )*( 10^3 ), error_relative_percent, '-', 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_error_percentage_comparison.png' ] )

% Create a surface that shows the decoding error percentage for the transmission subnetwork.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Decoding Error Percentage' ); hold on, grid on, xlabel( 'Input Decoding [-]' ), ylabel( 'Output Decoding Error Percentage [%]' ), title( 'Transmission Steady State Decoding Error Percentage' )
plot( Us_achieved_decoded_absolute( :, 1 ), error_decoded_absolute, '-', 'Linewidth', 3 )
plot( Us_achieved_decoded_relative( :, 1 ), error_decoded_relative, '-', 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_decoding_error_percentage_comparison.png' ] )

% Create a surface that shows the difference in error between the absolute and relative transmission subnetworks.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Error Difference' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Membrane Voltage Error Difference, dE [mV]' ), title( 'Transmission Steady State Error Difference' )
plot( Us_achieved_absolute( :, 1 )*( 10^3 ), error_difference*( 10^3 ), '-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'transmission_error_difference.png' ] )

% Create a surface that shows the difference in decoding error between the absolute and relative transmission subnetworks.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Decoding Error Difference' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Decoding Error Difference, dE [mV]' ), title( 'Transmission Steady State Decoding Error Difference' )
plot( Us_achieved_decoded_absolute( :, 1 ), error_decoded_difference, '-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'transmission_decoding_error_difference.png' ] )

% Create a surface that shows the difference in error between the absolute and relative percent transmission subnetworks.
fig = figure( 'Color', 'w', 'Name', 'Transmission Steady State Decoding Error Difference Percentage' ); hold on, grid on, xlabel( 'Membrane Voltage of Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage Error Difference Percentage, dE [%]' ), title( 'Transmission Steady State Decoding Error Difference Percentage' )
plot( Us_achieved_decoded_absolute( :, 1 ), error_decoded_difference_percent, 'b-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'transmission_decoding_error_percentage_difference.png' ] )

