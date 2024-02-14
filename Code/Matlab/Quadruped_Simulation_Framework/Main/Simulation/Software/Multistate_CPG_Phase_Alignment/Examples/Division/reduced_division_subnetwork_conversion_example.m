%% Reduced Division Subnetwork Conversion Example.

% Clear everything.
clear, close( 'all' ), clc


%% Set the Simulation Parameters.

% Define the network integration step size.
network_dt = 1.3e-4;                                                                                    % [s] Simulation Step Size.

% Define the network simulation duration.
network_tf = 3;                                                                                         % [s] Simulation Duration.

% Define the applied current state.
% current_state1 = 0;                                                                                 	% [0-1] Current Activation Ratio (Neuron 1). (0 = Input Current Completely Off, 1 = Input Current Completely On.)
% current_state1 = 0.25;                                                                              	% [0-1] Current Activation Ratio (Neuron 1). (0 = Input Current Completely Off, 1 = Input Current Completely On.)
current_state1 = 1;                                                                                   	% [0-1] Current Activation Ratio (Neuron 1). (0 = Input Current Completely Off, 1 = Input Current Completely On.)

% current_state2 = 0;                                                                                   % [0-1] Current Activation Ratio (Neuron 2). (0 = Input Current Completely Off, 1 = Input Current Completely On.)
% current_state2 = 0.25;                                                                                % [0-1] Current Activation Ratio (Neuron 2). (0 = Input Current Completely Off, 1 = Input Current Completely On.)
current_state2 = 1;                                                                                     % [0-1] Current Activation Ratio (Neuron 2).(0 = Input Current Completely Off, 1 = Input Current Completely On.)


%% Define the Fundamental Parameters of a Reduced Relative Division Subnetwork.

% Define the maximum membrane voltages.
R1_relative = 20e-3;                                                                                    % [V] Maximum Membrane Voltage (Neuron 1).
R2_relative = 20e-3;                                                                                    % [V] Maximum Membrane Voltage (Neuron 2).
R3_relative = 20e-3;                                                                                    % [V] Maximum Membrane Voltage (Neuron 3).

% Define the membrane conductances.
Gm1_relative = 1e-6;                                                                                    % [S] Membrane Conductance (Neuron 1).
Gm2_relative = 1e-6;                                                                                    % [S] Membrane Conductance (Neuron 2).
Gm3_relative = 1e-6;                                                                                    % [S] Membrane Conductance (Neuron 3).

% Define the membrane capacitances.
Cm1_relative = 5e-9;                                                                                    % [F] Membrane Conductance (Neuron 1).
Cm2_relative = 5e-9;                                                                                    % [F] Membrane Conductance (Neuron 2).
Cm3_relative = 5e-9;                                                                                    % [F] Membrane Conductance (Neuron 3).

% Define the sodium channel conductances.
Gna1_relative = 0;                                                                                    	% [S] Sodium Channel Conductance (Neuron 1).
Gna2_relative = 0;                                                                                   	% [S] Sodium Channel Conductance (Neuron 2).
Gna3_relative = 0;                                                                                   	% [S] Sodium Channel Conductance (Neuron 2).

% Define the synaptic conductances.
dEs31_relative = 194e-3;                                                                                % [V] Synaptic Reversal Potential (Synapse 31).
dEs32_relative = 0;                                                                                     % [V] Synaptic Reversal Potential (Synapse 32).

% Define the applied currents.
Ia1_relative = R1_relative*Gm1_relative;                                                            	% [A] Applied Current (Neuron 1).
Ia2_relative = R2_relative*Gm2_relative;                                                            	% [A] Applied Current (Neuron 2).
Ia3_relative = 0;                                                                                       % [A] Applied Current (Neuron 3).

% Define network design parameters.
delta_relative = 1e-3;                                                                                  % [V] Membrane Voltage Offset.


%% Compute the Derived Parameters of a Reduced Relative Division Subnetwork.

% Compute the network design parameters.
c1_relative = delta_relative/( R3_relative - delta_relative );                                                                                               	% [-] Relative Division Parameter 1.
c2_relative = c1_relative;                                                                                                                                      % [-] Relative Division Parameter 2.

% Compute the maximum synaptic conductances.
gs31_relative = ( R3_relative*Gm3_relative - Ia3_relative )/( dEs31_relative - R3_relative );                                                                   % [S] Maximum Synaptic Conductance (Synapse 31).
gs32_relative = ( ( dEs31_relative - delta_relative )*gs31_relative + Ia3_relative - delta_relative*Gm3_relative )/( delta_relative - dEs32_relative );         % [S] Maximum Synaptic Conductance (Synapse 32).


%% Print Relative Division Subnetwork Parameters.

% Print out a header.
fprintf( '\n------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )
fprintf( 'RELATIVE DIVISION SUBNETWORK PARAMETERS:\n' )
fprintf( '------------------------------------------------------------\n' )

% Print out neuron information.
fprintf( 'Neuron Parameters:\n' )
fprintf( 'R1 \t\t= \t%0.2f \t[mV]\n', R1_relative*( 10^3 ) )
fprintf( 'R2 \t\t= \t%0.2f \t[mV]\n', R2_relative*( 10^3 ) )
fprintf( 'R3 \t\t= \t%0.2f \t[mV]\n', R3_relative*( 10^3 ) )

fprintf( 'Gm1 \t= \t%0.2f \t[muS]\n', Gm1_relative*( 10^6 ) )
fprintf( 'Gm2 \t= \t%0.2f \t[muS]\n', Gm2_relative*( 10^6 ) )
fprintf( 'Gm3 \t= \t%0.2f \t[muS]\n', Gm3_relative*( 10^6 ) )

fprintf( 'Cm1 \t= \t%0.2f \t[nF]\n', Cm1_relative*( 10^9 ) )
fprintf( 'Cm2 \t= \t%0.2f \t[nF]\n', Cm2_relative*( 10^9 ) )
fprintf( 'Cm3 \t= \t%0.2f \t[nF]\n', Cm3_relative*( 10^9 ) )

fprintf( 'Gna1 \t= \t%0.2f \t[muS]\n', Gna1_relative*( 10^6 ) )
fprintf( 'Gna2 \t= \t%0.2f \t[muS]\n', Gna2_relative*( 10^6 ) )
fprintf( 'Gna3 \t= \t%0.2f \t[muS]\n', Gna3_relative*( 10^6 ) )
fprintf( '\n' )

% Print out synapse information.
fprintf( 'Synapse Parameters:\n' )
fprintf( 'dEs31 \t= \t%0.2f \t[mV]\n', dEs31_relative*( 10^3 ) )
fprintf( 'dEs32 \t= \t%0.2f \t[mV]\n', dEs32_relative*( 10^3 ) )

fprintf( 'gs31 \t= \t%0.2f \t[muS]\n', gs31_relative*( 10^6 ) )
fprintf( 'gs32 \t= \t%0.2f \t[muS]\n', gs32_relative*( 10^6 ) )
fprintf( '\n' )

% Print out the applied current information.
fprintf( 'Applied Current Parameters:\n' )
fprintf( 'Ia1 \t= \t%0.2f \t[nA]\n', current_state1*Ia1_relative*( 10^9 ) )
fprintf( 'Ia2 \t= \t%0.2f \t[nA]\n', current_state2*Ia2_relative*( 10^9 ) )
fprintf( 'Ia3 \t= \t%0.2f \t[nA]\n', Ia3_relative*( 10^9 ) )
fprintf( '\n' )

% Print out design parameters.
fprintf( 'Design Parameters:\n' )
fprintf( 'c1 \t\t= \t%0.2f \t[-]\n', c1_relative )
fprintf( 'c2 \t\t= \t%0.2f \t[-]\n', c2_relative )
fprintf( 'delta \t= \t%0.2f \t[mV]\n', delta_relative*( 10^3 ) )

% Print out ending information.
fprintf( '------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )


%% Create the Relative Division Subnetwork.

% Create an instance of the network class.
network_relative = network_class( network_dt, network_tf );

% Create the relative network components.
[ network_relative.neuron_manager, neuron_IDs_relative ] = network_relative.neuron_manager.create_neurons( 3 );
[ network_relative.synapse_manager, synapse_IDs_relative ] = network_relative.synapse_manager.create_synapses( 2 );
[ network_relative.applied_current_manager, applied_current_IDs_relative ] = network_relative.applied_current_manager.create_applied_currents( 3 );

% Set the relative network neuron parameters.
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs_relative, [ Gna1_relative, Gna2_relative, Gna3_relative ], 'Gna' );
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs_relative, [ R1_relative, R2_relative, R3_relative ], 'R' );
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs_relative, [ Gm1_relative, Gm2_relative, Gm3_relative ], 'Gm' );
network_relative.neuron_manager = network_relative.neuron_manager.set_neuron_property( neuron_IDs_relative, [ Cm1_relative, Cm2_relative, Cm3_relative ], 'Cm' );

% Set the relative network synapse parameters.
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs_relative, [ 1, 2 ], 'from_neuron_ID' );
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs_relative, [ 3, 3 ], 'to_neuron_ID' );
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs_relative, [ gs31_relative, gs32_relative ], 'g_syn_max' );
network_relative.synapse_manager = network_relative.synapse_manager.set_synapse_property( synapse_IDs_relative, [ dEs31_relative, dEs32_relative ], 'dE_syn' );

% Set the relative network applied current parameters.
network_relative.applied_current_manager = network_relative.applied_current_manager.set_applied_current_property( applied_current_IDs_relative, [ 1, 2, 3 ], 'neuron_ID' );
network_relative.applied_current_manager = network_relative.applied_current_manager.set_applied_current_property( applied_current_IDs_relative, [ current_state1*Ia1_relative, current_state2*Ia2_relative, Ia3_relative ], 'I_apps' );


%% Convert the Relative Division Parameters to Fundamental Absolute Division Parameters.

% Convert the maximum membrane voltages.
R1_absolute = R1_relative;                                                                                                  % [V] Maximum Membrane Voltage (Neuron 1).
R2_absolute = R2_relative;                                                                                                  % [V] Maximum Membrane Voltage (Neuron 2).
R3_absolute = R3_relative;                                                                                                  % [V] Maximum Membrane Voltage (Neuron 3).

% Convert the membrane conductances.
Gm1_absolute = Gm1_relative;                                                                                                % [S] Membrane Conductance (Neuron 1).
Gm2_absolute = Gm2_relative;                                                                                                % [S] Membrane Conductance (Neuron 2).
Gm3_absolute = Gm3_relative;                                                                                                % [S] Membrane Conductance (Neuron 3).

% Convert the membrane capacitances.
Cm1_absolute = Cm1_relative;                                                                                                % [F] Membrane Capacitance (Neuron 1).
Cm2_absolute = Cm2_relative;                                                                                                % [F] Membrane Capacitance (Neuron 2).
Cm3_absolute = Cm3_relative;                                                                                                % [F] Membrane Capacitance (Neuron 3).

% Convert the sodium channel conductances.
Gna1_absolute = Gna1_relative;                                                                                              % [S] Sodium Channel Conductance (Neuron 1).
Gna2_absolute = Gna2_relative;                                                                                              % [S] Sodium Channel Conductance (Neuron 2).
Gna3_absolute = Gna3_relative;                                                                                              % [S] Sodium Channel Conductance (Neuron 3).

% Define the synaptic conductances.
dEs31_absolute = 194e-3;                                                                                                    % [V] Synaptic Reversal Potential (Synapse 31).
dEs32_absolute = 0;                                                                                                         % [V] Synaptic Reversal Potential (Synapse 32).

% Convert the applied currents.
Ia1_absolute = Ia1_relative;                                                                                                % [A] Applied Current (Neuron 1).
Ia2_absolute = Ia2_relative;                                                                                                % [A] Applied Current (Neuron 2).
Ia3_absolute = 0;                                                                                                           % [A] Applied Current (Neuron 3).

% Set the network design parameters.
delta_absolute = 1e-3;                                                                                                    	% [V] Membrane Voltage Offset.
c1_absolute = ( delta_absolute*R2_absolute*R3_absolute )/( R1_absolute*R3_absolute - delta_absolute*R1_absolute );          % [V^2] Design Constant 1.


%% Compute the Derived Parameters of the Absolute Division Subnetwork.

% Compute the network design parameters.
c2_absolute = ( c1_absolute*R1_absolute - delta_absolute*R2_absolute )/( delta_absolute );                                                                      % [V] Absolute Division Parameter 2.

% Compute the maximum membrane voltages.
R3_absolute = c1_absolute*R1_absolute/c2_absolute;                                                                                                              % [V] Maximum Membrane Voltage (Neuron 3).

% Compute the synaptic conductances.
gs31_absolute = ( R3_absolute*Gm3_absolute - Ia3_absolute )/( dEs31_absolute - R3_absolute );                                                                  	% [S] Maximum Synaptic Conductance.
gs32_absolute = ( ( dEs31_absolute - delta_absolute )*gs31_absolute + Ia3_absolute - delta_absolute*Gm3_absolute )/( delta_absolute - dEs32_absolute );         % [S] Maximum Synaptic Conductance.


%% Print Absolute Division Subnetwork Parameters.

% Print out a header.
fprintf( '\n------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )
fprintf( 'ABSOLUTE DIVISION SUBNETWORK PARAMETERS:\n' )
fprintf( '------------------------------------------------------------\n' )

% Print out neuron information.
fprintf( 'Neuron Parameters:\n' )
fprintf( 'R1 \t\t= \t%0.2f \t[mV]\n', R1_absolute*( 10^3 ) )
fprintf( 'R2 \t\t= \t%0.2f \t[mV]\n', R2_absolute*( 10^3 ) )
fprintf( 'R3 \t\t= \t%0.2f \t[mV]\n', R3_absolute*( 10^3 ) )

fprintf( 'Gm1 \t= \t%0.2f \t[muS]\n', Gm1_absolute*( 10^6 ) )
fprintf( 'Gm2 \t= \t%0.2f \t[muS]\n', Gm2_absolute*( 10^6 ) )
fprintf( 'Gm3 \t= \t%0.2f \t[muS]\n', Gm3_absolute*( 10^6 ) )

fprintf( 'Cm1 \t= \t%0.2f \t[nF]\n', Cm1_absolute*( 10^9 ) )
fprintf( 'Cm2 \t= \t%0.2f \t[nF]\n', Cm2_absolute*( 10^9 ) )
fprintf( 'Cm3 \t= \t%0.2f \t[nF]\n', Cm3_absolute*( 10^9 ) )

fprintf( 'Gna1 \t= \t%0.2f \t[muS]\n', Gna1_absolute*( 10^6 ) )
fprintf( 'Gna2 \t= \t%0.2f \t[muS]\n', Gna2_absolute*( 10^6 ) )
fprintf( 'Gna3 \t= \t%0.2f \t[muS]\n', Gna3_absolute*( 10^6 ) )
fprintf( '\n' )

% Print out synapse information.
fprintf( 'Synapse Parameters:\n' )
fprintf( 'dEs31 \t= \t%0.2f \t[mV]\n', dEs31_absolute*( 10^3 ) )
fprintf( 'dEs32 \t= \t%0.2f \t[mV]\n', dEs32_absolute*( 10^3 ) )

fprintf( 'gs31 \t= \t%0.2f \t[muS]\n', gs31_absolute*( 10^6 ) )
fprintf( 'gs32 \t= \t%0.2f \t[muS]\n', gs32_absolute*( 10^6 ) )
fprintf( '\n' )

% Print out the applied current information.
fprintf( 'Applied Current Parameters:\n' )
fprintf( 'Ia1 \t= \t%0.2f \t[nA]\n', current_state1*Ia1_absolute*( 10^9 ) )
fprintf( 'Ia2 \t= \t%0.2f \t[nA]\n', current_state2*Ia2_absolute*( 10^9 ) )
fprintf( 'Ia3 \t= \t%0.2f \t[nA]\n', Ia3_absolute*( 10^9 ) )
fprintf( '\n' )

% Print out design parameters.
fprintf( 'Design Parameters:\n' )
fprintf( 'c1 \t\t= \t%0.2f \t[mV]\n', c1_absolute*( 10^3 ) )
fprintf( 'c2 \t\t= \t%0.2f \t[mV]\n', c2_absolute*( 10^3 ) )
fprintf( 'delta \t= \t%0.2f \t[mV]\n', delta_absolute*( 10^3 ) )

% Print out ending information.
fprintf( '------------------------------------------------------------\n' )
fprintf( '------------------------------------------------------------\n' )


%% Create the Absolute Division Subnetwork.

% Create an instance of the network_absolute class.
network_absolute = network_class( network_dt, network_tf );

% Create the network_absolute components.
[ network_absolute.neuron_manager, neuron_IDs_absolute ] = network_absolute.neuron_manager.create_neurons( 3 );
[ network_absolute.synapse_manager, synapse_IDs_absolute ] = network_absolute.synapse_manager.create_synapses( 2 );
[ network_absolute.applied_current_manager, applied_current_IDs_absolute ] = network_absolute.applied_current_manager.create_applied_currents( 3 );

% Set the absolute network neuron parameters.
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ R1_absolute, R2_absolute, R3_absolute ], 'R' );
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ Gm1_absolute, Gm2_absolute, Gm3_absolute ], 'Gm' );
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ Cm1_absolute, Cm2_absolute, Cm3_absolute ], 'Cm' );
network_absolute.neuron_manager = network_absolute.neuron_manager.set_neuron_property( neuron_IDs_absolute, [ Gna1_absolute, Gna2_absolute, Gna3_absolute ], 'Gna' );

% Set the absolute network synapse parameters.
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, [ 1, 2 ], 'from_neuron_ID' );
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, [ 3, 3 ], 'to_neuron_ID' );
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, [ gs31_absolute, gs32_absolute ], 'g_syn_max' );
network_absolute.synapse_manager = network_absolute.synapse_manager.set_synapse_property( synapse_IDs_absolute, [ dEs31_absolute, dEs32_absolute ], 'dE_syn' );

% Set the absolute network applied current parameters.
network_absolute.applied_current_manager = network_absolute.applied_current_manager.set_applied_current_property( applied_current_IDs_absolute, [ 1, 2, 3 ], 'neuron_ID' );
network_absolute.applied_current_manager = network_absolute.applied_current_manager.set_applied_current_property( applied_current_IDs_absolute, [ current_state1*Ia1_absolute, current_state2*Ia2_absolute, Ia3_absolute ], 'I_apps' );


%% Simulate the Relative Division Subnetwork.

% Start the timer.
tic

% Simulate the network.
[ network_relative, ts_relative, Us_relative, hs_relative, dUs_relative, dhs_relative, G_syns_relative, I_leaks_relative, I_syns_relative, I_nas_relative, I_apps_relative, I_totals_relative, m_infs_relative, h_infs_relative, tauhs_relative, neuron_IDs_relative ] = network_relative.compute_set_simulation(  );

% End the timer.
relative_simulation_duration = toc;


%% Simulate the Absolute Division Subnetwork.

% Start the timer.
tic

% Simulate the network.
[ network_absolute, ts_absolute, Us_absolute, hs_absolute, dUs_absolute, dhs_absolute, G_syns_absolute, I_leaks_absolute, I_syns_absolute, I_nas_absolute, I_apps_absolute, I_totals_absolute, m_infs_absolute, h_infs_absolute, tauhs_absolute, neuron_IDs_absolute ] = network_absolute.compute_set_simulation(  );

% End the timer.
absolute_simulation_duration = toc;


%% Plot the Division Subnetwork Results.

% Plot the network currents over time.
fig_relative_network_currents = network_relative.network_utilities.plot_network_currents( ts_relative, I_leaks_relative, I_syns_relative, I_nas_relative, I_apps_relative, I_totals_relative, neuron_IDs_relative );
fig_absolute_network_currents = network_absolute.network_utilities.plot_network_currents( ts_absolute, I_leaks_absolute, I_syns_absolute, I_nas_absolute, I_apps_absolute, I_totals_absolute, neuron_IDs_absolute );

% Plot the network states over time.
fig_relative_network_states = network_relative.network_utilities.plot_network_states( ts_relative, Us_relative, hs_relative, neuron_IDs_relative );
fig_absolute_network_states = network_absolute.network_utilities.plot_network_states( ts_absolute, Us_absolute, hs_absolute, neuron_IDs_absolute );

% Animate the network states over time.
fig_relative_network_animation = network_relative.network_utilities.animate_network_states( Us_relative, hs_relative, neuron_IDs_relative );
fig_absolute_network_animation = network_absolute.network_utilities.animate_network_states( Us_absolute, hs_absolute, neuron_IDs_absolute );

