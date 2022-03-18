%% Multistate CPG Phase Alignment

% Clear Everything.
clear, close('all'), clc


%% Initialize Project Options.

% Set the level of verbosity.
b_verbose = true;

% Define the path to the directory that contains the robot data.
robot_data_save_path = 'D:\Documents\GitHub\Quadruped_Robot\Code\Matlab\Quadruped_Simulation\Quadruped_Simulation_Framework\Main\Simulation\Software\Multistate_CPG_Phase_Alignment\Robot_Data\Save';
robot_data_load_path = 'D:\Documents\GitHub\Quadruped_Robot\Code\Matlab\Quadruped_Simulation\Quadruped_Simulation_Framework\Main\Simulation\Software\Multistate_CPG_Phase_Alignment\Robot_Data\Load';

% robot_data_save_path = 'C:\Users\Cody Scharzenberger\Documents\GitHub\Quadruped_Robot\Code\Matlab\Quadruped_Simulation\Quadruped_Simulation_Framework\Main\Simulation\Software\Multistate_CPG_Phase_Alignment\Robot_Data\Save';
% robot_data_load_path = 'C:\Users\Cody Scharzenberger\Documents\GitHub\Quadruped_Robot\Code\Matlab\Quadruped_Simulation\Quadruped_Simulation_Framework\Main\Simulation\Software\Multistate_CPG_Phase_Alignment\Robot_Data\Load';

% Define the network integration step size.
% network_dt = 1e-3;
network_dt = 1e-4;
% network_dt = 1e-6;
network_tf = 1;
% network_tf = 3;


%% Initialize the Data Loader Class.

% Determine whether to print status messages.
if b_verbose, fprintf( 'INITIALIZING DATA LOADER. Please Wait...\n' ), end

% Start a timer.
tic

% Create an instance of the data loader class.
data_loader = data_loader_utilities_class( robot_data_load_path );

% Retrieve the elapsed time.
elapsed_time = toc;

% Determine whether to print status messages.
if b_verbose, fprintf( 'INITIALIZING DATA LOADER. Please Wait... Done. %0.3f [s] \n\n', elapsed_time ), end


%% Create Multistate CPG Subnetworks.

% Create an instance of the network class.
network = network_class( network_dt, network_tf );

% Define the oscillatory and bistable delta CPG synapse design parameters.
delta_oscillatory = 0.01e-3;                    % [-] Relative Inhibition Parameter for Oscillatory Connections
delta_bistable = -10e-3;                        % [-] Relative Inhibition Parameter for Bistable Connections

% Create the first multistate cpg subnetwork.
[ network, neuron_IDs_cpg1, synapse_IDs_cpg1, applied_current_ID_cpg1 ] = network.create_multistate_cpg_subnetwork( 4, delta_oscillatory, delta_bistable );

% Create the second multistate cpg subnetwork.
[ network, neuron_IDs_cpg2, synapse_IDs_cpg2, applied_current_ID_cpg2 ] = network.create_multistate_cpg_subnetwork( 4, delta_oscillatory, delta_bistable );

% Disable the cpg subnetworks.
network.neuron_manager = network.neuron_manager.disable_neurons( neuron_IDs_cpg1 );
network.neuron_manager = network.neuron_manager.disable_neurons( neuron_IDs_cpg2 );


%% Create Addition Subnetwork.

% Create an addition subnetwork.
[ network, neuron_IDs_add, synapse_IDs_add ] = network.create_addition_subnetwork(  );

% Create applied currents.
[ network.applied_current_manager, applied_current_IDs_add ] = network.applied_current_manager.create_applied_currents( 2 );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_add(1), 5e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_add(2), 10e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_add, neuron_IDs_add(1:2), 'neuron_ID' );

% Disable the addition subnetwork.
network.neuron_manager = network.neuron_manager.disable_neurons( neuron_IDs_add );


%% Create Subtraction Subnetwork.

% Create a subtraction subnetwork.
[ network, neuron_IDs_sub, synapse_IDs_sub ] = network.create_subtraction_subnetwork(  );

% Create applied currents.
[ network.applied_current_manager, applied_current_IDs_sub ] = network.applied_current_manager.create_applied_currents( 2 );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_sub(1), 15e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_sub(2), 10e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_sub, neuron_IDs_sub(1:2), 'neuron_ID' );

% % Set the applied current properties.
% ts = 0:network.dt:network.tf;
% I_apps1 = (15e-9)*ts;
% I_apps2 = max( (15e-9)*( ts - 0.1 ), 0 );
% 
% % Create applied currents.
% [ network.applied_current_manager, applied_current_IDs_sub ] = network.applied_current_manager.create_applied_currents( 2 );
% network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_sub, { ts }, 'ts' );
% network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_sub, { I_apps1 I_apps1 }, 'I_apps' );
% % network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_sub, { I_apps1 I_apps2 }, 'I_apps' );
% network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_sub, neuron_IDs_sub(1:2), 'neuron_ID' );


% % Modify the subtraction subnetwork neurons.
% network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs_sub( 1:2 ), [ 0.95e-6, 0.95e-6 ], 'Gm' );
% network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs_sub( 1:2 ), [ 47.5e-9, 950e-9 ], 'Cm' );


% network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs_sub( 1:2 ), [ 0.1e-6, 0.1e-6 ], 'Gm' );
% network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs_sub( 1:3 ), [ 0.1e-6, 0.1e-6, 0.1e-6 ], 'Gm' );
% network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs_sub( 1:2 ), [ 0.95e-6, 0.95e-6 ], 'Gm' );
% network.neuron_manager = network.neuron_manager.set_neuron_property( neuron_IDs_sub( 1:3 ), [ 0.95e-6, 0.95e-6, 0.95e-6 ], 'Gm' );

% Disable the subtraction subnetwork.
network.neuron_manager = network.neuron_manager.disable_neurons( neuron_IDs_sub );


%% Create Division Subnetwork.

% Create a division subnetwork.
[ network, neuron_IDs_div, synapse_IDs_div ] = network.create_division_subnetwork(  );

% Create applied currents.
[ network.applied_current_manager, applied_current_IDs_div ] = network.applied_current_manager.create_applied_currents( 2 );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_div(1), 15e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_div(2), 5e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_div, neuron_IDs_div(1:2), 'neuron_ID' );

% Disable the division subnetwork.
network.neuron_manager = network.neuron_manager.disable_neurons( neuron_IDs_div );


%% Create Multiplication Subnetwork.

% Create a multiplication subnetwork.
[ network, neuron_IDs_mult, synapse_IDs_mult, applied_current_ID_mult ] = network.create_multiplication_subnetwork(  );

% Create applied currents for the multiplication subnetwork.
[ network.applied_current_manager, applied_current_IDs_mult ] = network.applied_current_manager.create_applied_currents( 2 );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_mult(1), 10e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_mult(2), 30e-9, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_mult, neuron_IDs_mult(1:2), 'neuron_ID' );

% Disable the multiplication subnetwork.
network.neuron_manager = network.neuron_manager.disable_neurons( neuron_IDs_mult );


%% Create Derivation Subnetwork.

% Create a derivation subnetwork.
[ network, neuron_IDs_der, synapse_IDs_der ] = network.create_derivation_subnetwork(  );

% Define the properties of the applied currents for the derivation subnetwork.
I_mag = 20e-9;
ts = ( 0:network.dt:network.tf )';
I_apps = I_mag*ts;
num_timesteps = length( ts );

% Create applied currents for the derivation subnetwork.
[ network.applied_current_manager, applied_current_IDs_der ] = network.applied_current_manager.create_applied_currents( 2 );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_der, { 'Der1', 'Der2' }, 'name' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_der, neuron_IDs_der( 1:2 ), 'neuron_ID' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_der, { ts }, 'ts' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_der, { I_apps }, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_der, num_timesteps, 'num_timesteps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_der, network.dt, 'dt' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_der, network.tf, 'tf' );

% Disable the derivation subnetwork.
network.neuron_manager = network.neuron_manager.disable_neurons( neuron_IDs_der );


%% Create an Integration Subnetwork.

% Define the integration constants.
ki_mean = 0.5e9;
ki_range = 0.5e9;

% Create an integration subnetwork.
[ network, neuron_IDs_int, synapse_IDs_int, applied_current_IDs_int ] = network.create_integration_subnetwork( ki_mean, ki_range );

% Define the properties of the applied currents for the derivation subnetwork.
Gm = cell2mat( network.neuron_manager.get_neuron_property( neuron_IDs_int(1), 'Gm' ) );
R = cell2mat( network.neuron_manager.get_neuron_property( neuron_IDs_int(1), 'R' ) );
ts = ( 0:network.dt:network.tf )';
num_timesteps = length( ts );
I_mag1 = 20e-9;
I_mag2 = 1e-9;
I_apps = I_mag1*ones( num_timesteps, 1 ) + I_mag2*( ts > 0.5 & ts < 0.75 );

% Create applied currents for the integration subnetwork.
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_int(1), { ts }, 'ts' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_int(1), { I_apps }, 'I_apps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_int(1), num_timesteps, 'num_timesteps' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_int(1), network.dt, 'dt' );
network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_int(1), network.tf, 'tf' );

% Disable the integration subnetwork.
network.neuron_manager = network.neuron_manager.disable_neurons( applied_current_IDs_int );



%% DEBUGGING CODE


%% Simulate the Network.

% Simulate the network.
[ network, ts, Us, hs, dUs, dhs, G_syns, I_leaks, I_syns, I_nas, I_apps, I_totals, m_infs, h_infs, tauhs, neuron_IDs ] = network.compute_set_simulation(  );


%% Plot the Network Results.

% Plot the network currents over time.
fig_network_currents = network.network_utilities.plot_network_currents( ts, I_leaks, I_syns, I_nas, I_apps, I_totals, neuron_IDs );

% Plot the network states over time.
fig_network_states = network.network_utilities.plot_network_states( ts, Us, hs, neuron_IDs );
% fig_network_states = network.network_utilities.plot_network_states( ts, Us(1:2, :), hs(1:2, :), neuron_IDs );

% dUs = Us( 1, : ) - Us( 2, : );
% figure( 'Color', 'w' ), plot( ts, dUs, '-m', 'Linewidth', 3 )

% Animate the network states over time.
fig_network_animation = network.network_utilities.animate_network_states( Us, hs, neuron_IDs );
% fig_network_animation = network.network_utilities.animate_network_states( Us(1:2, :), hs(1:2, :), neuron_IDs );


x = 1;

