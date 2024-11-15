%% Relative Addition Subnetwork Error

% Clear Everything.
clear, close('all'), clc


%% Initialize Project Options.

% Define the save and load directories.
save_directory = '.\Save';
load_directory = '.\Load';

% Set a flag to determine whether to simulate.
b_simulate = true;
% b_simulate = false;

% Set the level of verbosity.
b_verbose = true;

% Define the network integration step size.
network_dt = 1e-3;
network_tf = 3;


%% Create Relative Addition Subnetwork.

% Define the network parameters.
num_neurons = 3;                                            % [#] Number of Neurons.
c = 1;                                                      % [-] Subnetwork Gain.

% Create an instance of the network class.
network = network_class( network_dt, network_tf );

% Create an addition subnetwork.
[ network, neuron_IDs_add, synapse_IDs_add, applied_current_IDs_add ] = network.create_relative_addition_subnetwork( num_neurons, c );


%% Compute Desired & Achieved Relative Addition Formulations.

% Retrieve network information.
Rs = cell2mat( network.neuron_manager.get_neuron_property( 'all', 'R' ) );
Cms = cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Cm' ) );
Gms = cell2mat( network.neuron_manager.get_neuron_property( 'all', 'Gm' ) );
Ias = cell2mat( network.neuron_manager.get_neuron_property( 'all', 'I_tonic' ) );
gs = network.get_gsynmaxs( 'all' );
dEs = network.get_dEsyns( 'all' );
dt0 = 1e-6;

% Define the addition subnetwork inputs.
U1s = linspace( 0, Rs( 1 ), 20  );
U2s = linspace( 0, Rs( 2 ), 20  );

% Create an input grid.
[ U1s_grid, U2s_grid ] = meshgrid( U1s, U2s );

% Create the input points.
U1s_flat = reshape( U1s_grid, [ numel( U1s_grid ), 1 ] );
U2s_flat = reshape( U2s_grid, [ numel( U2s_grid ), 1 ] );

% Compute the desired and achieved relative addition steady state output.
U3s_flat_desired_relative = network.compute_desired_relative_addition_steady_state_output( [ U1s_flat, U2s_flat ], Rs, c );
[ U3s_flat_achieved_relative, As, dts, condition_numbers ] = network.achieved_addition_RK4_stability_analysis( U1s_flat, U2s_flat, Cms, Gms, Rs, Ias, gs, dEs, dt0 );

% Convert the flat steady state output results to grids.
dts_grid = reshape( dts, size( U1s_grid ) );
condition_numbers_grid = reshape( condition_numbers, size( U1s_grid ) );
U3s_grid_desired_relative = reshape( U3s_flat_desired_relative, size( U1s_grid ) );
U3s_grid_achieved_relative = reshape( U3s_flat_achieved_relative, size( U1s_grid ) );

% Retrieve the maximum RK4 step size and condition number.
[ dt_max, indexes_dt ] = max( dts );
[ condition_number_max, indexes_condition_number ] = max( condition_numbers );


%% Print the Desired and Achieved Relative Addition Formulation Results.

% Print out the stability information.
fprintf( 'STABILITY SUMMARY:\n' )
fprintf( 'Linearized System Matrix: A =\n\n' ), disp( As( :, :, indexes_condition_number ) )
fprintf( 'Max RK4 Step Size: \t\tdt_max = %0.3e [s] @ ( %0.2f [mV], %0.2f [mV] )\n', dt_max, U1s_flat( indexes_dt )*( 10^3 ), U2s_flat( indexes_dt )*( 10^3 ) )
fprintf( 'Proposed Step Size: \tdt = %0.3e [s]\n', network_dt )
fprintf( 'Condition Number: \t\tcond( A ) = %0.3e [-] @ ( %0.2f [mV], %0.2f [mV] )\n', condition_number_max, U1s_flat( indexes_condition_number )*( 10^3 ), U2s_flat( indexes_condition_number )*( 10^3 ) )


%% Plot the Desired and Achieved Relative Addition Formulation Results.

% Plot the desired and achieved relative addition formulation results.
fig = figure( 'Color', 'w', 'Name', 'Relative Addition Theory' ); hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage 1 (Input), U1 [mV]' ), ylabel( 'Membrane Voltage 2 (Input), U2 [mV]' ), zlabel( 'Membrane Voltage 3 (Output), U3 [mV]' ), title( 'Relative Addition Theory' )
surf( U1s_grid, U2s_grid, U3s_grid_desired_relative, 'Facecolor', 'b', 'Edgecolor', 'None' )
surf( U1s_grid, U2s_grid, U3s_grid_achieved_relative, 'Facecolor', 'r', 'Edgecolor', 'None' )
legend( 'Desired', 'Achieved' )
saveas( fig, [ save_directory, '\', 'relative_addition_theory' ] )

% Plot the RK4 maximum timestep.
fig = figure( 'Color', 'w', 'Name', 'Relative Addition RK4 Maximum Timestep' ); hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage 1 (Input), U1 [mV]' ), ylabel( 'Membrane Voltage 2 (Input), U2 [mV]' ), zlabel( 'Rk4 Maximum Timestep, dt [s]' ), title( 'Relative Addition RK4 Maximum Timestep' )
surf( U1s_grid, U2s_grid, dts_grid, 'Edgecolor', 'None' )
saveas( fig, [ save_directory, '\', 'relative_addition_rk4_maximum_timestep' ] )

% Plot the linearized system condition numbers.
fig = figure( 'Color', 'w', 'Name', 'Relative Addition Condition Numbers' ); hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage 1 (Input), U1 [mV]' ), ylabel( 'Membrane Voltage 2 (Input), U2 [mV]' ), zlabel( 'Condition Number [-]' ), title( 'Relative Addition Condition Number' )
surf( U1s_grid, U2s_grid, condition_numbers_grid, 'Edgecolor', 'None' )
saveas( fig, [ save_directory, '\', 'relative_addition_condition_numbers' ] )


%% Simulate the Network.

% Determine whether to simulate the network.
if b_simulate               % If we want to simulate the network....
    
    % Define the number of applied currents to use.
    n_applied_currents1 = 15;
    n_applied_currents2 = 10;
    
    % Create the applied currents.
    applied_currents1 = linspace( 0, 20e-9, n_applied_currents1 );
    applied_currents2 = linspace( 0, 20e-9, n_applied_currents2 );
    
    % Create a grid of the applied currents.
    [ Applied_Currents1, Applied_Currents2 ] = meshgrid( applied_currents1, applied_currents2 );
    
    % Create a matrix to store the membrane voltages.
    Us_achieved = zeros( n_applied_currents2, n_applied_currents1, 3 );
    
    % Simulate the network for each of the applied current combinations.
    for k1 = 1:n_applied_currents1                          % Iterate through each of the currents applied to the first neuron...
        for k2 = 1:n_applied_currents2                      % Iterate through each of the currents applied to the second neuron...
            
            % Create applied currents.
            network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_add( 1 ), Applied_Currents1( k2, k1 ), 'I_apps' );
            network.applied_current_manager = network.applied_current_manager.set_applied_current_property( applied_current_IDs_add( 2 ), Applied_Currents2( k2, k1 ), 'I_apps' );
            
            % Simulate the network.
            [ network, ts, Us, hs, dUs, dhs, G_syns, I_leaks, I_syns, I_nas, I_apps, I_totals, m_infs, h_infs, tauhs, neuron_IDs ] = network.compute_set_simulation(  );
            
            % Retrieve the final membrane voltages.
            Us_achieved( k2, k1, : ) = Us( :, end );
            
        end
    end

    % Save the simulation results.
    save( [ save_directory, '\', 'relative_addition_subnetwork_error' ], 'Applied_Currents1', 'Applied_Currents2', 'Us_achieved' )
    
else                % Otherwise... (We must want to load data from an existing simulation...
    
    % Load the simulation results.
    data = load( [ load_directory, '\', 'relative_addition_subnetwork_error' ] );
    
    % Store the simulation results in separate variables.
    Applied_Currents1 = data.Applied_Currents1;
    Applied_Currents2 = data.Applied_Currents2;
    Us_achieved = data.Us_achieved;

end


%% Plot the Network Results.

% Get the activation domains of the neurons.
R1 = network.neuron_manager.get_neuron_property( 1, 'R' ); R1 = R1{1};
R2 = network.neuron_manager.get_neuron_property( 2, 'R' ); R2 = R2{1};
R3 = network.neuron_manager.get_neuron_property( 3, 'R' ); R3 = R3{1};

% Compute the desired membrane voltage output.
Us_desired = cat( 3, Us_achieved( :, :, 1 ), Us_achieved( :, :, 2 ), ( R3/2 )*( Us_achieved( :, :, 1 )/R1 + Us_achieved( :, :, 2 )/R2 ) );

% Compute the error between the achieved and desired results.
error = Us_achieved( :, :, end ) - Us_desired( :, :, end );

% Compute the mean squared error summary statistic.
mse = sqrt( sum( error.^2, 'all' ) );

% Create a surface that shows the desired membrane voltage output.
fig = figure( 'Color', 'w', 'Name', 'Relative Addition Subnetwork Steady State Response (Desired)' ); hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [V]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [V]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [V]' ), title( 'Relative Addition Subnetwork Steady State Response (Desired)' )
surf( Us_desired( :, :, 1 ), Us_desired( :, :, 2 ), Us_desired( :, :, 3 ), 'Edgecolor', 'None' )
saveas( fig, [ save_directory, '\', 'relative_addition_ss_response_desired' ] )

% Create a surface that shows the achieved membrane voltage output.
fig = figure( 'Color', 'w', 'Name', 'Relative Addition Subnetwork Steady State Response (Achieved)' ); hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [V]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [V]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [V]' ), title( 'Relative Addition Subnetwork Steady State Response (Achieved)' )
surf( Us_achieved( :, :, 1 ), Us_achieved( :, :, 2 ), Us_achieved( :, :, 3 ), 'Edgecolor', 'None' )
saveas( fig, [ save_directory, '\', 'relative_addition_ss_response_achieved' ] )

% Create a figure that shows the differences between the achieved and desired membrane voltage outputs.
fig = figure( 'Color', 'w', 'Name', 'Relative Addition Subnetwork Steady State Response (Comparison)' ); hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [V]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [V]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [V]' ), title( 'Relative Addition Subnetwork Steady State Response (Comparison)' )
surf( Us_desired( :, :, 1 ), Us_desired( :, :, 2 ), Us_desired( :, :, 3 ), 'Edgecolor', 'None', 'Facecolor', 'b' )
surf( U1s_grid, U2s_grid, U3s_grid_achieved_relative, 'Edgecolor', 'None', 'Facecolor', 'g' )
surf( Us_achieved( :, :, 1 ), Us_achieved( :, :, 2 ), Us_achieved( :, :, 3 ), 'Edgecolor', 'None', 'Facecolor', 'r' )
legend( 'Desired', 'Achieved (Theory)', 'Achieved (Numerical)' )
saveas( fig, [ save_directory, '\', 'relative_addition_ss_response_comparison' ] )

% Create a surface that shows the membrane voltage error.
fig = figure( 'Color', 'w', 'Name', 'Relative Addition Subnetwork Steady State Error' ); hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [V]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [V]' ), zlabel( 'Membrane Voltage Error, E [V]' ), title( 'Relative Addition Subnetwork Steady State Error' )
surf( Us_achieved( :, :, 1 ), Us_achieved( :, :, 2 ), error, 'Edgecolor', 'None' )
saveas( fig, [ save_directory, '\', 'relative_addition_ss_response_error' ] )


% % Plot the network currents over time.
% fig_network_currents = network.network_utilities.plot_network_currents( ts, I_leaks, I_syns, I_nas, I_apps, I_totals, neuron_IDs );
%
% % Plot the network states over time.
% fig_network_states = network.network_utilities.plot_network_states( ts, Us, hs, neuron_IDs );
%
% % Animate the network states over time.
% fig_network_animation = network.network_utilities.animate_network_states( Us, hs, neuron_IDs );

