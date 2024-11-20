%% Relative Transmission Subnetwork Error.

% Clear Everything.
clear, close( 'all' ), clc


%% Define Simulation Parameters.

% Define the save and load directories.
save_directory = '.\Save';                      	% [str] Save Directory.
load_directory = '.\Load';                      	% [str] Load Directory.

% Set a flag to determine whether to simulate.
simulate_flag = true;                             	% [T/F] Simulation Flag. (Determines whether to create a new simulation of the steady state error or to load a previous simulation.)
% simulate_flag = false;                            % [T/F] Simulation Flag. (Determines whether to create a new simulation of the steady state error or to load a previous simulation.)

% Set the level of verbosity.
verbose_flag = true;                            	% [T/F] Printing Flag. (Determines whether to print out information.)

% Define the undetected option.
undetected_option = 'Error';                        % [str] Undetected Option.

% Define the network simulation timestep.
network_dt = 1e-3;                                 	% [s] Simulation Timestep.
% network_dt = 1e-4;                             	% [s] Simulation Timestep.

% Define the network simulation duration.
network_tf = 0.5;                                 	% [s] Simulation Duration.
% network_tf = 1;                                 	% [s] Simulation Duration.
% network_tf = 3;                                 	% [s] Simulation Duration.

% Define the number of neurons.
num_neurons = 2;                                 	% [#] Number of Neurons.

% Define the integration method.
integration_method = 'RK4';                         % [str] Integration Method (Either FE for Forward Euler or RK4 for Fourth Order Runge-Kutta).

% Define the encoding scheme.
encoding_scheme = 'relative';


%% Define Relative Transmission Subnetwork Parameters.

% Define the transmission subnetwork design parameters.
R1 = 20e-3;                                         % [V] Maximum Membrane Voltage (Neuron 1).
R2 = 20e-3;                                         % [V] Maximum Membrane Voltage (Neuron 2).
Gm1 = 1e-6;                                         % [S] Membrane Conductance (Neuron 1).
Gm2 = 1e-6;                                         % [S] Membrane Conductance (Neuron 2).
Cm1 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 1).
Cm2 = 5e-9;                                         % [F] Membrane Capacitance (Neuron 2).

% Store the transmission subnetwork design parameters in a cell.
transmission_parameters = { R1, R2, Gm1, Gm2, Cm1, Cm2 };


%% Define Encoding & Decoding Operations.

% Define the absolute transmission comparison example.
R1_absolute = 20e-3;                                % [V] Relative Maximum Membrane Voltage (Neuron 1).  (Used for decoding.)
R2_absolute = 20e-3;                                % [V] Relative Maximum Membrane Voltage (Neuron 2).  (Used for decoding.)

% Define the decoding operations.
f_decode1 = @( x ) ( R1_absolute/R1 )*x*( 10^3 );
f_decode2 = @( x ) ( R2_absolute/R2 )*x*( 10^3 );
f_decode = @( x ) [ f_decode1( x( :, 1 ) ), f_decode2( x( :, 2 ) ) ];


%% Define the Relative Transmission Subnetwork Input Current Parameters.

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


%% Create the Relative Transmission Subnetwork.

% Create an instance of the netwo5rk class.
network = network_class( network_dt, network_tf );

% Create a transmission subnetwork.
[ c, Gnas, R2, dEs21, gs21, Ia2, neurons, synapses, neuron_manager, synapse_manager, network ] = network.create_transmission_subnetwork( transmission_parameters, encoding_scheme, network.neuron_manager, network.synapse_manager, network.applied_current_manager, true, true, false, undetected_option );

% Create the input applied current.
[ ~, ~, ~, network.applied_current_manager ] = network.applied_current_manager.create_applied_current( input_current_ID, input_current_name, input_current_to_neuron_ID, ts, Ias1, true, network.applied_current_manager.applied_currents, true, false, network.applied_current_manager.array_utilities );


%% Print Relative Transmission Subnetwork Information.

% Print transmission subnetwork information.
network.print( network.neuron_manager, network.synapse_manager, network.applied_current_manager, verbose_flag );


%% Compute Desired & Achieved Relative Transmission Formulations.

% Define the property retrieval settings.
as_matrix_flag = true;

% Retrieve properties from the existing network.
Cms = network.neuron_manager.get_neuron_property( 'all', 'Cm', as_matrix_flag, network.neuron_manager.neurons, undetected_option );         % [F] Membrane Capacitance.
Gms = network.neuron_manager.get_neuron_property( 'all', 'Gm', as_matrix_flag, network.neuron_manager.neurons, undetected_option );         % [S] Membrane Conductance.
Rs = network.neuron_manager.get_neuron_property( 'all', 'R', as_matrix_flag, network.neuron_manager.neurons, undetected_option );           % [V] Maximum Membrane Voltage.
gs = network.get_gs( 'all', network.neuron_manager, network.synapse_manager );                                                              % [S] Synaptic Conductance.
dEs = network.get_dEs( 'all', network.neuron_manager, network.synapse_manager );                                                            % [V] Synaptic Reversal Potential.
Ias = network.neuron_manager.get_neuron_property( 'all', 'Itonic', as_matrix_flag, network.neuron_manager.neurons, undetected_option );     % [A] Applied Currents.
dt0 = 1e-6;                                                                                                                                 % [s] Numerical Stability Time Step.

% Define the transmission subnetwork inputs.
U1s = linspace( 0, Rs( 1 ), 100  );

% Create the input points.
U1s_flat = reshape( U1s, [ numel( U1s ), 1 ] );

% Compute the desired and achieved absolute transmission steady state output.
U2s_flat_desired = network.compute_dr_transmission_sso( U1s_flat, c, R1, R2, network.neuron_manager, undetected_option, network.network_utilities );
[ U2s_flat_achieved_theoretical, As, dts, condition_numbers ] = network.achieved_transmission_RK4_stability_analysis( U1s_flat, Cms, Gms, Rs, Ias, gs, dEs, dt0, network.neuron_manager, network.synapse_manager, network.applied_current_manager, undetected_option, network.network_utilities );

% Store the desired and theoretically achieved absolute transmission steady state results in arrays.
Us_flat_desired = [ U1s_flat, U2s_flat_desired ];
Us_flat_achieved_theoretical = [ U1s_flat, U2s_flat_achieved_theoretical ];

% Retrieve the maximum RK4 step size and condition number.
[ dt_max, indexes_dt ] = max( dts );
[ condition_number_max, indexes_condition_number ] = max( condition_numbers );


%% Print the Numerical Stability Information.

% Print out the stability information.
network.numerical_method_utilities.print_numerical_stability_info( As, dts, network_dt, condition_numbers );


%% Plot the Desired and Achieved Relative Transmission Formulation Results.

% Decode the input and output membrane voltages.
Us_flat_desired_decoded = f_decode( Us_flat_desired );
Us_flat_achieved_theoretical_decoded = f_decode( Us_flat_achieved_theoretical );

% Plot the desired and achieved relative transmission formulation results.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Theory' ); hold on, grid on, xlabel( 'Membrane Voltage 1 (Input), U1 [mV]' ), ylabel( 'Membrane Voltage 2 (Output), U2 [mV]' ), title( 'Relative Transmission Theory' )
plot( Us_flat_desired( :, 1 )*( 10^3 ), Us_flat_desired( :, 2 )*( 10^3 ), '-', 'Linewidth', 3 )
plot( Us_flat_achieved_theoretical( :, 1 )*( 10^3 ), Us_flat_achieved_theoretical( :, 2 )*( 10^3 ), '--', 'Linewidth', 3 )
legend( 'Desired', 'Achieved (Theory)' )
saveas( fig, [ save_directory, '\', 'relative_transmission_theory' ] )

% Plot the decoded desired and achieved relative transmission formulation results.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Theory Decoded' ); hold on, grid on, xlabel( 'Input, x [-]' ), ylabel( 'Output, y [-]' ), title( 'Relative Transmission Theory Decoded' )
plot( Us_flat_desired_decoded( :, 1 )*( 10^3 ), Us_flat_desired_decoded( :, 2 )*( 10^3 ), '-', 'Linewidth', 3 )
plot( Us_flat_achieved_theoretical_decoded( :, 1 )*( 10^3 ), Us_flat_achieved_theoretical_decoded( :, 2 )*( 10^3 ), '--', 'Linewidth', 3 )
legend( 'Desired', 'Achieved (Theory)' )
saveas( fig, [ save_directory, '\', 'relative_transmission_theory_decoded' ] )

% Plot the RK4 maximum timestep.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission RK4 Maximum Timestep' ); hold on, grid on, xlabel( 'Membrane Voltage 1 (Input), U1 [mV]' ), ylabel( 'RK4 Maximum Timestep, dt [s]' ), title( 'Relative Transmission RK4 Maximum Timestep' )
plot( Us_flat_desired( :, 1 )*( 10^3 ), dts, '-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'relative_transmission_rk4_maximum_timestep' ] )

% Plot the linearized system condition numbers.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Condition Numbers' ); hold on, grid on, xlabel( 'Membrane Voltage 1 (Input), U1 [mV]' ), ylabel( 'Condition Number [-]' ), title( 'Relative Transmission Condition Number' )
plot( Us_flat_desired( :, 1 )*( 10^3 ), condition_numbers, '-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'relative_transmission_condition_numbers' ] )


%% Simulate the Relative Transmission Network.

% Set additional simulation properties.
filter_disabled_flag = true;                % [T/F] Filter Disabled Flag.
set_flag = true;                            % [T/F] Set Flag.
process_option = 'None';                    % [str] Process Option.
undetected_option = 'Ignore';               % [str] Undetected Option.

% Determine whether to simulate the network.
if simulate_flag                            % If we want to simulate the network...
    
    % Define the number of applied currents to use.
    n_applied_currents = 20;             	% [#] Number of Applied Currents.
    
    % Create the applied currents.
    applied_currents = linspace( 0, network.neuron_manager.neurons( 1 ).R*network.neuron_manager.neurons( 1 ).Gm, n_applied_currents );

    % Create a matrix to store the membrane voltages.
    Us_achieved_numerical = zeros( n_applied_currents, num_neurons );
    
    % Simulate the network for each of the applied current combinations.
    for k = 1:n_applied_currents         	% Iterate through each of the currents applied to the input neuron...
            
            % Create applied currents.
            [ ~, network.applied_current_manager ] = network.applied_current_manager.set_applied_current_property( input_current_ID, applied_currents( k ), 'Ias', network.applied_current_manager.applied_currents, set_flag );

            % Simulate the network.            
            [ ts, Us, hs, dUs, dhs, Gs, I_leaks, I_syns, I_nas, I_apps, I_totals, m_infs, h_infs, tauhs, neurons, synapses, neuron_manager, synapse_manager, network ] = network.compute_simulation( network_dt, network_tf, integration_method, network.neuron_manager, network.synapse_manager, network.applied_current_manager, network.applied_voltage_manager, filter_disabled_flag, set_flag, process_option, undetected_option, network.network_utilities );

            % Retrieve the final membrane voltages.
            Us_achieved_numerical( k, : ) = Us( :, end );
            
    end

    % Save the simulation results.
    save( [ save_directory, '\', 'relative_transmission_subnetwork_error' ], 'applied_currents', 'Us_achieved_numerical' )
    
else                % Otherwise... ( We must want to load data from an existing simulation... )
    
    % Load the simulation results.
    data = load( [ load_directory, '\', 'relative_transmission_subnetwork_error' ] );
    
    % Store the simulation results in separate variables.
    applied_currents = data.applied_currents;
    Us_achieved_numerical = data.Us_achieved_numerical;

end


%% Compute the Relative Transmission Network Error.

% Compute the desired membrane voltage output.
Us_desired_output = network.compute_dr_transmission_sso( Us_achieved_numerical( :, 1 ), c, Rs( 1 ), Rs( 2 ), network.neuron_manager, undetected_option, network.network_utilities );
Us_achieved_theoretical_output = network.compute_achieved_transmission_sso( Us_achieved_numerical( :, 1 ), Rs( 1 ), Gms( 2 ), Ias( 2 ), gs( 2, 1 ), dEs( 2, 1 ), network.neuron_manager, network.synapse_manager, network.applied_current_manager, undetected_option, network.network_utilities );

% Compute the desired membrane voltage output.
Us_desired = Us_achieved_numerical; Us_desired( :, end ) = Us_desired_output;
Us_achieved_theoretical = Us_achieved_numerical; Us_achieved_theoretical( :, end ) = Us_achieved_theoretical_output;

% Decode the achieved and desired decoded membrane voltage output.
R2_decoded = f_decode2( R2 );
Us_desired_decoded = f_decode( Us_desired );
Us_achieved_theoretical_decoded = f_decode( Us_achieved_theoretical );
Us_achieved_numerical_decoded = f_decode( Us_achieved_numerical );

% Compute the summary statistics associated with the decoded and non-decoded theoretical and achieved results.
[ errors_theoretical, error_percentages_theoretical, error_rmse_theoretical, error_rmse_percentage_theoretical, error_std_theoretical, error_std_percentage_theoretical, error_min_theoretical, error_min_percentage_theoretical, index_min_theoretical, error_max_theoretical, error_max_percentage_theoretical, index_max_theoretical, error_range_theoretical, error_range_percentage_theoretical ] = network.numerical_method_utilities.compute_summary_statistics( Us_achieved_theoretical, Us_desired, R2 );
[ errors_numerical, error_percentages_numerical, error_rmse_numerical, error_rmse_percentage_numerical, error_std_numerical, error_std_percentage_numerical, error_min_numerical, error_min_percentage_numerical, index_min_numerical, error_max_numerical, error_max_percentage_numerical, index_max_numerical, error_range_numerical, error_range_percentage_numerical ] = network.numerical_method_utilities.compute_summary_statistics( Us_achieved_numerical, Us_desired, R2 );
[ errors_theoretical_decoded, error_percentages_theoretical_decoded, error_rmse_theoretical_decoded, error_rmse_percentage_theoretical_decoded, error_std_theoretical_decoded, error_std_percentage_theoretical_decoded, error_min_theoretical_decoded, error_min_percentage_theoretical_decoded, index_min_theoretical_decoded, error_max_theoretical_decoded, error_max_percentage_theoretical_decoded, index_max_theoretical_decoded, error_range_theoretical_decoded, error_range_percentage_theoretical_decoded ] = network.numerical_method_utilities.compute_summary_statistics( Us_achieved_theoretical_decoded, Us_desired_decoded, R2_decoded );
[ errors_numerical_decoded, error_percentages_numerical_decoded, error_rmse_numerical_decoded, error_rmse_percentage_numerical_decoded, error_std_numerical_decoded, error_std_percentage_numerical_decoded, error_min_numerical_decoded, error_min_percentage_numerical_decoded, index_min_numerical_decoded, error_max_numerical_decoded, error_max_percentage_numerical_decoded, index_max_numerical_decoded, error_range_numerical_decoded, error_range_percentage_numerical_decoded ] = network.numerical_method_utilities.compute_summary_statistics( Us_achieved_numerical_decoded, Us_desired_decoded, R2_decoded );


%% Print the Relative Tranmission Summary Statistics.

% Retrieve the membrane voltage associated min and max theoretical and numerical error.
Us_critmin_achieved_numerical_steady = Us_achieved_numerical( index_min_theoretical, : );
Us_critmax_achieved_numerical_steady = Us_achieved_numerical( index_max_theoretical, : );
Us_critmin_achieved_theoretical_steady = Us_achieved_theoretical( index_min_theoretical, : );
Us_critmax_achieved_theoretical_steady = Us_achieved_theoretical( index_max_theoretical, : );

% Retrieve the decoded result associated min and max theoretical and numerical error.
Us_critmin_achieved_numerical_steady_decoded = Us_achieved_numerical( index_min_theoretical_decoded, : );
Us_critmax_achieved_numerical_steady_decoded = Us_achieved_numerical( index_max_theoretical_decoded, : );
Us_critmin_achieved_theoretical_steady_decoded = Us_achieved_theoretical( index_min_theoretical_decoded, : );
Us_critmax_achieved_theoretical_steady_decoded = Us_achieved_theoretical( index_max_theoretical_decoded, : );

% Define the membrane voltage summary statistic printing information.
header_mv = 'Relative Transmission Summary Statistics (Membrane Voltages)\n';
unit_mv = 'mV';
scale_mv = 10^3;

% Print the summary statistics for the membrane voltage results.
network.numerical_method_utilities.print_summary_statistics( header_mv, unit_mv, scale_mv, error_rmse_theoretical, error_rmse_percentage_theoretical, error_rmse_numerical, error_rmse_percentage_numerical, error_std_theoretical, error_std_percentage_theoretical, error_std_numerical, error_std_percentage_numerical, error_min_theoretical, error_min_percentage_theoretical, Us_critmin_achieved_theoretical_steady, error_min_numerical, error_min_percentage_numerical, Us_critmin_achieved_numerical_steady, error_max_theoretical, error_max_percentage_theoretical, Us_critmax_achieved_theoretical_steady, error_max_numerical, error_max_percentage_numerical, Us_critmax_achieved_numerical_steady, error_range_theoretical, error_range_percentage_theoretical, error_range_numerical, error_range_percentage_numerical ) 

% Define the membrane voltage summary statistic printing information.
header_decoded = 'Relative Transmission Summary Statistics (Decoded)\n';
unit_decoded = '-';
scale_decoded = 1;

% Print the summary statistics for the decoded results.
network.numerical_method_utilities.print_summary_statistics( header_decoded, unit_decoded, scale_decoded, error_rmse_theoretical_decoded, error_rmse_percentage_theoretical_decoded, error_rmse_numerical_decoded, error_rmse_percentage_numerical_decoded, error_std_theoretical_decoded, error_std_percentage_theoretical_decoded, error_std_numerical_decoded, error_std_percentage_numerical_decoded, error_min_theoretical_decoded, error_min_percentage_theoretical_decoded, Us_critmin_achieved_theoretical_steady_decoded, error_min_numerical_decoded, error_min_percentage_numerical_decoded, Us_critmin_achieved_numerical_steady_decoded, error_max_theoretical_decoded, error_max_percentage_theoretical_decoded, Us_critmax_achieved_theoretical_steady_decoded, error_max_numerical_decoded, error_max_percentage_numerical_decoded, Us_critmax_achieved_numerical_steady_decoded, error_range_theoretical_decoded, error_range_percentage_theoretical_decoded, error_range_numerical_decoded, error_range_percentage_numerical_decoded ) 


%% Plot the Relative Transmission Network Results.

% Create a plot of the desired membrane voltage output.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Response (Desired)' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Output Neuron Membrane Voltage, U2 [mV]' ), title( 'Relative Transmission Steady State Response (Desired)' )
plot( Us_desired( :, 1 )*( 10^3 ), Us_desired( :, 2 )*( 10^3 ), '-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_desired' ] )

% Create a plot of the achieved numerical membrane voltage output.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Response (Achieved Theoretical)' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Output Neuron Membrane Voltage, U2 [mV]' ), title( 'Relative Transmission Steady State Response (Achieved Theoretical)' )
plot( Us_achieved_theoretical( :, 1 )*( 10^3 ), Us_achieved_theoretical( :, 2 )*( 10^3 ), '-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_achieved_theoretical' ] )

% Create a plot of the achieved numerical membrane voltage output.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Response (Achieved Numerical)' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Output Neuron Membrane Voltage, U2 [mV]' ), title( 'Relative Transmission Steady State Response (Achieved Numerical)' )
plot( Us_achieved_numerical( :, 1 )*( 10^3 ), Us_achieved_numerical( :, 2 )*( 10^3 ), '-', 'Linewidth', 3 )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_achieved_numerical' ] )

% Create a plot of the desired and achieved membrane voltage outputs.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Response (Comparison)' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Output Neuron Membrane Voltage, U2 [mV]' ), title( 'Relative Transmission Steady State Response (Comparison)' )
h1 = plot( Us_desired( :, 1 )*( 10^3 ), Us_desired( :, 2 )*( 10^3 ), '-', 'Linewidth', 3 );
h2 = plot( Us_achieved_theoretical( :, 1 )*( 10^3 ), Us_achieved_theoretical( :, 2 )*( 10^3 ), '-.', 'Linewidth', 3 );
h3 = plot( Us_achieved_numerical( :, 1 )*( 10^3 ), Us_achieved_numerical( :, 2 )*( 10^3 ), '--', 'Linewidth', 3 );
legend( [ h1, h2, h3 ], { 'Desired', 'Achieved (Theoretical)', 'Achieved (Numerical)' }, 'Location', 'Best' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_comparison' ] )

% Create a plot of the desired and achieved membrane voltage outputs.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Decoding (Comparison)' ); hold on, grid on, xlabel( 'Input, x [-]' ), ylabel( 'Output, y [-]' ), title( 'Relative Transmission Steady State Decoding (Comparison)' )
h1 = plot( Us_desired_decoded( :, 1 ), Us_desired_decoded( :, 2 ), '-', 'Linewidth', 3 );
h2 = plot( Us_achieved_theoretical_decoded( :, 1 ), Us_achieved_theoretical_decoded( :, 2 ), '-.', 'Linewidth', 3 );
h3 = plot( Us_achieved_numerical_decoded( :, 1 ), Us_achieved_numerical_decoded( :, 2 ), '--', 'Linewidth', 3 );
legend( [ h1, h2, h3 ], { 'Desired', 'Achieved (Theoretical)', 'Achieved (Numerical)' }, 'Location', 'Best' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_decoding_comparison' ] )

% Create a surface that shows the membrane voltage error.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Error' ); hold on, grid on, xlabel( 'Input Neuron Membrane Voltage, U1 [mV]' ), ylabel( 'Membrane Voltage Error, E [mV]' ), title( 'Relative Transmission Steady State Error' )
plot( Us_achieved_theoretical( :, 1 )*( 10^3 ), errors_theoretical*( 10^3 ), '-', 'Linewidth', 3 )
plot( Us_achieved_numerical( :, 1 )*( 10^3 ), errors_numerical*( 10^3 ), '--', 'Linewidth', 3 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Best', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_error' ] )

% Create a surface that shows the decoding error.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Decoding Error' ); hold on, grid on, xlabel( 'Input, x [-]' ), ylabel( 'Decoding Error, E [-]' ), title( 'Relative Transmission Steady State Decoding Error' )
plot( Us_achieved_theoretical_decoded( :, 1 ), errors_theoretical_decoded, '-', 'Linewidth', 3 )
plot( Us_achieved_numerical_decoded( :, 1 ), errors_numerical_decoded, '--', 'Linewidth', 3 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Best', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_decoding_error' ] )

% Create a surface that shows the decoding error.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Error Percentage' ); hold on, grid on, xlabel( 'Input, x [-]' ), ylabel( 'Membrane Voltage Error Percentage, E [%]' ), title( 'Relative Transmission Steady State Error Percentage' )
plot( Us_achieved_theoretical( :, 1 ), error_percentages_theoretical, '-', 'Linewidth', 3 )
plot( Us_achieved_numerical( :, 1 ), error_percentages_numerical, '--', 'Linewidth', 3 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Best', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_error_percentage' ] )

% Create a surface that shows the decoding error.
fig = figure( 'Color', 'w', 'Name', 'Relative Transmission Steady State Decoding Error Percentage' ); hold on, grid on, xlabel( 'Input, x [-]' ), ylabel( 'Membrane Voltage Decoding Error Percentage, E [%]' ), title( 'Relative Transmission Steady State Decoding Error Percentage' )
plot( Us_achieved_theoretical_decoded( :, 1 ), error_percentages_theoretical_decoded, '-', 'Linewidth', 3 )
plot( Us_achieved_numerical_decoded( :, 1 ), error_percentages_numerical_decoded, '--', 'Linewidth', 3 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Best', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'relative_transmission_ss_response_decoding_error_percentage' ] )


