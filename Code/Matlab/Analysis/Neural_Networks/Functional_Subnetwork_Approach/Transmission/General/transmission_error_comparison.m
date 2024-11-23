%% Transmission Subnetwork Encoding Comparison.

% Clear Everything.
clear, close( 'all' ), clc


%% Define Simulation Parameters.

% Define the save and load directories.
save_directory = '.\Save';                         	% [str] Save Directory.
load_directory = '.\Load';                        	% [str] Load Directory.

% Set a flag to determine whether to simulate.
simulate_flag = true;                             	% [T/F] Simulation Flag. (Determines whether to create a new simulation of the steady state error or to load a previous simulation.)
% simulate_flag = false;                            % [T/F] Simulation Flag. (Determines whether to create a new simulation of the steady state error or to load a previous simulation.)

% Set the level of verbosity.
verbose_flag = true;                            	% [T/F] Printing Flag. (Determines whether to print out information.)

% Define the undetected option.
undetected_option = 'Error';                        % [str] Undetected Option.

% Define the network simulation time step.
network_dt = 1e-3;                                 	% [s] Simulation Time Step.
% network_dt = 1e-4;                             	% [s] Simulation Timestep.

% Define the network simulation duration.
network_tf = 0.5;                                 	% [s] Simulation Duration.
% network_tf = 1;                                 	% [s] Simulation Duration.
% network_tf = 3;                                 	% [s] Simulation Duration.

% Compute the number of simulation timesteps.
n_timesteps = floor( network_tf / network_dt ) + 1; % [#] Number of Simulation Timesteps.

% Construct the simulation times associated with the input currents.
ts = ( 0:network_dt:network_tf )';                 	% [s] Simulation Times.

% Define the integration method.
integration_method = 'RK4';                         % [str] Integration Method (Either FE for Forward Euler or RK4 for Fourth Order Runge-Kutta).


%% Define the Desired Transmission Subnetwork Parameters.

% Create an instance of the network utilities class.
network_utilities = network_utilities_class(  );

% Define the transmission subnetwork parameters.
c = 1.0;                                            % [-] Subnetwork Gain.

% Define the desired mapping operation.
f_desired = @( x ) network_utilities.compute_desired_transmission_sso( x, c );

% Define the domain of the input and output signals.
x_max = 20;
y_max = f_desired( x_max );


%% Define the Encoding & Decoding Operations.

% Define the encoding operation.
f_encode_absolute = @( x ) x*( 10^( -3 ) );
f_encode_relative = @( x, R_encode, R_decode ) ( R_encode./R_decode ).*x;

% Define the decoding operations.
f_decode_absolute = @( U ) U*( 10^3 );
f_decode_relative = @( U, R_encode, R_decode ) ( R_decode./R_encode ).*U;


%% Define Transmission Subnetwork Parameters.

% Define the absolute transmission subnetwork design parameters.
R1_absolute = 20e-3;                                        % [V] Maximum Membrane Voltage (Neuron 1).
Gm1_absolute = 1e-6;                                        % [S] Membrane Conductance (Neuron 1).
Gm2_absolute = 1e-6;                                        % [S] Membrane Conductance (Neuron 2).
Cm1_absolute = 5e-9;                                        % [F] Membrane Capacitance (Neuron 1).
Cm2_absolute = 5e-9;                                        % [F] Membrane Capacitance (Neuron 2).

% Define the relative transmission subnetwork design parameters.
R1_relative = 20e-3;                                     	% [V] Maximum Membrane Voltage (Neuron 1).
R2_relative = 20e-3;                                      	% [V] Maximum Membrane Voltage (Neuron 2).
Gm1_relative = 1e-6;                                       	% [S] Membrane Conductance (Neuron 1).
Gm2_relative = 1e-6;                                      	% [S] Membrane Conductance (Neuron 2).
Cm1_relative = 5e-9;                                       	% [F] Membrane Capacitance (Neuron 1).
Cm2_relative = 5e-9;                                       	% [F] Membrane Capacitance (Neuron 2).

% Store the transmission subnetwork design parameters in a cell.
absolute_transmission_parameters = { c, R1_absolute, Gm1_absolute, Gm2_absolute, Cm1_absolute, Cm2_absolute };
relative_transmission_parameters = { R1_relative, R2_relative, Gm1_relative, Gm2_relative, Cm1_relative, Cm2_relative };


%% Define the Absolute Transmission Subnetwork Input Current Parameters.

% Define the applied current ID.
input_current_ID_absolute = 1;                                  % [#] Absolute Input Current ID.
input_current_ID_relative = 1;                                  % [#] Relative Input Current ID.

% Define the applied current name.
input_current_name_absolute = 'Applied Current 1 (Absolute)';   % [str] Absolute Input Current Name.
input_current_name_relative = 'Applied Current 1 (Relative)';  	% [str] Relative Input Current Name.

% Define the IDs of the neurons to which the currents are applied.
input_current_to_neuron_ID_absolute = 1;                        % [#] Absolute Neuron ID to Which Input Current is Applied.
input_current_to_neuron_ID_relative = 1;                        % [#] Relative Neuron ID to Which Input Current is Applied.

% Define the desired input signal.
xs_desired = x_max*ones( n_timesteps, 1 );

% Encode the input signal.
Us1_desired_absolute = f_encode_absolute( xs_desired );
Us1_desired_relative = f_encode_relative( xs_desired, R1_relative, x_max );

% Define the applied current magnitudes.
Ias1_absolute = Us1_desired_absolute*Gm1_absolute;            	% [A] Applied Current Magnitude.
Ias1_relative = Us1_desired_relative*Gm1_relative;           	% [A] Applied Current Magnitude.


%% Create the Relative Transmission Subnetwork.

% Create an instance of the network class.
network_absolute = network_class( network_dt, network_tf );
network_relative = network_class( network_dt, network_tf );

% Create a transmission subnetwork.
[ c_absolute, Gnas_absolute, R2_absolute, dEs21_absolute, gs21_absolute, Ia2_absolute, neurons_absolute, synapses_absolute, neuron_manager_absolute, synapse_manager_absolute, network_absolute ] = network_absolute.create_transmission_subnetwork( absolute_transmission_parameters, 'absolute', network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, true, true, false, undetected_option );
[ c_relative, Gnas_relative, R2_relative, dEs21_relative, gs21_relative, Ia2_relative, neurons_relative, synapses_relative, neuron_manager_relative, synapse_manager_relative, network_relative ] = network_relative.create_transmission_subnetwork( relative_transmission_parameters, 'relative', network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, true, true, false, undetected_option );

% Create the input applied current.
[ ~, ~, ~, network_absolute.applied_current_manager ] = network_absolute.applied_current_manager.create_applied_current( input_current_ID_absolute, input_current_name_absolute, input_current_to_neuron_ID_absolute, ts, Ias1_absolute, true, network_absolute.applied_current_manager.applied_currents, true, false, network_absolute.applied_current_manager.array_utilities );
[ ~, ~, ~, network_relative.applied_current_manager ] = network_relative.applied_current_manager.create_applied_current( input_current_ID_relative, input_current_name_relative, input_current_to_neuron_ID_relative, ts, Ias1_relative, true, network_relative.applied_current_manager.applied_currents, true, false, network_relative.applied_current_manager.array_utilities );


%% Print Transmission Subnetwork Information.

% Print absolute transmission subnetwork information.
fprintf( '----------------------------------- ABSOLUTE TRANSMISSION SUBNETWORK -----------------------------------\n\n' )
network_absolute.print( network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, verbose_flag );
fprintf( '---------------------------------------------------------------------------------------------------------\n\n\n' )

% Print the relative transmission subnetwork information.
fprintf( '----------------------------------- RELATIVE TRANSMISSION SUBNETWORK -----------------------------------\n\n' )
network_relative.print( network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, verbose_flag );
fprintf( '---------------------------------------------------------------------------------------------------------\n\n\n' )


%% Load the Absolute & Relative Transmission Subnetworks.

% Load the simulation results.
data_absolute = load( [ load_directory, '\', 'absolute_transmission_subnetwork_error' ] );
data_relative = load( [ load_directory, '\', 'relative_transmission_subnetwork_error' ] );

% Unpack the absolute transmission simulation data.
xs_achieved_numerical_absolute = data_absolute.xs_achieved_numerical;
Us_input_absolute = data_absolute.Us_input;
applied_currents_absolute = data_absolute.applied_currents;
Us_achieved_numerical_absolute = data_absolute.Us_achieved_numerical;
ys_achieved_numerical_absolute = data_absolute.ys_achieved_numerical;

% Unpack the absolute transmission simulation data.
xs_achieved_numerical_relative = data_relative.xs_achieved_numerical;
Us_input_relative = data_relative.Us_input;
applied_currents_relative = data_relative.applied_currents;
Us_achieved_numerical_relative = data_relative.Us_achieved_numerical;
ys_achieved_numerical_relative = data_relative.ys_achieved_numerical;


%% Compute the Absolute & Relative Transmission Desired & Achieved (Theory) Network Output.

% Compute the absolute and relative desired subnetwork output.
Us_desired_output_absolute = network_absolute.compute_da_transmission_sso( Us_achieved_numerical_absolute( :, 1 ), c, network_absolute.neuron_manager, undetected_option, network_absolute.network_utilities );
Us_desired_output_relative = network_relative.compute_dr_transmission_sso( Us_achieved_numerical_relative( :, 1 ), c, R1_relative, R2_relative, network_relative.neuron_manager, undetected_option, network_relative.network_utilities );

% Compute the absolute and relative achieved theoretical subnetwork output.
Us_achieved_theoretical_output_absolute = network_absolute.compute_achieved_transmission_sso( Us_achieved_numerical_absolute( :, 1 ), R1_absolute, Gm2_absolute, Ia2_absolute, gs21_absolute, dEs21_absolute, network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, undetected_option, network_absolute.network_utilities );
Us_achieved_theoretical_output_relative = network_relative.compute_achieved_transmission_sso( Us_achieved_numerical_relative( :, 1 ), R1_relative, Gm2_relative, Ia2_relative, gs21_relative, dEs21_relative, network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, undetected_option, network_relative.network_utilities );

% Compute the absolute and relative desired subnetwork states.
Us_desired_absolute = Us_achieved_numerical_absolute; Us_desired_absolute( :, end ) = Us_desired_output_absolute;
Us_desired_relative = Us_achieved_numerical_relative; Us_desired_relative( :, end ) = Us_desired_output_relative;

% Compute the absolute and relative achieved theoretical subnetwork states.
Us_achieved_theoretical_absolute = Us_achieved_numerical_absolute; Us_achieved_theoretical_absolute( :, end ) = Us_achieved_theoretical_output_absolute;
Us_achieved_theoretical_relative = Us_achieved_numerical_relative; Us_achieved_theoretical_relative( :, end ) = Us_achieved_theoretical_output_relative;

% Compute the decoded desired absolute and relative network inputs.
xs_desired_absolute = f_decode_absolute( Us_desired_absolute( :, 1 ) );
xs_desired_relative = f_decode_relative( Us_desired_relative( :, 1 ), R1_relative, x_max );

% Compute the decoded desired absolute and relative network outputs.
ys_desired_absolute = f_decode_absolute( Us_desired_absolute( :, 2 ) );
ys_desired_relative = f_decode_relative( Us_desired_relative( :, 2 ), R2_relative, y_max );

% Compute the decoded achieved theoretical absolute and relative network inputs.
xs_achieved_theoretical_absolute = f_decode_absolute( Us_achieved_theoretical_absolute( :, 1 ) );
xs_achieved_theoretical_relative = f_decode_relative( Us_achieved_theoretical_relative( :, 1 ), R1_relative, x_max );

% Compute the decoded achieved theoretical absolute and relative network outputs.
ys_achieved_theoretical_absolute = f_decode_absolute( Us_achieved_theoretical_absolute( :, 2 ) );
ys_achieved_theoretical_relative = f_decode_relative( Us_achieved_theoretical_relative( :, 2 ), R2_relative, y_max );


%% Compute the Absolute & Relative Transmission Network Error.

% Compute the error between the encoded theoretical output and the desired output.
[ errors_theoretical_encoded_absolute, error_percentages_theoretical_encoded_absolute, error_rmse_theoretical_encoded_absolute, error_rmse_percentage_theoretical_encoded_absolute, error_std_theoretical_encoded_absolute, error_std_percentage_theoretical_encoded_absolute, error_min_theoretical_encoded_absolute, error_min_percentage_theoretical_encoded_absolute, index_min_theoretical_encoded_absolute, error_max_theoretical_encoded_absolute, error_max_percentage_theoretical_encoded_absolute, index_max_theoretical_encoded_absolute, error_range_theoretical_encoded_absolute, error_range_percentage_theoretical_encoded_absolute ] = network_absolute.numerical_method_utilities.compute_error_statistics( Us_achieved_theoretical_absolute, Us_desired_absolute, R2_absolute );
[ errors_theoretical_encoded_relative, error_percentages_theoretical_encoded_relative, error_rmse_theoretical_encoded_relative, error_rmse_percentage_theoretical_encoded_relative, error_std_theoretical_encoded_relative, error_std_percentage_theoretical_encoded_relative, error_min_theoretical_encoded_relative, error_min_percentage_theoretical_encoded_relative, index_min_theoretical_encoded_relative, error_max_theoretical_encoded_relative, error_max_percentage_theoretical_encoded_relative, index_max_theoretical_encoded_relative, error_range_theoretical_encoded_relative, error_range_percentage_theoretical_encoded_relative ] = network_relative.numerical_method_utilities.compute_error_statistics( Us_achieved_theoretical_relative, Us_desired_relative, R2_relative );

% Compute the error between the encoded numerical output and the desired output.
[ errors_numerical_encoded_absolute, error_percentages_numerical_encoded_absolute, error_rmse_numerical_encoded_absolute, error_rmse_percentage_numerical_encoded_absolute, error_std_numerical_encoded_absolute, error_std_percentage_numerical_encoded_absolute, error_min_numerical_encoded_absolute, error_min_percentage_numerical_encoded_absolute, index_min_numerical_encoded_absolute, error_max_numerical_encoded_absolute, error_max_percentage_numerical_encoded_absolute, index_max_numerical_encoded_absolute, error_range_numerical_encoded_absolute, error_range_percentage_numerical_encoded_absolute ] = network_absolute.numerical_method_utilities.compute_error_statistics( Us_achieved_numerical_absolute, Us_desired_absolute, R2_absolute );
[ errors_numerical_encoded_relative, error_percentages_numerical_encoded_relative, error_rmse_numerical_encoded_relative, error_rmse_percentage_numerical_encoded_relative, error_std_numerical_encoded_relative, error_std_percentage_numerical_encoded_relative, error_min_numerical_encoded_relative, error_min_percentage_numerical_encoded_relative, index_min_numerical_encoded_relative, error_max_numerical_encoded_relative, error_max_percentage_numerical_encoded_relative, index_max_numerical_encoded_relative, error_range_numerical_encoded_relative, error_range_percentage_numerical_encoded_relative ] = network_relative.numerical_method_utilities.compute_error_statistics( Us_achieved_numerical_relative, Us_desired_relative, R2_relative );

% Compute the error between the decoded theoretical output and the desired output.
[ errors_theoretical_decoded_absolute, error_percentages_theoretical_decoded_absolute, error_rmse_theoretical_decoded_absolute, error_rmse_percentage_theoretical_decoded_absolute, error_std_theoretical_decoded_absolute, error_std_percentage_theoretical_decoded_absolute, error_min_theoretical_decoded_absolute, error_min_percentage_theoretical_decoded_absolute, index_min_theoretical_decoded_absolute, error_max_theoretical_decoded_absolute, error_max_percentage_theoretical_decoded_absolute, index_max_theoretical_decoded_absolute, error_range_theoretical_decoded_absolute, error_range_percentage_theoretical_decoded_absolute ] = network_absolute.numerical_method_utilities.compute_error_statistics( ys_achieved_theoretical_absolute, ys_desired_absolute, y_max );
[ errors_theoretical_decoded_relative, error_percentages_theoretical_decoded_relative, error_rmse_theoretical_decoded_relative, error_rmse_percentage_theoretical_decoded_relative, error_std_theoretical_decoded_relative, error_std_percentage_theoretical_decoded_relative, error_min_theoretical_decoded_relative, error_min_percentage_theoretical_decoded_relative, index_min_theoretical_decoded_relative, error_max_theoretical_decoded_relative, error_max_percentage_theoretical_decoded_relative, index_max_theoretical_decoded_relative, error_range_theoretical_decoded_relative, error_range_percentage_theoretical_decoded_relative ] = network_relative.numerical_method_utilities.compute_error_statistics( ys_achieved_theoretical_relative, ys_desired_relative, y_max );

% Compute the error between the decoded numerical output and the desired output.
[ errors_numerical_decoded_absolute, error_percentages_numerical_decoded_absolute, error_rmse_numerical_decoded_absolute, error_rmse_percentage_numerical_decoded_absolute, error_std_numerical_decoded_absolute, error_std_percentage_numerical_decoded_absolute, error_min_numerical_decoded_absolute, error_min_percentage_numerical_decoded_absolute, index_min_numerical_decoded_absolute, error_max_numerical_decoded_absolute, error_max_percentage_numerical_decoded_absolute, index_max_numerical_decoded_absolute, error_range_numerical_decoded_absolute, error_range_percentage_numerical_decoded_absolute ] = network_absolute.numerical_method_utilities.compute_error_statistics( ys_achieved_numerical_absolute, ys_desired_absolute, y_max );
[ errors_numerical_decoded_relative, error_percentages_numerical_decoded_relative, error_rmse_numerical_decoded_relative, error_rmse_percentage_numerical_decoded_relative, error_std_numerical_decoded_relative, error_std_percentage_numerical_decoded_relative, error_min_numerical_decoded_relative, error_min_percentage_numerical_decoded_relative, index_min_numerical_decoded_relative, error_max_numerical_decoded_relative, error_max_percentage_numerical_decoded_relative, index_max_numerical_decoded_relative, error_range_numerical_decoded_relative, error_range_percentage_numerical_decoded_relative ] = network_relative.numerical_method_utilities.compute_error_statistics( ys_achieved_numerical_relative, ys_desired_relative, y_max );


%% Print the Absolute & Relative Transmission Summary Statistics.

% Define the absolute header strings.
header_str_encoded_absolute = 'Absolute Transmission Encoded Error Statistics';
header_str_decoded_absolute = 'Absolute Transmission Decoded Error Statistics';

% Define the relative header strings.
header_str_encoded_relative = 'Relative Transmission Encoded Error Statistics';
header_str_decoded_relative = 'Relative Transmission Decoded Error Statistics';

% Define the unit strings.
unit_str_encoded = 'mV';
unit_str_decoded = '-';

% Retrieve the minimum and maximum encoded theoretical and numerical absolute network results.
Us_critmin_achieved_theoretical_steady_absolute = Us_achieved_theoretical_absolute( index_min_theoretical_encoded_absolute, : );
Us_critmin_achieved_numerical_steady_absolute = Us_achieved_numerical_absolute( index_min_numerical_encoded_absolute, : );
Us_critmax_achieved_theoretical_steady_absolute = Us_achieved_theoretical_absolute( index_max_theoretical_encoded_absolute, : );
Us_critmax_achieved_numerical_steady_absolute = Us_achieved_numerical_absolute( index_max_numerical_encoded_absolute, : );

% Retrieve the minimum and maximum encoded theoretical and numerical relative network results.
Us_critmin_achieved_theoretical_steady_relative = Us_achieved_theoretical_relative( index_min_theoretical_encoded_relative, : );
Us_critmin_achieved_numerical_steady_relative = Us_achieved_numerical_relative( index_min_numerical_encoded_relative, : );
Us_critmax_achieved_theoretical_steady_relative = Us_achieved_theoretical_relative( index_max_theoretical_encoded_relative, : );
Us_critmax_achieved_numerical_steady_relative = Us_achieved_numerical_relative( index_max_numerical_encoded_relative, : );

% Retrieve the minimum and maximum decoded theoretical and numerical absolute network results.
ys_critmin_achieved_theoretical_steady_absolute = f_decode_absolute( Us_critmin_achieved_theoretical_steady_absolute );
ys_critmin_achieved_numerical_steady_absolute = f_decode_absolute( Us_critmin_achieved_numerical_steady_absolute );
ys_critmax_achieved_theoretical_steady_absolute = f_decode_absolute( Us_critmax_achieved_theoretical_steady_absolute );
ys_critmax_achieved_numerical_steady_absolute = f_decode_absolute( Us_critmax_achieved_numerical_steady_absolute );

% Retrieve the minimum and maximum decoded theoretical and numerical relative network results.
ys_critmin_achieved_theoretical_steady_relative = f_decode_relative( Us_critmin_achieved_theoretical_steady_relative, [ R1_relative, R2_relative ], [ x_max, y_max ] );
ys_critmin_achieved_numerical_steady_relative = f_decode_relative( Us_critmin_achieved_numerical_steady_relative, [ R1_relative, R2_relative ], [ x_max, y_max ] );
ys_critmax_achieved_theoretical_steady_relative = f_decode_relative( Us_critmax_achieved_theoretical_steady_relative, [ R1_relative, R2_relative ], [ x_max, y_max ] );
ys_critmax_achieved_numerical_steady_relative = f_decode_relative( Us_critmax_achieved_numerical_steady_relative, [ R1_relative, R2_relative ], [ x_max, y_max ] );

% Print the absolute transmission summary statistics.
network_absolute.numerical_method_utilities.print_error_statistics( header_str_encoded_absolute, unit_str_encoded, 10^( -3 ), error_rmse_theoretical_encoded_absolute, error_rmse_percentage_theoretical_encoded_absolute, error_rmse_numerical_encoded_absolute, error_rmse_percentage_numerical_encoded_absolute, error_std_theoretical_encoded_absolute, error_std_percentage_theoretical_encoded_absolute, error_std_numerical_encoded_absolute, error_std_percentage_numerical_encoded_absolute, error_min_theoretical_encoded_absolute, error_min_percentage_theoretical_encoded_absolute, Us_critmin_achieved_theoretical_steady_absolute, error_min_numerical_encoded_absolute, error_min_percentage_numerical_encoded_absolute, Us_critmin_achieved_numerical_steady_absolute, error_max_theoretical_encoded_absolute, error_max_percentage_theoretical_encoded_absolute, Us_critmax_achieved_theoretical_steady_absolute, error_max_numerical_encoded_absolute, error_max_percentage_numerical_encoded_absolute, Us_critmax_achieved_numerical_steady_absolute, error_range_theoretical_encoded_absolute, error_range_percentage_theoretical_encoded_absolute, error_range_numerical_encoded_absolute, error_range_percentage_numerical_encoded_absolute )    
network_absolute.numerical_method_utilities.print_error_statistics( header_str_decoded_absolute, unit_str_decoded, 1, error_rmse_theoretical_decoded_absolute, error_rmse_percentage_theoretical_decoded_absolute, error_rmse_numerical_decoded_absolute, error_rmse_percentage_numerical_decoded_absolute, error_std_theoretical_decoded_absolute, error_std_percentage_theoretical_decoded_absolute, error_std_numerical_decoded_absolute, error_std_percentage_numerical_decoded_absolute, error_min_theoretical_decoded_absolute, error_min_percentage_theoretical_decoded_absolute, ys_critmin_achieved_theoretical_steady_absolute, error_min_numerical_decoded_absolute, error_min_percentage_numerical_decoded_absolute, ys_critmin_achieved_numerical_steady_absolute, error_max_theoretical_decoded_absolute, error_max_percentage_theoretical_decoded_absolute, ys_critmax_achieved_theoretical_steady_absolute, error_max_numerical_decoded_absolute, error_max_percentage_numerical_decoded_absolute, ys_critmax_achieved_numerical_steady_absolute, error_range_theoretical_decoded_absolute, error_range_percentage_theoretical_decoded_absolute, error_range_numerical_decoded_absolute, error_range_percentage_numerical_decoded_absolute )    

% Print the relative transmission summary statistics.
network_relative.numerical_method_utilities.print_error_statistics( header_str_encoded_relative, unit_str_encoded, 10^( -3 ), error_rmse_theoretical_encoded_relative, error_rmse_percentage_theoretical_encoded_relative, error_rmse_numerical_encoded_relative, error_rmse_percentage_numerical_encoded_relative, error_std_theoretical_encoded_relative, error_std_percentage_theoretical_encoded_relative, error_std_numerical_encoded_relative, error_std_percentage_numerical_encoded_relative, error_min_theoretical_encoded_relative, error_min_percentage_theoretical_encoded_relative, Us_critmin_achieved_theoretical_steady_relative, error_min_numerical_encoded_relative, error_min_percentage_numerical_encoded_relative, Us_critmin_achieved_numerical_steady_relative, error_max_theoretical_encoded_relative, error_max_percentage_theoretical_encoded_relative, Us_critmax_achieved_theoretical_steady_relative, error_max_numerical_encoded_relative, error_max_percentage_numerical_encoded_relative, Us_critmax_achieved_numerical_steady_relative, error_range_theoretical_encoded_relative, error_range_percentage_theoretical_encoded_relative, error_range_numerical_encoded_relative, error_range_percentage_numerical_encoded_relative )    
network_relative.numerical_method_utilities.print_error_statistics( header_str_decoded_relative, unit_str_decoded, 1, error_rmse_theoretical_decoded_relative, error_rmse_percentage_theoretical_decoded_relative, error_rmse_numerical_decoded_relative, error_rmse_percentage_numerical_decoded_relative, error_std_theoretical_decoded_relative, error_std_percentage_theoretical_decoded_relative, error_std_numerical_decoded_relative, error_std_percentage_numerical_decoded_relative, error_min_theoretical_decoded_relative, error_min_percentage_theoretical_decoded_relative, ys_critmin_achieved_theoretical_steady_relative, error_min_numerical_decoded_relative, error_min_percentage_numerical_decoded_relative, ys_critmin_achieved_numerical_steady_relative, error_max_theoretical_decoded_relative, error_max_percentage_theoretical_decoded_relative, ys_critmax_achieved_theoretical_steady_relative, error_max_numerical_decoded_relative, error_max_percentage_numerical_decoded_relative, ys_critmax_achieved_numerical_steady_relative, error_range_theoretical_decoded_relative, error_range_percentage_theoretical_decoded_relative, error_range_numerical_decoded_relative, error_range_percentage_numerical_decoded_relative )    


%% Compute the Difference between the Absolute & Relative Transmission Network Errors.


% DEBUGGING THIS SECTION.


% Compute the difference between the theoretical absolute and relative network errors.
[ error_diff_theoretical_encoded, error_percent_diff_theoretical_encoded, error_mse_diff_theoretical_encoded, error_mse_percent_diff_theoretical_encoded, error_std_diff_theoretical_encoded, error_std_percent_diff_theoretical_encoded, error_min_diff_theoretical_encoded, error_min_percent_diff_theoretical_encoded, error_max_diff_theoretical_encoded, error_max_percent_diff_theoretical_encoded ] = compute_error_difference_statistics( errors_theoretical_encoded_absolute, errors_theoretical_encoded_relative, error_percentages_theoretical_encoded_absolute, error_percentages_theoretical_encoded_relative, error_rmse_theoretical_encoded_absolute, error_rmse_theoretical_encoded_relative, error_rmse_percentage_theoretical_encoded_absolute, error_rmse_percentage_theoretical_encoded_relative, error_std_theoretical_encoded_absolute, error_std_theoretical_encoded_relative, error_std_percentage_theoretical_encoded_absolute, error_std_percentage_theoretical_encoded_relative, error_min_theoretical_encoded_absolute, error_min_theoretical_encoded_relative, error_min_percentage_theoretical_encoded_absolute, error_min_percentage_theoretical_encoded_relative, error_max_theoretical_encoded_absolute, error_max_theoretical_encoded_relative, error_max_percentage_theoretical_encoded_absolute, error_max_percentage_theoretical_encoded_relative );
[ error_diff_theoretical_decoded, error_percent_diff_theoretical_decoded, error_mse_diff_theoretical_decoded, error_mse_percent_diff_theoretical_decoded, error_std_diff_theoretical_decoded, error_std_percent_diff_theoretical_decoded, error_min_diff_theoretical_decoded, error_min_percent_diff_theoretical_decoded, error_max_diff_theoretical_decoded, error_max_percent_diff_theoretical_decoded ] = compute_error_difference_statistics( errors_theoretical_decoded_absolute, errors_theoretical_decoded_relative, error_percentages_theoretical_decoded_absolute, error_percentages_theoretical_decoded_relative, error_rmse_theoretical_decoded_absolute, error_rmse_theoretical_decoded_relative, error_rmse_percentage_theoretical_decoded_absolute, error_rmse_percentage_theoretical_decoded_relative, error_std_theoretical_decoded_absolute, error_std_theoretical_decoded_relative, error_std_percentage_theoretical_decoded_absolute, error_std_percentage_theoretical_decoded_relative, error_min_theoretical_decoded_absolute, error_min_theoretical_decoded_relative, error_min_percentage_theoretical_decoded_absolute, error_min_percentage_theoretical_decoded_relative, error_max_theoretical_decoded_absolute, error_max_theoretical_decoded_relative, error_max_percentage_theoretical_decoded_absolute, error_max_percentage_theoretical_decoded_relative );

% Compute the difference between the numerical absolute and relative network errors.
[ error_diff_numerical_encoded, error_percent_diff_numerical_encoded, error_mse_diff_numerical_encoded, error_mse_percent_diff_numerical_encoded, error_std_diff_numerical_encoded, error_std_percent_diff_numerical_encoded, error_min_diff_numerical_encoded, error_min_percent_diff_numerical_encoded, error_max_diff_numerical_encoded, error_max_percent_diff_numerical_encoded ] = compute_error_difference_statistics( errors_numerical_encoded_absolute, errors_numerical_encoded_relative, error_percentages_numerical_encoded_absolute, error_percentages_numerical_encoded_relative, error_rmse_numerical_encoded_absolute, error_rmse_numerical_encoded_relative, error_rmse_percentage_numerical_encoded_absolute, error_rmse_percentage_numerical_encoded_relative, error_std_numerical_encoded_absolute, error_std_numerical_encoded_relative, error_std_percentage_numerical_encoded_absolute, error_std_percentage_numerical_encoded_relative, error_min_numerical_encoded_absolute, error_min_numerical_encoded_relative, error_min_percentage_numerical_encoded_absolute, error_min_percentage_numerical_encoded_relative, error_max_numerical_encoded_absolute, error_max_numerical_encoded_relative, error_max_percentage_numerical_encoded_absolute, error_max_percentage_numerical_encoded_relative );
[ error_diff_numerical_decoded, error_percent_diff_numerical_decoded, error_mse_diff_numerical_decoded, error_mse_percent_diff_numerical_decoded, error_std_diff_numerical_decoded, error_std_percent_diff_numerical_decoded, error_min_diff_numerical_decoded, error_min_percent_diff_numerical_decoded, error_max_diff_numerical_decoded, error_max_percent_diff_numerical_decoded ] = compute_error_difference_statistics( errors_numerical_decoded_absolute, errors_numerical_decoded_relative, error_percentages_numerical_decoded_absolute, error_percentages_numerical_decoded_relative, error_rmse_numerical_decoded_absolute, error_rmse_numerical_decoded_relative, error_rmse_percentage_numerical_decoded_absolute, error_rmse_percentage_numerical_decoded_relative, error_std_numerical_decoded_absolute, error_std_numerical_decoded_relative, error_std_percentage_numerical_decoded_absolute, error_std_percentage_numerical_decoded_relative, error_min_numerical_decoded_absolute, error_min_numerical_decoded_relative, error_min_percentage_numerical_decoded_absolute, error_min_percentage_numerical_decoded_relative, error_max_numerical_decoded_absolute, error_max_numerical_decoded_relative, error_max_percentage_numerical_decoded_absolute, error_max_percentage_numerical_decoded_relative );

% Compute the improvement between the theoretical absolute and relative network errors.
[ error_improv_theoretical_encoded, error_percent_improv_theoretical_encoded, error_mse_improv_theoretical_encoded, error_mse_percent_improv_theoretical_encoded, error_std_improv_theoretical_encoded, error_std_percent_improv_theoretical_encoded, error_min_improv_theoretical_encoded, error_min_percent_improv_theoretical_encoded, error_max_improv_theoretical_encoded, error_max_percent_improv_theoretical_encoded ] = compute_error_improvement_statistics( errors_theoretical_encoded_absolute, errors_theoretical_encoded_relative, error_percentages_theoretical_encoded_absolute, error_percentages_theoretical_encoded_relative, error_rmse_theoretical_encoded_absolute, error_rmse_theoretical_encoded_relative, error_rmse_percentage_theoretical_encoded_absolute, error_rmse_percentage_theoretical_encoded_relative, error_std_theoretical_encoded_absolute, error_std_theoretical_encoded_relative, error_std_percentage_theoretical_encoded_absolute, error_std_percentage_theoretical_encoded_relative, error_min_theoretical_encoded_absolute, error_min_theoretical_encoded_relative, error_min_percentage_theoretical_encoded_absolute, error_min_percentage_theoretical_encoded_relative, error_max_theoretical_encoded_absolute, error_max_theoretical_encoded_relative, error_max_percentage_theoretical_encoded_absolute, error_max_percentage_theoretical_encoded_relative );
[ error_improv_theoretical_decoded, error_percent_improv_theoretical_decoded, error_mse_improv_theoretical_decoded, error_mse_percent_improv_theoretical_decoded, error_std_improv_theoretical_decoded, error_std_percent_improv_theoretical_decoded, error_min_improv_theoretical_decoded, error_min_percent_improv_theoretical_decoded, error_max_improv_theoretical_decoded, error_max_percent_improv_theoretical_decoded ] = compute_error_improvement_statistics( errors_theoretical_decoded_absolute, errors_theoretical_decoded_relative, error_percentages_theoretical_decoded_absolute, error_percentages_theoretical_decoded_relative, error_rmse_theoretical_decoded_absolute, error_rmse_theoretical_decoded_relative, error_rmse_percentage_theoretical_decoded_absolute, error_rmse_percentage_theoretical_decoded_relative, error_std_theoretical_decoded_absolute, error_std_theoretical_decoded_relative, error_std_percentage_theoretical_decoded_absolute, error_std_percentage_theoretical_decoded_relative, error_min_theoretical_decoded_absolute, error_min_theoretical_decoded_relative, error_min_percentage_theoretical_decoded_absolute, error_min_percentage_theoretical_decoded_relative, error_max_theoretical_decoded_absolute, error_max_theoretical_decoded_relative, error_max_percentage_theoretical_decoded_absolute, error_max_percentage_theoretical_decoded_relative );

% Compute the improvement between the numerical absolute and relative network errors.
[ error_improv_numerical_encoded, error_percent_improv_numerical_encoded, error_mse_improv_numerical_encoded, error_mse_percent_improv_numerical_encoded, error_std_improv_numerical_encoded, error_std_percent_improv_numerical_encoded, error_min_improv_numerical_encoded, error_min_percent_improv_numerical_encoded, error_max_improv_numerical_encoded, error_max_percent_improv_numerical_encoded ] = compute_error_improvement_statistics( errors_numerical_encoded_absolute, errors_numerical_encoded_relative, error_percentages_numerical_encoded_absolute, error_percentages_numerical_encoded_relative, error_rmse_numerical_encoded_absolute, error_rmse_numerical_encoded_relative, error_rmse_percentage_numerical_encoded_absolute, error_rmse_percentage_numerical_encoded_relative, error_std_numerical_encoded_absolute, error_std_numerical_encoded_relative, error_std_percentage_numerical_encoded_absolute, error_std_percentage_numerical_encoded_relative, error_min_numerical_encoded_absolute, error_min_numerical_encoded_relative, error_min_percentage_numerical_encoded_absolute, error_min_percentage_numerical_encoded_relative, error_max_numerical_encoded_absolute, error_max_numerical_encoded_relative, error_max_percentage_numerical_encoded_absolute, error_max_percentage_numerical_encoded_relative );
[ error_improv_numerical_decoded, error_percent_improv_numerical_decoded, error_mse_improv_numerical_decoded, error_mse_percent_improv_numerical_decoded, error_std_improv_numerical_decoded, error_std_percent_improv_numerical_decoded, error_min_improv_numerical_decoded, error_min_percent_improv_numerical_decoded, error_max_improv_numerical_decoded, error_max_percent_improv_numerical_decoded ] = compute_error_improvement_statistics( errors_numerical_decoded_absolute, errors_numerical_decoded_relative, error_percentages_numerical_decoded_absolute, error_percentages_numerical_decoded_relative, error_rmse_numerical_decoded_absolute, error_rmse_numerical_decoded_relative, error_rmse_percentage_numerical_decoded_absolute, error_rmse_percentage_numerical_decoded_relative, error_std_numerical_decoded_absolute, error_std_numerical_decoded_relative, error_std_percentage_numerical_decoded_absolute, error_std_percentage_numerical_decoded_relative, error_min_numerical_decoded_absolute, error_min_numerical_decoded_relative, error_min_percentage_numerical_decoded_absolute, error_min_percentage_numerical_decoded_relative, error_max_numerical_decoded_absolute, error_max_numerical_decoded_relative, error_max_percentage_numerical_decoded_absolute, error_max_percentage_numerical_decoded_relative );


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

