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
n_timesteps = floor( network_tf/network_dt ) + 1;   % [#] Number of Simulation Timesteps.

% Construct the simulation times associated with the input currents.
ts = ( 0:network_dt:network_tf )';                 	% [s] Simulation Times.

% Define the integration method.
integration_method = 'RK4';                         % [str] Integration Method (Either FE for Forward Euler or RK4 for Fourth Order Runge-Kutta).


%% Define the Desired Transmission Subnetwork Parameters.

% Create an instance of the network utilities class.
network_utilities = network_utilities_class(  );

% Define the transmission subnetwork parameters.
c = 1.0;            % [-] Subnetwork Gain.

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
R1_relative = 20e-3;                                         % [V] Maximum Membrane Voltage (Neuron 1).
R2_relative = 20e-3;                                         % [V] Maximum Membrane Voltage (Neuron 2).
Gm1_relative = 1e-6;                                         % [S] Membrane Conductance (Neuron 1).
Gm2_relative = 1e-6;                                         % [S] Membrane Conductance (Neuron 2).
Cm1_relative = 5e-9;                                         % [F] Membrane Capacitance (Neuron 1).
Cm2_relative = 5e-9;                                         % [F] Membrane Capacitance (Neuron 2).

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

