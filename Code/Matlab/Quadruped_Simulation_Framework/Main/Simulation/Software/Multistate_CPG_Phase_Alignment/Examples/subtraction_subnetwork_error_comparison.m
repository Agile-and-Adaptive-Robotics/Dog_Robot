% %% Subtraction Subnetwork Encoding Comparison
% 
% % Clear Everything.
% clear, close( 'all' ), clc
% 
% 
% %% Initialize Project Options
% 
% % Define the save and load directories.
% save_directory = '.\Save';
% load_directory = '.\Load';
% 
% % Define the network integration step size.
% network_dt = 1e-3;
% network_tf = 3;
% 
% % Set the subtraction subnetwork properties.
% num_subtraction_neurons = 3;
% c = 1;
% npm_k = [ 1, 1 ];
% s_ks = [ 1, -1 ];
% % s_ks = [ 1, 1 ];
% 
% 
% %% Create Absolute Subtraction Subnetwork.
% 
% % Create an instance of the network class.
% network_absolute = network_class( network_dt, network_tf );
% network_relative = network_class( network_dt, network_tf );
% 
% % Create absolute and relative subtraction subnetwork.
% [ network_absolute, neuron_IDs_absolute, synapse_IDs_absolute, applied_current_IDs_absolute ] = network_absolute.create_absolute_subtraction_subnetwork( num_subtraction_neurons, c, s_ks );
% [ network_relative, neuron_IDs_relative, synapse_IDs_relative, applied_current_IDs_relative ] = network_relative.create_relative_subtraction_subnetwork( num_subtraction_neurons, c, npm_k, s_ks );
% 
% 
% %% Load the Absolute & Relative Subtraction Subnetworks
% 
% % Load the simulation results.
% absolute_subtraction_simulation_data = load( [ load_directory, '\', 'absolute_subtraction_subnetwork_error' ] );
% relative_subtraction_simulation_data = load( [ load_directory, '\', 'relative_subtraction_subnetwork_error' ] );
% 
% % Store the absolute simulation results in separate variables.
% Absolute_Applied_Currents1 = absolute_subtraction_simulation_data.Applied_Currents1;
% Absolute_Applied_Currents2 = absolute_subtraction_simulation_data.Applied_Currents2;
% Us_achieved_absolute = absolute_subtraction_simulation_data.Us_achieved;
% 
% % Store the relative simulation results in separate variables.
% Relative_Applied_Currents1 = relative_subtraction_simulation_data.Applied_Currents1;
% Relative_Applied_Currents2 = relative_subtraction_simulation_data.Applied_Currents2;
% Us_achieved_relative = relative_subtraction_simulation_data.Us_achieved;
% 
% 
% %% Compute the Error in the Steady State Subtraction Subnetwork Responses
% 
% % Get the absolute activation domains of the neurons.
% R3_absolute = network_absolute.neuron_manager.get_neuron_property( 3, 'R' ); R3_absolute = R3_absolute{ 1 };
% 
% % Get the relative activation domains of the neurons.
% R1_relative = network_relative.neuron_manager.get_neuron_property( 1, 'R' ); R1_relative = R1_relative{ 1 };
% R2_relative = network_relative.neuron_manager.get_neuron_property( 2, 'R' ); R2_relative = R2_relative{ 1 };
% R3_relative = network_relative.neuron_manager.get_neuron_property( 3, 'R' ); R3_relative = R3_relative{ 1 };
% 
% % Compute the desired steady state output membrane voltage.
% Us_desired_absolute_output = Us_achieved_absolute( :, :, 1 ) - Us_achieved_absolute( :, :, 2 );
% Us_desired_relative_output = c*R3_relative*( ( 1/npm_k( 1 ) )*( Us_achieved_relative( :, :, 1 )/R1_relative ) - ( 1/npm_k( 2 ) )*( Us_achieved_relative( :, :, 2 )/R2_relative ) );
% 
% % Generate desired steady state membrane voltage matrices.
% Us_desired_absolute = Us_achieved_absolute; Us_desired_absolute( :, :, end ) = Us_desired_absolute_output;
% Us_desired_relative = Us_achieved_relative; Us_desired_relative( :, :, end ) = Us_desired_relative_output;
% 
% % Compute the error between the achieved and desired results.
% error_absolute = Us_achieved_absolute( :, :, end ) - Us_desired_absolute( :, :, end );
% error_relative = Us_achieved_relative( :, :, end ) - Us_desired_relative( :, :, end );
% 
% % Compute the percent error between the achieve and desired results.
% error_absolute_percent = 100*( error_absolute/R3_absolute );
% error_relative_percent = 100*( error_relative/R3_relative );
% 
% % Compute the mean error.
% mse_absolute = ( 1/numel( error_absolute ) )*sqrt( sum( error_absolute.^2, 'all' ) );
% mse_relative = ( 1/numel( error_relative ) )*sqrt( sum( error_relative.^2, 'all' ) );
% 
% % Compute the mean error percentage.
% mse_absolute_percent = 100*( mse_absolute/R3_absolute );
% mse_relative_percent = 100*( mse_relative/R3_relative );
% % mse_absolute_percent = ( 1/numel( error_absolute_percent ) )*sqrt( sum( error_absolute_percent.^2, 'all' ) );
% % mse_relative_percent = ( 1/numel( error_relative_percent ) )*sqrt( sum( error_relative_percent.^2, 'all' ) );
% 
% % Compute the standard deviation of the error.
% std_absolute = std( error_absolute, 0, 'all' );
% std_relative = std( error_relative, 0, 'all' );
% 
% % Compute the standard deviation of the error percentage.
% std_absolute_percent = 100*( std_absolute/R3_absolute );
% std_relative_percent = 100*( std_relative/R3_relative );
% % std_absolute_percent = std( error_absolute_percent, 0, 'all' );
% % std_relative_percent = std( error_relative_percent, 0, 'all' );
% 
% % Compute the maximum errors.
% [ error_absolute_max, index_absolute_max ] = max( abs( error_absolute ), [  ], 'all', 'linear' );
% [ error_relative_max, index_relative_max ] = max( abs( error_relative ), [  ], 'all', 'linear' );
% 
% % Compute the maximum error percentages.
% error_absolute_max_percent = 100*( error_absolute_max/R3_absolute );
% error_relative_max_percent = 100*( error_relative_max/R3_relative );
% % error_absolute_max_percent = max( abs( error_absolute_percent ), [  ], 'all' );
% % error_relative_max_percent = max( abs( error_relative_percent ), [  ], 'all' );
% 
% % Compute the minimum errors.
% [ error_absolute_min, index_absolute_min ] = min( abs( error_absolute ), [  ], 'all', 'linear' );
% [ error_relative_min, index_relative_min ] = min( abs( error_relative ), [  ], 'all', 'linear' );
% 
% % Compute the minimum error percentages.
% error_absolute_min_percent = 100*( error_absolute_min/R3_absolute );
% error_relative_min_percent = 100*( error_relative_min/R3_relative );
% % error_absolute_min_percent = min( abs( error_absolute_percent ), [  ], 'all' );
% % error_relative_min_percent = min( abs( error_relative_percent ), [  ], 'all' );
% 
% % Compute the range of the error.
% error_absolute_range = error_absolute_max - error_absolute_min;
% error_relative_range = error_relative_max - error_relative_min;
% 
% % Compute the range of the error percentages.
% error_absolute_range_percent = 100*( error_absolute_range/R3_absolute );
% error_relative_range_percent = 100*( error_relative_range/R3_relative );
% % error_absolute_range_percent = error_absolute_max_percent - error_absolute_min_percent;
% % error_relative_range_percent = error_relative_max_percent - error_relative_min_percent;
% 
% % Compute the difference in error between the absolute and relative encoding schemes.
% error_difference = abs( error_relative ) - abs( error_absolute );
% error_difference_percent = abs( error_relative_percent ) - abs( error_absolute_percent );
% 
% % Compute the mean squared error difference.
% error_difference_mse = abs( mse_relative ) - abs( mse_absolute );
% error_difference_mse_percent = abs( mse_relative_percent ) - abs( mse_absolute_percent );
% % error_difference_mse = ( 1/numel( error_difference ) )*sqrt( sum( error_difference.^2, 'all' ) );
% % error_difference_mse_percent = ( 1/numel( error_difference_percent ) )*sqrt( sum( error_difference_percent.^2, 'all' ) );
% 
% % Compute the standard deviation difference.
% error_difference_std = abs( std_relative ) - abs( std_absolute );
% error_difference_std_percent = abs( std_relative_percent ) - abs( std_absolute_percent );
% 
% % Compute the maximum error difference.
% error_difference_max = abs( error_relative_max ) - abs( error_absolute_max );
% error_difference_max_percent = abs( error_relative_max_percent ) - abs( error_absolute_max_percent );
% % error_difference_max = max( abs( error_difference ), [  ], 'all' );
% % error_difference_max_percent = max( abs( error_difference_percent ), [  ], 'all' );
% 
% 
% %% Print Out the Summary Information
% 
% % Retrieve the absolute input voltage matrices.
% Us1_achieved_absolute = Us_achieved_absolute( :, :, 1 );
% Us2_achieved_absolute = Us_achieved_absolute( :, :, 2 );
% 
% % Retrieve the relative input voltage matrices.
% Us1_achieved_relative = Us_achieved_relative( :, :, 1 );
% Us2_achieved_relative = Us_achieved_relative( :, :, 2 );
% 
% % Print out the absolute subtraction summary statistics.
% fprintf( 'Absolute Subtraction Summary Statistics\n' )
% fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_absolute, mse_absolute_percent )
% fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_absolute, std_absolute_percent )
% fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_absolute_max, error_absolute_max_percent, Us1_achieved_absolute( index_absolute_max ), Us2_achieved_absolute( index_absolute_max ), 20e-3 )
% fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_absolute_min, error_absolute_min_percent, Us1_achieved_absolute( index_absolute_min ), Us2_achieved_absolute( index_absolute_min ), 20e-3 )
% fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_absolute_range, error_absolute_range_percent )
% 
% fprintf( '\n' )
% fprintf( 'Relative Subtraction Summary Statistics\n' )
% fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_relative, mse_relative_percent )
% fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_relative, std_relative_percent )
% fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_relative_max, error_relative_max_percent, Us1_achieved_relative( index_relative_max ), Us2_achieved_relative( index_relative_max ), 20e-3 )
% fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_relative_min, error_relative_min_percent, Us1_achieved_relative( index_relative_min ), Us2_achieved_relative( index_relative_min ), 20e-3 )
% fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_relative_range, error_relative_range_percent )
% 
% fprintf( '\n' )
% fprintf( 'Absolute vs Relative Subtraction Summary Statistics:\n' )
% fprintf( 'delta MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', error_difference_mse, error_difference_mse_percent )
% fprintf( 'delta STD:\t%9.3e [V] (%6.2f [%%])\n', error_difference_std, error_difference_std_percent )
% fprintf( 'delta Max Error:\t%9.3e [mV] (%6.2f [%%])\n', error_difference_max, error_difference_max_percent )
% 
% 
% %% Plot the Steady State Subtraction Error Surfaces
% 
% % Create a figure that shows the differences between the achieved and desired membrane voltage outputs for the absolute subtraction subnetwork.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [mV]' ), title( 'Absolute Subtraction Subnetwork Steady State Response (Comparison)' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_desired_absolute( :, :, 1 )*(10^3), Us_desired_absolute( :, :, 2 )*(10^3), Us_desired_absolute( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'k', 'FaceAlpha', 0.25 )
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), Us_achieved_absolute( :, :, end )*(10^3), 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
% legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% view( -45, 15 )
% colormap( get_bichromatic_colormap(  ) )
% saveas( fig, [ save_directory, '\', 'Absolute_Subtraction_Subnetwork_Steady_State_Response.png' ] )
% 
% % Create a figure that shows the differences between the achieved and desired membrane voltage outputs for the relative subtraction subnetwork.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [mV]' ), title( 'Relative Subtraction Subnetwork Steady State Response (Comparison)' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_desired_relative( :, :, 1 )*(10^3), Us_desired_relative( :, :, 2 )*(10^3), Us_desired_relative( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'k', 'FaceAlpha', 0.25 )
% surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), Us_achieved_relative( :, :, end )*(10^3), 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
% legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% view( -45, 30 )
% colormap( get_bichromatic_colormap(  ) )
% saveas( fig, [ save_directory, '\', 'Relative_Subtraction_Subnetwork_Steady_State_Response.png' ] )
% 
% % Create a surface that shows the membrane voltage error.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error, E [mV]' ), title( 'Subtraction Subnetwork Steady State Error' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% % surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute*(10^3), 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
% % surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative*(10^3), 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute*(10^3), 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
% surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative*(10^3), 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
% legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% view( 45, 15 )
% saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Comparison.png' ] )
% 
% % Create a surface that shows the membrane voltage error percentage.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error Percentage, E [%]' ), title( 'Subtraction Subnetwork Steady State Error Percentage' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% % surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute_percent, 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
% % surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative_percent, 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute_percent, 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
% surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative_percent, 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
% legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% view( 45, 15 )
% saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Percentage_Comparison.png' ] )
% 
% % Create a surface that shows the difference in error between the absolute and relative subtraction subnetworks.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error Difference, dE [mV]' ), title( 'Subtraction Subnetwork Steady State Error Difference' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_difference*(10^3), 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
% view( 45, 15 )
% colormap( get_bichromatic_colormap(  ) )
% saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Difference.png' ] )
% 
% % Create a surface that shows the difference in error between the absolute and relative percent subtraction subnetworks.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error Difference Percentage, dE [%]' ), title( 'Subtraction Subnetwork Steady State Error Percentage Difference' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_difference_percent, 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
% view( 45, 15 )
% colormap( get_bichromatic_colormap(  ) )
% saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Percentage_Difference.png' ] )





%% Subtraction Subnetwork Encoding Comparison

% Clear Everything.
clear, close( 'all' ), clc


%% Initialize Project Options

% Define the save and load directories.
save_directory = '.\Save';
load_directory = '.\Load';

% Define the network integration step size.
network_dt = 1e-3;
network_tf = 3;

% Set the subtraction subnetwork properties.
num_subtraction_neurons = 4;
c = 1;
npm_k = [ 2, 1 ];
s_ks = [ 1, -1, 1 ];
% s_ks = [ 1, 1 ];


%% Create Absolute Subtraction Subnetwork.

% Create an instance of the network class.
network_absolute = network_class( network_dt, network_tf );
network_relative = network_class( network_dt, network_tf );

% Create absolute and relative subtraction subnetwork.
[ network_absolute, neuron_IDs_absolute, synapse_IDs_absolute, applied_current_IDs_absolute ] = network_absolute.create_absolute_subtraction_subnetwork( num_subtraction_neurons, c, s_ks );
[ network_relative, neuron_IDs_relative, synapse_IDs_relative, applied_current_IDs_relative ] = network_relative.create_relative_subtraction_subnetwork( num_subtraction_neurons, c, npm_k, s_ks );


%% Load the Absolute & Relative Subtraction Subnetworks

% Load the simulation results.
absolute_subtraction_simulation_data = load( [ load_directory, '\', 'absolute_subtraction_subnetwork_error' ] );
relative_subtraction_simulation_data = load( [ load_directory, '\', 'relative_subtraction_subnetwork_error' ] );

% Store the absolute simulation results in separate variables.
Absolute_Applied_Currents1 = absolute_subtraction_simulation_data.Applied_Currents1;
Absolute_Applied_Currents2 = absolute_subtraction_simulation_data.Applied_Currents2;
Us_achieved_absolute = absolute_subtraction_simulation_data.Us_achieved;

% Store the relative simulation results in separate variables.
Relative_Applied_Currents1 = relative_subtraction_simulation_data.Applied_Currents1;
Relative_Applied_Currents2 = relative_subtraction_simulation_data.Applied_Currents2;
Us_achieved_relative = relative_subtraction_simulation_data.Us_achieved;


%% Compute the Error in the Steady State Subtraction Subnetwork Responses

% Get the absolute activation domains of the neurons.
R4_absolute = network_absolute.neuron_manager.get_neuron_property( 4, 'R' ); R4_absolute = R4_absolute{ 1 };

% Get the relative activation domains of the neurons.
R1_relative = network_relative.neuron_manager.get_neuron_property( 1, 'R' ); R1_relative = R1_relative{ 1 };
R2_relative = network_relative.neuron_manager.get_neuron_property( 2, 'R' ); R2_relative = R2_relative{ 1 };
R3_relative = network_relative.neuron_manager.get_neuron_property( 3, 'R' ); R3_relative = R3_relative{ 1 };
R4_relative = network_relative.neuron_manager.get_neuron_property( 4, 'R' ); R4_relative = R4_relative{ 1 };

% Compute the desired steady state output membrane voltage.
Us_desired_absolute_output = Us_achieved_absolute( :, :, 1 ) - Us_achieved_absolute( :, :, 2 ) + Us_achieved_absolute( :, :, 3 );
Us_desired_relative_output = c*R4_relative*( ( 1/npm_k( 1 ) )*( Us_achieved_relative( :, :, 1 )/R1_relative + Us_achieved_relative( :, :, 3 )/R3_relative ) - ( 1/npm_k( 2 ) )*( Us_achieved_relative( :, :, 2 )/R2_relative ) );

% Generate desired steady state membrane voltage matrices.
Us_desired_absolute = Us_achieved_absolute; Us_desired_absolute( :, :, end ) = Us_desired_absolute_output;
Us_desired_relative = Us_achieved_relative; Us_desired_relative( :, :, end ) = Us_desired_relative_output;

% Compute the error between the achieved and desired results.
error_absolute = Us_achieved_absolute( :, :, end ) - Us_desired_absolute( :, :, end );
error_relative = Us_achieved_relative( :, :, end ) - Us_desired_relative( :, :, end );

% Compute the percent error between the achieve and desired results.
error_absolute_percent = 100*( error_absolute/R4_absolute );
error_relative_percent = 100*( error_relative/R4_relative );

% Compute the mean error.
mse_absolute = ( 1/numel( error_absolute ) )*sqrt( sum( error_absolute.^2, 'all' ) );
mse_relative = ( 1/numel( error_relative ) )*sqrt( sum( error_relative.^2, 'all' ) );

% Compute the mean error percentage.
mse_absolute_percent = 100*( mse_absolute/R4_absolute );
mse_relative_percent = 100*( mse_relative/R4_relative );
% mse_absolute_percent = ( 1/numel( error_absolute_percent ) )*sqrt( sum( error_absolute_percent.^2, 'all' ) );
% mse_relative_percent = ( 1/numel( error_relative_percent ) )*sqrt( sum( error_relative_percent.^2, 'all' ) );

% Compute the standard deviation of the error.
std_absolute = std( error_absolute, 0, 'all' );
std_relative = std( error_relative, 0, 'all' );

% Compute the standard deviation of the error percentage.
std_absolute_percent = 100*( std_absolute/R4_absolute );
std_relative_percent = 100*( std_relative/R4_relative );
% std_absolute_percent = std( error_absolute_percent, 0, 'all' );
% std_relative_percent = std( error_relative_percent, 0, 'all' );

% Compute the maximum errors.
[ error_absolute_max, index_absolute_max ] = max( abs( error_absolute ), [  ], 'all', 'linear' );
[ error_relative_max, index_relative_max ] = max( abs( error_relative ), [  ], 'all', 'linear' );

% Compute the maximum error percentages.
error_absolute_max_percent = 100*( error_absolute_max/R4_absolute );
error_relative_max_percent = 100*( error_relative_max/R4_relative );
% error_absolute_max_percent = max( abs( error_absolute_percent ), [  ], 'all' );
% error_relative_max_percent = max( abs( error_relative_percent ), [  ], 'all' );

% Compute the minimum errors.
[ error_absolute_min, index_absolute_min ] = min( abs( error_absolute ), [  ], 'all', 'linear' );
[ error_relative_min, index_relative_min ] = min( abs( error_relative ), [  ], 'all', 'linear' );

% Compute the minimum error percentages.
error_absolute_min_percent = 100*( error_absolute_min/R4_absolute );
error_relative_min_percent = 100*( error_relative_min/R4_relative );
% error_absolute_min_percent = min( abs( error_absolute_percent ), [  ], 'all' );
% error_relative_min_percent = min( abs( error_relative_percent ), [  ], 'all' );

% Compute the range of the error.
error_absolute_range = error_absolute_max - error_absolute_min;
error_relative_range = error_relative_max - error_relative_min;

% Compute the range of the error percentages.
error_absolute_range_percent = 100*( error_absolute_range/R4_absolute );
error_relative_range_percent = 100*( error_relative_range/R4_relative );
% error_absolute_range_percent = error_absolute_max_percent - error_absolute_min_percent;
% error_relative_range_percent = error_relative_max_percent - error_relative_min_percent;

% Compute the difference in error between the absolute and relative encoding schemes.
error_difference = abs( error_relative ) - abs( error_absolute );
error_difference_percent = abs( error_relative_percent ) - abs( error_absolute_percent );

% Compute the mean squared error difference.
error_difference_mse = abs( mse_relative ) - abs( mse_absolute );
error_difference_mse_percent = abs( mse_relative_percent ) - abs( mse_absolute_percent );
% error_difference_mse = ( 1/numel( error_difference ) )*sqrt( sum( error_difference.^2, 'all' ) );
% error_difference_mse_percent = ( 1/numel( error_difference_percent ) )*sqrt( sum( error_difference_percent.^2, 'all' ) );

% Compute the standard deviation difference.
error_difference_std = abs( std_relative ) - abs( std_absolute );
error_difference_std_percent = abs( std_relative_percent ) - abs( std_absolute_percent );

% Compute the maximum error difference.
error_difference_max = abs( error_relative_max ) - abs( error_absolute_max );
error_difference_max_percent = abs( error_relative_max_percent ) - abs( error_absolute_max_percent );
% error_difference_max = max( abs( error_difference ), [  ], 'all' );
% error_difference_max_percent = max( abs( error_difference_percent ), [  ], 'all' );


%% Print Out the Summary Information

% Retrieve the absolute input voltage matrices.
Us1_achieved_absolute = Us_achieved_absolute( :, :, 1 );
Us2_achieved_absolute = Us_achieved_absolute( :, :, 2 );

% Retrieve the relative input voltage matrices.
Us1_achieved_relative = Us_achieved_relative( :, :, 1 );
Us2_achieved_relative = Us_achieved_relative( :, :, 2 );

% Print out the absolute subtraction summary statistics.
fprintf( 'Absolute Subtraction Summary Statistics\n' )
fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_absolute, mse_absolute_percent )
fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_absolute, std_absolute_percent )
fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_absolute_max, error_absolute_max_percent, Us1_achieved_absolute( index_absolute_max ), Us2_achieved_absolute( index_absolute_max ), 20e-3 )
fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_absolute_min, error_absolute_min_percent, Us1_achieved_absolute( index_absolute_min ), Us2_achieved_absolute( index_absolute_min ), 20e-3 )
fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_absolute_range, error_absolute_range_percent )

fprintf( '\n' )
fprintf( 'Relative Subtraction Summary Statistics\n' )
fprintf( 'MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', mse_relative, mse_relative_percent )
fprintf( 'STD: \t\t\t%9.3e [mV] (%6.2f [%%])\n', std_relative, std_relative_percent )
fprintf( 'Max Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_relative_max, error_relative_max_percent, Us1_achieved_relative( index_relative_max ), Us2_achieved_relative( index_relative_max ), 20e-3 )
fprintf( 'Min Error: \t\t%9.3e [mV] (%6.2f [%%]) @ (%9.3e [mV], %9.3e [mV], %9.3e [mV])\n', error_relative_min, error_relative_min_percent, Us1_achieved_relative( index_relative_min ), Us2_achieved_relative( index_relative_min ), 20e-3 )
fprintf( 'Range Error: \t%0.3e [mV] (%6.2f [%%])\n', error_relative_range, error_relative_range_percent )

fprintf( '\n' )
fprintf( 'Absolute vs Relative Subtraction Summary Statistics:\n' )
fprintf( 'delta MSE: \t\t\t%9.3e [mV] (%6.2f [%%])\n', error_difference_mse, error_difference_mse_percent )
fprintf( 'delta STD:\t%9.3e [V] (%6.2f [%%])\n', error_difference_std, error_difference_std_percent )
fprintf( 'delta Max Error:\t%9.3e [mV] (%6.2f [%%])\n', error_difference_max, error_difference_max_percent )


%% Plot the Steady State Subtraction Error Surfaces

% Create a figure that shows the differences between the achieved and desired membrane voltage outputs for the absolute subtraction subnetwork.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [mV]' ), title( 'Absolute Subtraction Subnetwork Steady State Response (Comparison)' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_desired_absolute( :, :, 1 )*(10^3), Us_desired_absolute( :, :, 2 )*(10^3), Us_desired_absolute( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'k', 'FaceAlpha', 0.25 )
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), Us_achieved_absolute( :, :, end )*(10^3), 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
% legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% view( -45, 15 )
% colormap( get_bichromatic_colormap(  ) )
% saveas( fig, [ save_directory, '\', 'Absolute_Subtraction_Subnetwork_Steady_State_Response.png' ] )

% figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [mV]' ), title( 'Absolute Subtraction Subnetwork Steady State Response (Comparison)' )
fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
surf( Us_desired_absolute( :, :, 1 )*(10^3), Us_desired_absolute( :, :, 2 )*(10^3), Us_desired_absolute( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'k', 'FaceAlpha', 0.25 )
surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), Us_achieved_absolute( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
view( -45, 15 )
% colormap( get_bichromatic_colormap(  ) )
saveas( fig, [ save_directory, '\', 'Absolute_Subtraction_Subnetwork_Steady_State_Response.png' ] )

% Create a figure that shows the differences between the achieved and desired membrane voltage outputs for the relative subtraction subnetwork.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [mV]' ), title( 'Relative Subtraction Subnetwork Steady State Response (Comparison)' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_desired_relative( :, :, 1 )*(10^3), Us_desired_relative( :, :, 2 )*(10^3), Us_desired_relative( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'k', 'FaceAlpha', 0.25 )
% surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), Us_achieved_relative( :, :, end )*(10^3), 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
% legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% view( -45, 30 )
% colormap( get_bichromatic_colormap(  ) )
% % colorbar(  )
% saveas( fig, [ save_directory, '\', 'Relative_Subtraction_Subnetwork_Steady_State_Response.png' ] )

% figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage of Output Neuron, U3 [mV]' ), title( 'Relative Subtraction Subnetwork Steady State Response (Comparison)' )
fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
surf( Us_desired_relative( :, :, 1 )*(10^3), Us_desired_relative( :, :, 2 )*(10^3), Us_desired_relative( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'k', 'FaceAlpha', 0.25 )
surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), Us_achieved_relative( :, :, end )*(10^3), 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
legend( { 'Desired', 'Achieved' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
view( -45, 30 )
% colormap( get_bichromatic_colormap(  ) )
% colorbar(  )
saveas( fig, [ save_directory, '\', 'Relative_Subtraction_Subnetwork_Steady_State_Response.png' ] )

% Create a surface that shows the membrane voltage error.
% figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error, E [mV]' ), title( 'Subtraction Subnetwork Steady State Error' )
fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute*(10^3), 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
% surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative*(10^3), 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute*(10^3), 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative*(10^3), 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
view( 45, 15 )
saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Comparison.png' ] )

% Create a surface that shows the membrane voltage error percentage.
% figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error Percentage, E [%]' ), title( 'Subtraction Subnetwork Steady State Error Percentage' )
fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute_percent, 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
% surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative_percent, 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_absolute_percent, 'Edgecolor', 'None', 'Facecolor', 'r', 'FaceAlpha', 0.75 )
surf( Us_achieved_relative( :, :, 1 )*(10^3), Us_achieved_relative( :, :, 2 )*(10^3), error_relative_percent, 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
view( 45, 15 )
saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Percentage_Comparison.png' ] )

% Create a surface that shows the difference in error between the absolute and relative subtraction subnetworks.
% figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error Difference, dE [mV]' ), title( 'Subtraction Subnetwork Steady State Error Difference' )
fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_difference*(10^3), 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
view( 45, 15 )
colormap( get_bichromatic_colormap(  ) )
saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Difference.png' ] )

% Create a surface that shows the difference in error between the absolute and relative percent subtraction subnetworks.
% % figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error Difference Percentage, dE [%]' ), title( 'Subtraction Subnetwork Steady State Error Percentage Difference' )
% fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
% surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_difference_percent, 'Edgecolor', 'Interp', 'Facecolor', 'Interp' )
% view( 45, 15 )
% colormap( get_bichromatic_colormap(  ) )
% saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Percentage_Difference.png' ] )

% figure( 'color', 'w' ), hold on, grid on, rotate3d on, xlabel( 'Membrane Voltage of First Input Neuron, U1 [mV]' ), ylabel( 'Membrane Voltage of Second Input Neuron, U2 [mV]' ), zlabel( 'Membrane Voltage Error Difference Percentage, dE [%]' ), title( 'Subtraction Subnetwork Steady State Error Percentage Difference' )
fig = figure( 'color', 'w' ); hold on, grid on, rotate3d on
surf( Us_achieved_absolute( :, :, 1 )*(10^3), Us_achieved_absolute( :, :, 2 )*(10^3), error_difference_percent, 'Edgecolor', 'None', 'Facecolor', 'b', 'FaceAlpha', 0.75 )
view( 45, 15 )
% colormap( get_bichromatic_colormap(  ) )
saveas( fig, [ save_directory, '\', 'Subtraction_Subnetwork_Approximation_Error_Percentage_Difference.png' ] )
