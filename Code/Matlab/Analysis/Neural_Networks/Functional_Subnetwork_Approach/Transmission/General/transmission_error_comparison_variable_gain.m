%% Transmission Subnetwork Encoding Comparison.

% Clear Everything.
clear, close( 'all' ), clc


%% Define Simulation Parameters.

% Define the save and load directories.
save_directory = '.\Save';                         	% [str] Save Directory.
load_directory = '.\Load';                        	% [str] Load Directory.

% Set a flag to determine whether to simulate.
% simulate_flag = false;                            % [T/F] Simulation Flag. (Determines whether to create a new simulation of the steady state error or to load a previous simulation.)

% Set the level of verbosity.

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

% Define the number of input signals.
n_input_signals = 20;                               % [#] Number of Input Signals.

% Define whether to save simulation data.
simulate_flag = true;                             	% [T/F] Simulation Flag. (Determines whether to create a new simulation of the steady state error or to load a previous simulation.)
save_flag = true;                                   % [T/F] Save Flag.  (Determine whether to save simulation data.
verbose_flag = true;                            	% [T/F] Printing Flag. (Determines whether to print out information.)

% Create an instance of the network utilities class.
network_utilities = network_utilities_class(  );
numerical_method_utilities = numerical_method_utilities_class(  );
plotting_utilities = plotting_utilities_class(  );


%% Define the Desired Transmission Subnetwork Parameters.

% Define the number of gains.
n_gains = 10;

% Define the minimum and maximum gains.
c_min = 1;
c_max = 8;

% Define the transmission subnetwork parameters.
cs = linspace( c_min, c_max, n_gains );                                            % [-] Subnetwork Gain.

% Create arrays to store the encoded steady state output information.
Us_desired_absolute_output = zeros( n_input_signals, n_gains );
Us_theoretical_absolute_output = zeros( n_input_signals, n_gains );
Us_numerical_absolute_output = zeros( n_input_signals, n_gains );

Us_desired_relative_output = zeros( n_input_signals, n_gains );
Us_theoretical_relative_output = zeros( n_input_signals, n_gains );
Us_numerical_relative_output = zeros( n_input_signals, n_gains );

% Create arrays to store the decoded steady state output information.
Xs_desired_absolute_output = zeros( n_input_signals, n_gains );
Xs_theoretical_absolute_output = zeros( n_input_signals, n_gains );
Xs_numerical_absolute_output = zeros( n_input_signals, n_gains );

Xs_desired_relative_output = zeros( n_input_signals, n_gains );
Xs_theoretical_relative_output = zeros( n_input_signals, n_gains );
Xs_numerical_relative_output = zeros( n_input_signals, n_gains );

% Create arrays to store the error information.
errors_theoretical_encoded_absolute = zeros( n_input_signals, n_gains );
errors_percentage_theoretical_encoded_absolute = zeros( n_input_signals, n_gains );
errors_rmse_theoretical_encoded_absolute = zeros( n_gains, 1 ); 
errors_rmse_percentage_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_std_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_std_percentage_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_min_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_min_percentage_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_max_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_max_percentage_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_range_theoretical_encoded_absolute = zeros( n_gains, 1 );
errors_range_percentage_theoretical_encoded_absolute = zeros( n_gains, 1 );

errors_theoretical_encoded_relative = zeros( n_input_signals, n_gains );
errors_percentage_theoretical_encoded_relative = zeros( n_input_signals, n_gains );
errors_rmse_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_rmse_percentage_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_std_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_std_percentage_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_min_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_min_percentage_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_max_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_max_percentage_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_range_theoretical_encoded_relative = zeros( n_gains, 1 );
errors_range_percentage_theoretical_encoded_relative = zeros( n_gains, 1 );

errors_numerical_encoded_absolute = zeros( n_input_signals, n_gains );
errors_percentage_numerical_encoded_absolute = zeros( n_input_signals, n_gains );
errors_rmse_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_rmse_percentage_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_std_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_std_percentage_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_min_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_min_percentage_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_max_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_max_percentage_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_range_numerical_encoded_absolute = zeros( n_gains, 1 );
errors_range_percentage_numerical_encoded_absolute = zeros( n_gains, 1 );

errors_numerical_encoded_relative = zeros( n_input_signals, n_gains );
errors_percentage_numerical_encoded_relative = zeros( n_input_signals, n_gains );
errors_rmse_numerical_encoded_relative = zeros( n_gains, 1 );
errors_rmse_percentage_numerical_encoded_relative = zeros( n_gains, 1 );
errors_std_numerical_encoded_relative = zeros( n_gains, 1 );
errors_std_percentage_numerical_encoded_relative = zeros( n_gains, 1 );
errors_min_numerical_encoded_relative = zeros( n_gains, 1 );
errors_min_percentage_numerical_encoded_relative = zeros( n_gains, 1 );
errors_max_numerical_encoded_relative = zeros( n_gains, 1 );
errors_max_percentage_numerical_encoded_relative = zeros( n_gains, 1 );
errors_range_numerical_encoded_relative = zeros( n_gains, 1 );
errors_range_percentage_numerical_encoded_relative = zeros( n_gains, 1 );

errors_theoretical_decoded_absolute = zeros( n_input_signals, n_gains );
errors_percentage_theoretical_decoded_absolute = zeros( n_input_signals, n_gains );
errors_rmse_theoretical_decoded_absolute = zeros( n_gains, 1 ); 
errors_rmse_percentage_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_std_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_std_percentage_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_min_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_min_percentage_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_max_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_max_percentage_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_range_theoretical_decoded_absolute = zeros( n_gains, 1 );
errors_range_percentage_theoretical_decoded_absolute = zeros( n_gains, 1 );

errors_theoretical_decoded_relative = zeros( n_input_signals, n_gains );
errors_percentage_theoretical_decoded_relative = zeros( n_input_signals, n_gains );
errors_rmse_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_rmse_percentage_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_std_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_std_percentage_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_min_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_min_percentage_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_max_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_max_percentage_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_range_theoretical_decoded_relative = zeros( n_gains, 1 );
errors_range_percentage_theoretical_decoded_relative = zeros( n_gains, 1 );

errors_numerical_decoded_absolute = zeros( n_input_signals, n_gains );
errors_percentage_numerical_decoded_absolute = zeros( n_input_signals, n_gains );
errors_rmse_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_rmse_percentage_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_std_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_std_percentage_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_min_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_min_percentage_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_max_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_max_percentage_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_range_numerical_decoded_absolute = zeros( n_gains, 1 );
errors_range_percentage_numerical_decoded_absolute = zeros( n_gains, 1 );

errors_numerical_decoded_relative = zeros( n_input_signals, n_gains );
errors_percentage_numerical_decoded_relative = zeros( n_input_signals, n_gains );
errors_rmse_numerical_decoded_relative = zeros( n_gains, 1 );
errors_rmse_percentage_numerical_decoded_relative = zeros( n_gains, 1 );
errors_std_numerical_decoded_relative = zeros( n_gains, 1 );
errors_std_percentage_numerical_decoded_relative = zeros( n_gains, 1 );
errors_min_numerical_decoded_relative = zeros( n_gains, 1 );
errors_min_percentage_numerical_decoded_relative = zeros( n_gains, 1 );
errors_max_numerical_decoded_relative = zeros( n_gains, 1 );
errors_max_percentage_numerical_decoded_relative = zeros( n_gains, 1 );
errors_range_numerical_decoded_relative = zeros( n_gains, 1 );
errors_range_percentage_numerical_decoded_relative = zeros( n_gains, 1 );


errors_diff_theoretical_encoded = zeros( n_input_signals, n_gains );
errors_percent_diff_theoretical_encoded = zeros( n_input_signals, n_gains );
errors_mse_diff_theoretical_encoded = zeros( n_gains, 1 );
errors_mse_percent_diff_theoretical_encoded = zeros( n_gains, 1 );
errors_std_diff_theoretical_encoded = zeros( n_gains, 1 );
errors_std_percent_diff_theoretical_encoded = zeros( n_gains, 1 );
errors_min_diff_theoretical_encoded = zeros( n_gains, 1 );
errors_min_percent_diff_theoretical_encoded = zeros( n_gains, 1 );
errors_max_diff_theoretical_encoded = zeros( n_gains, 1 );
errors_max_percent_diff_theoretical_encoded = zeros( n_gains, 1 );

errors_diff_numerical_encoded = zeros( n_input_signals, n_gains );
errors_percent_diff_numerical_encoded = zeros( n_input_signals, n_gains );
errors_mse_diff_numerical_encoded = zeros( n_gains, 1 );
errors_mse_percent_diff_numerical_encoded = zeros( n_gains, 1 );
errors_std_diff_numerical_encoded = zeros( n_gains, 1 );
errors_std_percent_diff_numerical_encoded = zeros( n_gains, 1 );
errors_min_diff_numerical_encoded = zeros( n_gains, 1 );
errors_min_percent_diff_numerical_encoded = zeros( n_gains, 1 );
errors_max_diff_numerical_encoded = zeros( n_gains, 1 );
errors_max_percent_diff_numerical_encoded = zeros( n_gains, 1 );

errors_diff_theoretical_decoded = zeros( n_input_signals, n_gains );
errors_percent_diff_theoretical_decoded = zeros( n_input_signals, n_gains );
errors_mse_diff_theoretical_decoded = zeros( n_gains, 1 );
errors_mse_percent_diff_theoretical_decoded = zeros( n_gains, 1 );
errors_std_diff_theoretical_decoded = zeros( n_gains, 1 );
errors_std_percent_diff_theoretical_decoded = zeros( n_gains, 1 );
errors_min_diff_theoretical_decoded = zeros( n_gains, 1 );
errors_min_percent_diff_theoretical_decoded = zeros( n_gains, 1 );
errors_max_diff_theoretical_decoded = zeros( n_gains, 1 );
errors_max_percent_diff_theoretical_decoded = zeros( n_gains, 1 );

errors_diff_numerical_decoded = zeros( n_input_signals, n_gains );
errors_percent_diff_numerical_decoded = zeros( n_input_signals, n_gains );
errors_mse_diff_numerical_decoded = zeros( n_gains, 1 );
errors_mse_percent_diff_numerical_decoded = zeros( n_gains, 1 );
errors_std_diff_numerical_decoded = zeros( n_gains, 1 );
errors_std_percent_diff_numerical_decoded = zeros( n_gains, 1 );
errors_min_diff_numerical_decoded = zeros( n_gains, 1 );
errors_min_percent_diff_numerical_decoded = zeros( n_gains, 1 );
errors_max_diff_numerical_decoded = zeros( n_gains, 1 );
errors_max_percent_diff_numerical_decoded = zeros( n_gains, 1 );

errors_improv_theoretical_encoded = zeros( n_input_signals, n_gains );
errors_percent_improv_theoretical_encoded = zeros( n_input_signals, n_gains );
errors_mse_improv_theoretical_encoded = zeros( n_gains, 1 );
errors_mse_percent_improv_theoretical_encoded = zeros( n_gains, 1 );
errors_std_improv_theoretical_encoded = zeros( n_gains, 1 );
errors_std_percent_improv_theoretical_encoded = zeros( n_gains, 1 );
errors_min_improv_theoretical_encoded = zeros( n_gains, 1 );
errors_min_percent_improv_theoretical_encoded = zeros( n_gains, 1 );
errors_max_improv_theoretical_encoded = zeros( n_gains, 1 );
errors_max_percent_improv_theoretical_encoded = zeros( n_gains, 1 );

errors_improv_numerical_encoded = zeros( n_input_signals, n_gains );
errors_percent_improv_numerical_encoded = zeros( n_input_signals, n_gains );
errors_mse_improv_numerical_encoded = zeros( n_gains, 1 );
errors_mse_percent_improv_numerical_encoded = zeros( n_gains, 1 );
errors_std_improv_numerical_encoded = zeros( n_gains, 1 );
errors_std_percent_improv_numerical_encoded = zeros( n_gains, 1 );
errors_min_improv_numerical_encoded = zeros( n_gains, 1 );
errors_min_percent_improv_numerical_encoded = zeros( n_gains, 1 );
errors_max_improv_numerical_encoded = zeros( n_gains, 1 );
errors_max_percent_improv_numerical_encoded = zeros( n_gains, 1 );

errors_improv_theoretical_decoded = zeros( n_input_signals, n_gains );
errors_percent_improv_theoretical_decoded = zeros( n_input_signals, n_gains );
errors_mse_improv_theoretical_decoded = zeros( n_gains, 1 );
errors_mse_percent_improv_theoretical_decoded = zeros( n_gains, 1 );
errors_std_improv_theoretical_decoded = zeros( n_gains, 1 );
errors_std_percent_improv_theoretical_decoded = zeros( n_gains, 1 );
errors_min_improv_theoretical_decoded = zeros( n_gains, 1 );
errors_min_percent_improv_theoretical_decoded = zeros( n_gains, 1 );
errors_max_improv_theoretical_decoded = zeros( n_gains, 1 );
errors_max_percent_improv_theoretical_decoded = zeros( n_gains, 1 );

errors_improv_numerical_decoded = zeros( n_input_signals, n_gains );
errors_percent_improv_numerical_decoded = zeros( n_input_signals, n_gains );
errors_mse_improv_numerical_decoded = zeros( n_gains, 1 );
errors_mse_percent_improv_numerical_decoded = zeros( n_gains, 1 );
errors_std_improv_numerical_decoded = zeros( n_gains, 1 );
errors_std_percent_improv_numerical_decoded = zeros( n_gains, 1 );
errors_min_improv_numerical_decoded = zeros( n_gains, 1 );
errors_min_percent_improv_numerical_decoded = zeros( n_gains, 1 );
errors_max_improv_numerical_decoded = zeros( n_gains, 1 );
errors_max_percent_improv_numerical_decoded = zeros( n_gains, 1 );

% Create arrays to store the maximum RK4 step size.
dts_max_absolute = zeros( n_gains, 1 );
dts_max_relative = zeros( n_gains, 1 );

% Create arrays to store the maximum condition numbers.
condition_numbers_max_absolute = zeros( n_gains, 1 );
condition_numbers_max_relative = zeros( n_gains, 1 );

% Create arrays to store the network parameters.
Gnas_absolute = zeros( n_gains, 2 );
Rs2_absolute = zeros( n_gains, 1 );
dEs21_absolute = zeros( n_gains, 1 );
gs21_absolute = zeros( n_gains, 1 );
Ias2_absolute = zeros( n_gains, 1 );

Gnas_relative = zeros( n_gains, 2 );
Rs2_relative = zeros( n_gains, 1 );
dEs21_relative = zeros( n_gains, 1 );
gs21_relative = zeros( n_gains, 1 );
Ias2_relative = zeros( n_gains, 1 );

% Perform the following analysis given each gain value.
for k = 1:n_gains               % Iterate through each of the gains...
    
    % Retrieve the gain.
    c = cs( k );
    
    % Define the desired mapping operation.
    f_desired = @( x, c ) network_utilities.compute_desired_transmission_sso( x, c );
    
    % Define the domain of the input and output signals.
    x_max_input = 20;
    x_max_output = f_desired( x_max_input, c );
    
    
    %% Define the Encoding & Decoding Operations.
    
    % Define the encoding operations.
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
    
    
    %% Define the Absolute & Relative Transmission Subnetwork Input Currents.
    
    % Define the applied current ID.
    input_current_ID_absolute = 1;                                  % [#] Absolute Input Current ID.
    input_current_ID_relative = 1;                                  % [#] Relative Input Current ID.
    
    % Define the applied current name.
    input_current_name_absolute = 'Applied Current 1 (Absolute)';   % [str] Absolute Input Current Name.
    input_current_name_relative = 'Applied Current 1 (Relative)';  	% [str] Relative Input Current Name.
    
    % Define the IDs of the neurons to which the currents are applied.
    input_current_to_neuron_ID_absolute = 1;                        % [#] Absolute Neuron ID to Which Input Current is Applied.
    input_current_to_neuron_ID_relative = 1;                        % [#] Relative Neuron ID to Which Input Current is Applied.
    
    % Define the applied current magnitudes.
    Ias1_absolute = zeros( n_timesteps, 1 );                        % [A] Applied Current Magnitude.
    Ias1_relative = zeros( n_timesteps, 1 );                        % [A] Applied Current Magnitude.
    
    
    %% Create the Relative Transmission Subnetwork.
    
    % Create an instance of the network class.
    network_absolute = network_class( network_dt, network_tf );
    network_relative = network_class( network_dt, network_tf );
    
    % Create a transmission subnetwork.
    [ c_absolute, Gnas_absolute( k, : ), Rs2_absolute( k ), dEs21_absolute( k ), gs21_absolute( k ), Ias2_absolute( k ), neurons_absolute, synapses_absolute, neuron_manager_absolute, synapse_manager_absolute, network_absolute ] = network_absolute.create_transmission_subnetwork( absolute_transmission_parameters, 'absolute', network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, true, true, false, undetected_option );
    [ c_relative, Gnas_relative( k, : ), Rs2_relative( k ), dEs21_relative( k ), gs21_relative( k ), Ias2_relative( k ), neurons_relative, synapses_relative, neuron_manager_relative, synapse_manager_relative, network_relative ] = network_relative.create_transmission_subnetwork( relative_transmission_parameters, 'relative', network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, true, true, false, undetected_option );
    
    % Create the input applied current.
    [ ~, ~, ~, network_absolute.applied_current_manager ] = network_absolute.applied_current_manager.create_applied_current( input_current_ID_absolute, input_current_name_absolute, input_current_to_neuron_ID_absolute, ts, Ias1_absolute, true, network_absolute.applied_current_manager.applied_currents, true, false, network_absolute.applied_current_manager.array_utilities );
    [ ~, ~, ~, network_relative.applied_current_manager ] = network_relative.applied_current_manager.create_applied_current( input_current_ID_relative, input_current_name_relative, input_current_to_neuron_ID_relative, ts, Ias1_relative, true, network_relative.applied_current_manager.applied_currents, true, false, network_relative.applied_current_manager.array_utilities );
    
    
%     %% Print Transmission Subnetwork Information.
%     
%     % Print absolute transmission subnetwork information.
%     fprintf( '----------------------------------- ABSOLUTE TRANSMISSION SUBNETWORK -----------------------------------\n\n' )
%     network_absolute.print( network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, verbose_flag );
%     fprintf( '---------------------------------------------------------------------------------------------------------\n\n\n' )
%     
%     % Print the relative transmission subnetwork information.
%     fprintf( '----------------------------------- RELATIVE TRANSMISSION SUBNETWORK -----------------------------------\n\n' )
%     network_relative.print( network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, verbose_flag );
%     fprintf( '---------------------------------------------------------------------------------------------------------\n\n\n' )
    
    
    %% Simulate the Transmission Network.
    
    % Set additional simulation properties.
    filter_disabled_flag = true;                % [T/F] Filter Disabled Flag.
    set_flag = true;                            % [T/F] Set Flag.
    process_option = 'None';                    % [str] Process Option.
    undetected_option = 'Ignore';               % [str] Undetected Option.
    
    % Determine whether to simulate the network.
    if simulate_flag                            % If we want to simulate the network...
        
        % Define the decoded input signals.
        xs_numerical_input = linspace( 0, x_max_input, n_input_signals )';
        
        % Define the general encoding operations.
        f_general_encode_absolute = @( x, parameters ) f_encode_absolute( x );
        f_general_encode_relative = @( x, parameters ) f_encode_relative( x, parameters{ 1 }, parameters{ 2 } );
        
        % Define the general decoding operations.
        f_general_decode_absolute = @( U, parameters ) f_decode_absolute( U );
        f_general_decode_relative = @( U, parameters ) f_decode_relative( U, parameters{ 1 }, parameters{ 2 } );
        
        % Define the encoding parameters.
        encode_parameters_absolute = {  };
        encode_parameters_relative = { R1_relative, x_max_input };
        
        % Define the decoding parameters.
        decode_parameters_absolute = {  };
        decode_parameters_relative = { R2_relative, x_max_output };
        
        % Compute the decoded steady state simulation results.
        [ xs_numerical_absolute, Us_numerical_absolute, Ias_magnitude_absolute ] = network_absolute.compute_steady_state_simulation_decoded( network_dt, network_tf, integration_method, input_current_ID_absolute, xs_numerical_input, f_general_encode_absolute, encode_parameters_absolute, f_general_decode_absolute, decode_parameters_absolute, network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, network_absolute.applied_voltage_manager, filter_disabled_flag, process_option, undetected_option, network_absolute.network_utilities );
        [ xs_numerical_relative, Us_numerical_relative, Ias_magnitude_relative ] = network_relative.compute_steady_state_simulation_decoded( network_dt, network_tf, integration_method, input_current_ID_absolute, xs_numerical_input, f_general_encode_relative, encode_parameters_relative, f_general_decode_relative, decode_parameters_relative, network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, network_relative.applied_voltage_manager, filter_disabled_flag, process_option, undetected_option, network_relative.network_utilities );
        
        % Determine whether to save the simulation data.
        if save_flag                    % If we want to save the simulation data...
            
            data_absolute.Ias_magnitude = Ias_magnitude_absolute;
            data_absolute.Us_numerical = Us_numerical_absolute;
            data_absolute.xs_numerical = xs_numerical_absolute;
            
            data_relative.Ias_magnitude = Ias_magnitude_relative;
            data_relative.Us_numerical = Us_numerical_relative;
            data_relative.xs_numerical = xs_numerical_relative;
            
            % Define the save file names.
            file_name_absolute = sprintf( 'absolute_transmission_subnetwork_error_gain_%0.0f', c );
            file_name_relative = sprintf( 'relative_transmission_subnetwork_error_gain_%0.0f', c );
            
            % Save the simulation results.
            save( [ save_directory, '\', file_name_absolute ], 'data_absolute' )
            save( [ save_directory, '\', file_name_relative ], 'data_relative' )
            % save( [ save_directory, '\', file_name_absolute ], 'Ias_magnitude_absolute', 'Us_numerical_absolute', 'xs_numerical_absolute' )
            % save( [ save_directory, '\', file_name_relative ], 'Ias_magnitude_relative', 'Us_numerical_relative', 'xs_numerical_relative' )

        end
        
    else                % Otherwise... ( We must want to load data from an existing simulation... )
        
        % Define the load file names.
        file_name_absolute = sprintf( 'absolute_transmission_subnetwork_error_gain_%0.0f', c );
        file_name_relative = sprintf( 'relative_transmission_subnetwork_error_gain_%0.0f', c );
        
        % Load the simulation results.
        data_absolute = load( [ load_directory, '\', file_name_absolute ] );
        data_relative = load( [ load_directory, '\', file_name_relative ] );
        
        % Unpack the steady state simulation data.
        [ xs_numerical_absolute, Us_numerical_absolute, Ias_magnitude_absolute ] = network_absolute.unpack_steady_state_simulation_data( data_absolute.data_absolute );
        [ xs_numerical_relative, Us_numerical_relative, Ias_magnitude_relative ] = network_relative.unpack_steady_state_simulation_data( data_relative.data_relative );
        % [ xs_numerical_absolute, Us_numerical_absolute, Ias_magnitude_absolute ] = network_absolute.unpack_steady_state_simulation_data( data_absolute );
        % [ xs_numerical_relative, Us_numerical_relative, Ias_magnitude_relative ] = network_relative.unpack_steady_state_simulation_data( data_relative );
        
    end
    
    
    %% Compute the Absolute & Relative Transmission Desired & Achieved (Theory) Network Output.
    
    % Initialize the desired decoded steady state response.
    xs_desired_absolute = [ xs_numerical_absolute( :, 1 ), zeros( size( xs_numerical_absolute, 1 ), 1 ) ];
    xs_desired_relative = [ xs_numerical_relative( :, 1 ), zeros( size( xs_numerical_relative, 1 ), 1 ) ];
    
    % Initialize the theoretically achieved decoded steady state response.
    xs_theoretical_absolute = [ xs_numerical_absolute( :, 1 ), zeros( size( xs_numerical_absolute, 1 ), 1 ) ];
    xs_theoretical_relative = [ xs_numerical_relative( :, 1 ), zeros( size( xs_numerical_relative, 1 ), 1 ) ];
    
    % Initialize the desired encoded steady state response.
    Us_desired_absolute = [ Us_numerical_absolute( :, 1 ), zeros( size( Us_numerical_absolute, 1 ), 1 ) ];
    Us_desired_relative = [ Us_numerical_relative( :, 1 ), zeros( size( Us_numerical_relative, 1 ), 1 ) ];
    
    % Initialize the theoretically achieved encoded stady state response.
    Us_theoretical_absolute = [ Us_numerical_absolute( :, 1 ), zeros( size( Us_numerical_absolute, 1 ), 1 ) ];
    Us_theoretical_relative = [ Us_numerical_relative( :, 1 ), zeros( size( Us_numerical_relative, 1 ), 1 ) ];
    
    % Compute the absolute and relative desired subnetwork output.
    Us_desired_absolute( :, 2 ) = network_absolute.compute_da_transmission_sso( Us_desired_absolute( :, 1 ), c, network_absolute.neuron_manager, undetected_option, network_absolute.network_utilities );
    Us_desired_relative( :, 2 ) = network_relative.compute_dr_transmission_sso( Us_desired_relative( :, 1 ), 1.0, R1_relative, R2_relative, network_relative.neuron_manager, undetected_option, network_relative.network_utilities );
    
    % Compute the absolute and relative achieved theoretical subnetwork output.
    Us_theoretical_absolute( :, 2 ) = network_absolute.compute_achieved_transmission_sso( Us_theoretical_absolute( :, 1 ), R1_absolute, Gm2_absolute, Ias2_absolute( k ), gs21_absolute( k ), dEs21_absolute( k ), network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, undetected_option, network_absolute.network_utilities );
    Us_theoretical_relative( :, 2 ) = network_relative.compute_achieved_transmission_sso( Us_theoretical_relative( :, 1 ), R1_relative, Gm2_relative, Ias2_relative( k ), gs21_relative( k ), dEs21_absolute( k ), network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, undetected_option, network_relative.network_utilities );
    
    % Compute the decoded desired absolute and relative network outputs.
    xs_desired_absolute( :, 2 ) = f_decode_absolute( Us_desired_absolute( :, 2 ) );
    xs_desired_relative( :, 2 ) = f_decode_relative( Us_desired_relative( :, 2 ), R2_relative, x_max_output );
    
    % Compute the decoded achieved theoretical absolute and relative network outputs.
    xs_theoretical_absolute( :, 2 ) = f_decode_absolute( Us_theoretical_absolute( :, 2 ) );
    xs_theoretical_relative( :, 2 ) = f_decode_relative( Us_theoretical_relative( :, 2 ), R2_relative, x_max_output );
    
    
    %% Store the Encoded & Decoded Network Outputs.
    
    % Store the encoded network output for plotting.
    Us_desired_absolute_output( :, k ) = Us_desired_absolute( :, 2 );
    Us_theoretical_absolute_output( :, k ) = Us_theoretical_absolute( :, 2 );
    Us_numerical_absolute_output( :, k ) = Us_numerical_absolute( :, 2 );
    
    Us_desired_relative_output( :, k ) = Us_desired_relative( :, 2 );
    Us_theoretical_relative_output( :, k ) = Us_theoretical_relative( :, 2 );
    Us_numerical_relative_output( :, k ) = Us_numerical_relative( :, 2 );
    
    % Store the decoded network output for plotting.
    Xs_desired_absolute_output( :, k ) = xs_desired_absolute( :, 2 );
    Xs_theoretical_absolute_output( :, k ) = xs_theoretical_absolute( :, 2 );
    Xs_numerical_absolute_output( :, k ) = xs_numerical_absolute( :, 2 );
    
    Xs_desired_relative_output( :, k ) = xs_desired_relative( :, 2 );
    Xs_theoretical_relative_output( :, k ) = xs_theoretical_relative( :, 2 );
    Xs_numerical_relative_output( :, k ) = xs_numerical_relative( :, 2 );
    
    
    %% Compute the Absolute & Relative Transmission Network Error.
    
    % Compute the error between the encoded theoretical output and the desired output.
    [ errors_theoretical_encoded_absolute( :, k ), errors_percentage_theoretical_encoded_absolute( :, k ), errors_rmse_theoretical_encoded_absolute( k ), errors_rmse_percentage_theoretical_encoded_absolute( k ), errors_std_theoretical_encoded_absolute( k ), errors_std_percentage_theoretical_encoded_absolute( k ), errors_min_theoretical_encoded_absolute( k ), errors_min_percentage_theoretical_encoded_absolute( k ), index_min_theoretical_encoded_absolute, errors_max_theoretical_encoded_absolute( k ), errors_max_percentage_theoretical_encoded_absolute( k ), index_max_theoretical_encoded_absolute, errors_range_theoretical_encoded_absolute( k ), errors_range_percentage_theoretical_encoded_absolute( k ) ] = numerical_method_utilities.compute_error_statistics( Us_theoretical_absolute, Us_desired_absolute, Rs2_absolute( k ) );
    [ errors_theoretical_encoded_relative( :, k ), errors_percentage_theoretical_encoded_relative( :, k ), errors_rmse_theoretical_encoded_relative( k ), errors_rmse_percentage_theoretical_encoded_relative( k ), errors_std_theoretical_encoded_relative( k ), errors_std_percentage_theoretical_encoded_relative( k ), errors_min_theoretical_encoded_relative( k ), errors_min_percentage_theoretical_encoded_relative( k ), index_min_theoretical_encoded_relative, errors_max_theoretical_encoded_relative( k ), errors_max_percentage_theoretical_encoded_relative( k ), index_max_theoretical_encoded_relative, errors_range_theoretical_encoded_relative( k ), errors_range_percentage_theoretical_encoded_relative( k ) ] = numerical_method_utilities.compute_error_statistics( Us_theoretical_relative, Us_desired_relative, Rs2_relative( k ) );
    
    % Compute the error between the encoded numerical output and the desired output.
    [ errors_numerical_encoded_absolute( :, k ), errors_percentage_numerical_encoded_absolute( :, k ), errors_rmse_numerical_encoded_absolute( k ), errors_rmse_percentage_numerical_encoded_absolute( k ), errors_std_numerical_encoded_absolute( k ), errors_std_percentage_numerical_encoded_absolute( k ), errors_min_numerical_encoded_absolute( k ), errors_min_percentage_numerical_encoded_absolute( k ), index_min_numerical_encoded_absolute, errors_max_numerical_encoded_absolute( k ), errors_max_percentage_numerical_encoded_absolute( k ), index_max_numerical_encoded_absolute, errors_range_numerical_encoded_absolute( k ), errors_range_percentage_numerical_encoded_absolute( k ) ] = numerical_method_utilities.compute_error_statistics( Us_numerical_absolute, Us_desired_absolute, Rs2_absolute( k ) );
    [ errors_numerical_encoded_relative( :, k ), errors_percentage_numerical_encoded_relative( :, k ), errors_rmse_numerical_encoded_relative( k ), errors_rmse_percentage_numerical_encoded_relative( k ), errors_std_numerical_encoded_relative( k ), errors_std_percentage_numerical_encoded_relative( k ), errors_min_numerical_encoded_relative( k ), errors_min_percentage_numerical_encoded_relative( k ), index_min_numerical_encoded_relative, errors_max_numerical_encoded_relative( k ), errors_max_percentage_numerical_encoded_relative( k ), index_max_numerical_encoded_relative, errors_range_numerical_encoded_relative( k ), errors_range_percentage_numerical_encoded_relative( k ) ] = numerical_method_utilities.compute_error_statistics( Us_numerical_relative, Us_desired_relative, Rs2_relative( k ) );
    
    % Compute the error between the decoded theoretical output and the desired output.
    [ errors_theoretical_decoded_absolute( :, k ), errors_percentage_theoretical_decoded_absolute( :, k ), errors_rmse_theoretical_decoded_absolute( k ), errors_rmse_percentage_theoretical_decoded_absolute( k ), errors_std_theoretical_decoded_absolute( k ), errors_std_percentage_theoretical_decoded_absolute( k ), errors_min_theoretical_decoded_absolute( k ), errors_min_percentage_theoretical_decoded_absolute( k ), index_min_theoretical_decoded_absolute, errors_max_theoretical_decoded_absolute( k ), errors_max_percentage_theoretical_decoded_absolute( k ), index_max_theoretical_decoded_absolute, errors_range_theoretical_decoded_absolute( k ), errors_range_percentage_theoretical_decoded_absolute( k ) ] = numerical_method_utilities.compute_error_statistics( xs_theoretical_absolute, xs_desired_absolute, x_max_output );
    [ errors_theoretical_decoded_relative( :, k ), errors_percentage_theoretical_decoded_relative( :, k ), errors_rmse_theoretical_decoded_relative( k ), errors_rmse_percentage_theoretical_decoded_relative( k ), errors_std_theoretical_decoded_relative( k ), errors_std_percentage_theoretical_decoded_relative( k ), errors_min_theoretical_decoded_relative( k ), errors_min_percentage_theoretical_decoded_relative( k ), index_min_theoretical_decoded_relative, errors_max_theoretical_decoded_relative( k ), errors_max_percentage_theoretical_decoded_relative( k ), index_max_theoretical_decoded_relative, errors_range_theoretical_decoded_relative( k ), errors_range_percentage_theoretical_decoded_relative( k ) ] = numerical_method_utilities.compute_error_statistics( xs_theoretical_relative, xs_desired_relative, x_max_output );
    
    % Compute the error between the decoded numerical output and the desired output.
    [ errors_numerical_decoded_absolute( :, k ), errors_percentage_numerical_decoded_absolute( :, k ), errors_rmse_numerical_decoded_absolute( k ), errors_rmse_percentage_numerical_decoded_absolute( k ), errors_std_numerical_decoded_absolute( k ), errors_std_percentage_numerical_decoded_absolute( k ), errors_min_numerical_decoded_absolute( k ), errors_min_percentage_numerical_decoded_absolute( k ), index_min_numerical_decoded_absolute, errors_max_numerical_decoded_absolute( k ), errors_max_percentage_numerical_decoded_absolute( k ), index_max_numerical_decoded_absolute, errors_range_numerical_decoded_absolute( k ), errors_range_percentage_numerical_decoded_absolute( k ) ] = numerical_method_utilities.compute_error_statistics( xs_numerical_absolute, xs_desired_absolute, x_max_output );
    [ errors_numerical_decoded_relative( :, k ), errors_percentage_numerical_decoded_relative( :, k ), errors_rmse_numerical_decoded_relative( k ), errors_rmse_percentage_numerical_decoded_relative( k ), errors_std_numerical_decoded_relative( k ), errors_std_percentage_numerical_decoded_relative( k ), errors_min_numerical_decoded_relative( k ), errors_min_percentage_numerical_decoded_relative( k ), index_min_numerical_decoded_relative, errors_max_numerical_decoded_relative( k ), errors_max_percentage_numerical_decoded_relative( k ), index_max_numerical_decoded_relative, errors_range_numerical_decoded_relative( k ), errors_range_percentage_numerical_decoded_relative( k ) ] = numerical_method_utilities.compute_error_statistics( xs_numerical_relative, xs_desired_relative, x_max_output );
    
    
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
    Us_critmin_theoretical_absolute = Us_theoretical_absolute( index_min_theoretical_encoded_absolute, : );
    Us_critmin_numerical_absolute = Us_numerical_absolute( index_min_numerical_encoded_absolute, : );
    Us_critmax_theoretical_absolute = Us_theoretical_absolute( index_max_theoretical_encoded_absolute, : );
    Us_critmax_numerical_absolute = Us_numerical_absolute( index_max_numerical_encoded_absolute, : );
    
    % Retrieve the minimum and maximum encoded theoretical and numerical relative network results.
    Us_critmin_theoretical_relative = Us_theoretical_relative( index_min_theoretical_encoded_relative, : );
    Us_critmin_numerical_relative = Us_numerical_relative( index_min_numerical_encoded_relative, : );
    Us_critmax_theoretical_relative = Us_theoretical_relative( index_max_theoretical_encoded_relative, : );
    Us_critmax_numerical_relative = Us_numerical_relative( index_max_numerical_encoded_relative, : );
    
    % Retrieve the minimum and maximum decoded theoretical and numerical absolute network results.
    xs_critmin_theoretical_absolute = f_decode_absolute( Us_critmin_theoretical_absolute );
    xs_critmin_numerical_absolute = f_decode_absolute( Us_critmin_numerical_absolute );
    xs_critmax_theoretical_absolute = f_decode_absolute( Us_critmax_theoretical_absolute );
    xs_critmax_numerical_absolute = f_decode_absolute( Us_critmax_numerical_absolute );
    
    % Retrieve the minimum and maximum decoded theoretical and numerical relative network results.
    xs_critmin_theoretical_relative = f_decode_relative( Us_critmin_theoretical_relative, [ R1_relative, R2_relative ], [ x_max_input, x_max_output ] );
    xs_critmin_numerical_relative = f_decode_relative( Us_critmin_numerical_relative, [ R1_relative, R2_relative ], [ x_max_input, x_max_output ] );
    xs_critmax_theoretical_relative = f_decode_relative( Us_critmax_theoretical_relative, [ R1_relative, R2_relative ], [ x_max_input, x_max_output ] );
    xs_critmax_numerical_relative = f_decode_relative( Us_critmax_numerical_relative, [ R1_relative, R2_relative ], [ x_max_input, x_max_output ] );
    
%     % Print the absolute transmission summary statistics.
%     network_absolute.numerical_method_utilities.print_error_statistics( header_str_encoded_absolute, unit_str_encoded, 10^( -3 ), error_rmse_theoretical_encoded_absolute, error_rmse_percentage_theoretical_encoded_absolute, error_rmse_numerical_encoded_absolute, error_rmse_percentage_numerical_encoded_absolute, error_std_theoretical_encoded_absolute, error_std_percentage_theoretical_encoded_absolute, error_std_numerical_encoded_absolute, error_std_percentage_numerical_encoded_absolute, error_min_theoretical_encoded_absolute, error_min_percentage_theoretical_encoded_absolute, Us_critmin_theoretical_absolute, error_min_numerical_encoded_absolute, error_min_percentage_numerical_encoded_absolute, Us_critmin_numerical_absolute, error_max_theoretical_encoded_absolute, error_max_percentage_theoretical_encoded_absolute, Us_critmax_theoretical_absolute, error_max_numerical_encoded_absolute, error_max_percentage_numerical_encoded_absolute, Us_critmax_numerical_absolute, error_range_theoretical_encoded_absolute, error_range_percentage_theoretical_encoded_absolute, error_range_numerical_encoded_absolute, error_range_percentage_numerical_encoded_absolute )
%     network_absolute.numerical_method_utilities.print_error_statistics( header_str_decoded_absolute, unit_str_decoded, 1, error_rmse_theoretical_decoded_absolute, error_rmse_percentage_theoretical_decoded_absolute, error_rmse_numerical_decoded_absolute, error_rmse_percentage_numerical_decoded_absolute, error_std_theoretical_decoded_absolute, error_std_percentage_theoretical_decoded_absolute, error_std_numerical_decoded_absolute, error_std_percentage_numerical_decoded_absolute, error_min_theoretical_decoded_absolute, error_min_percentage_theoretical_decoded_absolute, xs_critmin_theoretical_absolute, error_min_numerical_decoded_absolute, error_min_percentage_numerical_decoded_absolute, xs_critmin_numerical_absolute, error_max_theoretical_decoded_absolute, error_max_percentage_theoretical_decoded_absolute, xs_critmax_theoretical_absolute, error_max_numerical_decoded_absolute, error_max_percentage_numerical_decoded_absolute, xs_critmax_numerical_absolute, error_range_theoretical_decoded_absolute, error_range_percentage_theoretical_decoded_absolute, error_range_numerical_decoded_absolute, error_range_percentage_numerical_decoded_absolute )
%     
%     % Print the relative transmission summary statistics.
%     network_relative.numerical_method_utilities.print_error_statistics( header_str_encoded_relative, unit_str_encoded, 10^( -3 ), error_rmse_theoretical_encoded_relative, error_rmse_percentage_theoretical_encoded_relative, error_rmse_numerical_encoded_relative, error_rmse_percentage_numerical_encoded_relative, error_std_theoretical_encoded_relative, error_std_percentage_theoretical_encoded_relative, error_std_numerical_encoded_relative, error_std_percentage_numerical_encoded_relative, error_min_theoretical_encoded_relative, error_min_percentage_theoretical_encoded_relative, Us_critmin_theoretical_relative, error_min_numerical_encoded_relative, error_min_percentage_numerical_encoded_relative, Us_critmin_numerical_relative, error_max_theoretical_encoded_relative, error_max_percentage_theoretical_encoded_relative, Us_critmax_theoretical_relative, error_max_numerical_encoded_relative, error_max_percentage_numerical_encoded_relative, Us_critmax_numerical_relative, error_range_theoretical_encoded_relative, error_range_percentage_theoretical_encoded_relative, error_range_numerical_encoded_relative, error_range_percentage_numerical_encoded_relative )
%     network_relative.numerical_method_utilities.print_error_statistics( header_str_decoded_relative, unit_str_decoded, 1, error_rmse_theoretical_decoded_relative, error_rmse_percentage_theoretical_decoded_relative, error_rmse_numerical_decoded_relative, error_rmse_percentage_numerical_decoded_relative, error_std_theoretical_decoded_relative, error_std_percentage_theoretical_decoded_relative, error_std_numerical_decoded_relative, error_std_percentage_numerical_decoded_relative, error_min_theoretical_decoded_relative, error_min_percentage_theoretical_decoded_relative, xs_critmin_theoretical_relative, error_min_numerical_decoded_relative, error_min_percentage_numerical_decoded_relative, xs_critmin_numerical_relative, error_max_theoretical_decoded_relative, error_max_percentage_theoretical_decoded_relative, xs_critmax_theoretical_relative, error_max_numerical_decoded_relative, error_max_percentage_numerical_decoded_relative, xs_critmax_numerical_relative, error_range_theoretical_decoded_relative, error_range_percentage_theoretical_decoded_relative, error_range_numerical_decoded_relative, error_range_percentage_numerical_decoded_relative )
    
    
    %% Compute the Difference between the Absolute & Relative Transmission Network Errors.
    
    % Compute the difference between the theoretical absolute and relative network errors.
    [ errors_diff_theoretical_encoded( :, k ), errors_percent_diff_theoretical_encoded( :, k ), errors_mse_diff_theoretical_encoded( k ), errors_mse_percent_diff_theoretical_encoded( k ), errors_std_diff_theoretical_encoded( k ), errors_std_percent_diff_theoretical_encoded( k ), errors_min_diff_theoretical_encoded( k ), errors_min_percent_diff_theoretical_encoded( k ), errors_max_diff_theoretical_encoded( k ), errors_max_percent_diff_theoretical_encoded( k ) ] = numerical_method_utilities.compute_error_difference_statistics( errors_theoretical_encoded_absolute( :, k ), errors_theoretical_encoded_relative( :, k ), errors_percentage_theoretical_encoded_absolute( :, k ), errors_percentage_theoretical_encoded_relative( :, k ), errors_rmse_theoretical_encoded_absolute( k ), errors_rmse_theoretical_encoded_relative( k ), errors_rmse_percentage_theoretical_encoded_absolute( k ), errors_rmse_percentage_theoretical_encoded_relative( k ), errors_std_theoretical_encoded_absolute( k ), errors_std_theoretical_encoded_relative( k ), errors_std_percentage_theoretical_encoded_absolute( k ), errors_std_percentage_theoretical_encoded_relative( k ), errors_min_theoretical_encoded_absolute( k ), errors_min_theoretical_encoded_relative( k ), errors_min_percentage_theoretical_encoded_absolute( k ), errors_min_percentage_theoretical_encoded_relative( k ), errors_max_theoretical_encoded_absolute( k ), errors_max_theoretical_encoded_relative( k ), errors_max_percentage_theoretical_encoded_absolute( k ), errors_max_percentage_theoretical_encoded_relative( k ) );
    [ errors_diff_theoretical_decoded( :, k ), errors_percent_diff_theoretical_decoded( :, k ), errors_mse_diff_theoretical_decoded( k ), errors_mse_percent_diff_theoretical_decoded( k ), errors_std_diff_theoretical_decoded( k ), errors_std_percent_diff_theoretical_decoded( k ), errors_min_diff_theoretical_decoded( k ), errors_min_percent_diff_theoretical_decoded( k ), errors_max_diff_theoretical_decoded( k ), errors_max_percent_diff_theoretical_decoded( k ) ] = numerical_method_utilities.compute_error_difference_statistics( errors_theoretical_decoded_absolute( :, k ), errors_theoretical_decoded_relative( :, k ), errors_percentage_theoretical_decoded_absolute( :, k ), errors_percentage_theoretical_decoded_relative( :, k ), errors_rmse_theoretical_decoded_absolute( k ), errors_rmse_theoretical_decoded_relative( k ), errors_rmse_percentage_theoretical_decoded_absolute( k ), errors_rmse_percentage_theoretical_decoded_relative( k ), errors_std_theoretical_decoded_absolute( k ), errors_std_theoretical_decoded_relative( k ), errors_std_percentage_theoretical_decoded_absolute( k ), errors_std_percentage_theoretical_decoded_relative( k ), errors_min_theoretical_decoded_absolute( k ), errors_min_theoretical_decoded_relative( k ), errors_min_percentage_theoretical_decoded_absolute( k ), errors_min_percentage_theoretical_decoded_relative( k ), errors_max_theoretical_decoded_absolute( k ), errors_max_theoretical_decoded_relative( k ), errors_max_percentage_theoretical_decoded_absolute( k ), errors_max_percentage_theoretical_decoded_relative( k ) );
    
    % Compute the difference between the numerical absolute and relative network errors.
    [ errors_diff_numerical_encoded( :, k ), errors_percent_diff_numerical_encoded( :, k ), errors_mse_diff_numerical_encoded( k ), errors_mse_percent_diff_numerical_encoded( k ), errors_std_diff_numerical_encoded( k ), errors_std_percent_diff_numerical_encoded( k ), errors_min_diff_numerical_encoded( k ), errors_min_percent_diff_numerical_encoded( k ), errors_max_diff_numerical_encoded( k ), errors_max_percent_diff_numerical_encoded( k ) ] = numerical_method_utilities.compute_error_difference_statistics( errors_numerical_encoded_absolute( :, k ), errors_numerical_encoded_relative( :, k ), errors_percentage_numerical_encoded_absolute( :, k ), errors_percentage_numerical_encoded_relative( :, k ), errors_rmse_numerical_encoded_absolute( k ), errors_rmse_numerical_encoded_relative( k ), errors_rmse_percentage_numerical_encoded_absolute( k ), errors_rmse_percentage_numerical_encoded_relative( k ), errors_std_numerical_encoded_absolute( k ), errors_std_numerical_encoded_relative( k ), errors_std_percentage_numerical_encoded_absolute( k ), errors_std_percentage_numerical_encoded_relative( k ), errors_min_numerical_encoded_absolute( k ), errors_min_numerical_encoded_relative( k ), errors_min_percentage_numerical_encoded_absolute( k ), errors_min_percentage_numerical_encoded_relative( k ), errors_max_numerical_encoded_absolute( k ), errors_max_numerical_encoded_relative( k ), errors_max_percentage_numerical_encoded_absolute( k ), errors_max_percentage_numerical_encoded_relative( k ) );
    [ errors_diff_numerical_decoded( :, k ), errors_percent_diff_numerical_decoded( :, k ), errors_mse_diff_numerical_decoded( k ), errors_mse_percent_diff_numerical_decoded( k ), errors_std_diff_numerical_decoded( k ), errors_std_percent_diff_numerical_decoded( k ), errors_min_diff_numerical_decoded( k ), errors_min_percent_diff_numerical_decoded( k ), errors_max_diff_numerical_decoded( k ), errors_max_percent_diff_numerical_decoded( k ) ] = numerical_method_utilities.compute_error_difference_statistics( errors_numerical_decoded_absolute( :, k ), errors_numerical_decoded_relative( :, k ), errors_percentage_numerical_decoded_absolute( :, k ), errors_percentage_numerical_decoded_relative( :, k ), errors_rmse_numerical_decoded_absolute( k ), errors_rmse_numerical_decoded_relative( k ), errors_rmse_percentage_numerical_decoded_absolute( k ), errors_rmse_percentage_numerical_decoded_relative( k ), errors_std_numerical_decoded_absolute( k ), errors_std_numerical_decoded_relative( k ), errors_std_percentage_numerical_decoded_absolute( k ), errors_std_percentage_numerical_decoded_relative( k ), errors_min_numerical_decoded_absolute( k ), errors_min_numerical_decoded_relative( k ), errors_min_percentage_numerical_decoded_absolute( k ), errors_min_percentage_numerical_decoded_relative( k ), errors_max_numerical_decoded_absolute( k ), errors_max_numerical_decoded_relative( k ), errors_max_percentage_numerical_decoded_absolute( k ), errors_max_percentage_numerical_decoded_relative( k ) );
    
    % Compute the improvement between the theoretical absolute and relative network errors.
    [ errors_improv_theoretical_encoded( :, k ), errors_percent_improv_theoretical_encoded( :, k ), errors_mse_improv_theoretical_encoded( k ), errors_mse_percent_improv_theoretical_encoded( k ), errors_std_improv_theoretical_encoded( k ), errors_std_percent_improv_theoretical_encoded( k ), errors_min_improv_theoretical_encoded( k ), errors_min_percent_improv_theoretical_encoded( k ), errors_max_improv_theoretical_encoded( k ), errors_max_percent_improv_theoretical_encoded( k ) ] = numerical_method_utilities.compute_error_improvement_statistics( errors_theoretical_encoded_absolute( :, k ), errors_theoretical_encoded_relative( :, k ), errors_percentage_theoretical_encoded_absolute( :, k ), errors_percentage_theoretical_encoded_relative( :, k ), errors_rmse_theoretical_encoded_absolute( k ), errors_rmse_theoretical_encoded_relative( k ), errors_rmse_percentage_theoretical_encoded_absolute( k ), errors_rmse_percentage_theoretical_encoded_relative( k ), errors_std_theoretical_encoded_absolute( k ), errors_std_theoretical_encoded_relative( k ), errors_std_percentage_theoretical_encoded_absolute( k ), errors_std_percentage_theoretical_encoded_relative( k ), errors_min_theoretical_encoded_absolute( k ), errors_min_theoretical_encoded_relative( k ), errors_min_percentage_theoretical_encoded_absolute( k ), errors_min_percentage_theoretical_encoded_relative( k ), errors_max_theoretical_encoded_absolute( k ), errors_max_theoretical_encoded_relative( k ), errors_max_percentage_theoretical_encoded_absolute( k ), errors_max_percentage_theoretical_encoded_relative( k ) );
    [ errors_improv_theoretical_decoded( :, k ), errors_percent_improv_theoretical_decoded( :, k ), errors_mse_improv_theoretical_decoded( k ), errors_mse_percent_improv_theoretical_decoded( k ), errors_std_improv_theoretical_decoded( k ), errors_std_percent_improv_theoretical_decoded( k ), errors_min_improv_theoretical_decoded( k ), errors_min_percent_improv_theoretical_decoded( k ), errors_max_improv_theoretical_decoded( k ), errors_max_percent_improv_theoretical_decoded( k ) ] = numerical_method_utilities.compute_error_improvement_statistics( errors_theoretical_decoded_absolute( :, k ), errors_theoretical_decoded_relative( :, k ), errors_percentage_theoretical_decoded_absolute( :, k ), errors_percentage_theoretical_decoded_relative( :, k ), errors_rmse_theoretical_decoded_absolute( k ), errors_rmse_theoretical_decoded_relative( k ), errors_rmse_percentage_theoretical_decoded_absolute( k ), errors_rmse_percentage_theoretical_decoded_relative( k ), errors_std_theoretical_decoded_absolute( k ), errors_std_theoretical_decoded_relative( k ), errors_std_percentage_theoretical_decoded_absolute( k ), errors_std_percentage_theoretical_decoded_relative( k ), errors_min_theoretical_decoded_absolute( k ), errors_min_theoretical_decoded_relative( k ), errors_min_percentage_theoretical_decoded_absolute( k ), errors_min_percentage_theoretical_decoded_relative( k ), errors_max_theoretical_decoded_absolute( k ), errors_max_theoretical_decoded_relative( k ), errors_max_percentage_theoretical_decoded_absolute( k ), errors_max_percentage_theoretical_decoded_relative( k ) );
    
    % Compute the improvement between the numerical absolute and relative network errors.
    [ errors_improv_numerical_encoded( :, k ), errors_percent_improv_numerical_encoded( :, k ), errors_mse_improv_numerical_encoded( k ), errors_mse_percent_improv_numerical_encoded( k ), errors_std_improv_numerical_encoded( k ), errors_std_percent_improv_numerical_encoded( k ), errors_min_improv_numerical_encoded( k ), errors_min_percent_improv_numerical_encoded( k ), errors_max_improv_numerical_encoded( k ), errors_max_percent_improv_numerical_encoded( k ) ] = numerical_method_utilities.compute_error_improvement_statistics( errors_numerical_encoded_absolute( :, k ), errors_numerical_encoded_relative( :, k ), errors_percentage_numerical_encoded_absolute( :, k ), errors_percentage_numerical_encoded_relative( :, k ), errors_rmse_numerical_encoded_absolute( k ), errors_rmse_numerical_encoded_relative( k ), errors_rmse_percentage_numerical_encoded_absolute( k ), errors_rmse_percentage_numerical_encoded_relative( k ), errors_std_numerical_encoded_absolute( k ), errors_std_numerical_encoded_relative( k ), errors_std_percentage_numerical_encoded_absolute( k ), errors_std_percentage_numerical_encoded_relative( k ), errors_min_numerical_encoded_absolute( k ), errors_min_numerical_encoded_relative( k ), errors_min_percentage_numerical_encoded_absolute( k ), errors_min_percentage_numerical_encoded_relative( k ), errors_max_numerical_encoded_absolute( k ), errors_max_numerical_encoded_relative( k ), errors_max_percentage_numerical_encoded_absolute( k ), errors_max_percentage_numerical_encoded_relative( k ) );
    [ errors_improv_numerical_decoded( :, k ), errors_percent_improv_numerical_decoded( :, k ), errors_mse_improv_numerical_decoded( k ), errors_mse_percent_improv_numerical_decoded( k ), errors_std_improv_numerical_decoded( k ), errors_std_percent_improv_numerical_decoded( k ), errors_min_improv_numerical_decoded( k ), errors_min_percent_improv_numerical_decoded( k ), errors_max_improv_numerical_decoded( k ), errors_max_percent_improv_numerical_decoded( k ) ] = numerical_method_utilities.compute_error_improvement_statistics( errors_numerical_decoded_absolute( :, k ), errors_numerical_decoded_relative( :, k ), errors_percentage_numerical_decoded_absolute( :, k ), errors_percentage_numerical_decoded_relative( :, k ), errors_rmse_numerical_decoded_absolute( k ), errors_rmse_numerical_decoded_relative( k ), errors_rmse_percentage_numerical_decoded_absolute( k ), errors_rmse_percentage_numerical_decoded_relative( k ), errors_std_numerical_decoded_absolute( k ), errors_std_numerical_decoded_relative( k ), errors_std_percentage_numerical_decoded_absolute( k ), errors_std_percentage_numerical_decoded_relative( k ), errors_min_numerical_decoded_absolute( k ), errors_min_numerical_decoded_relative( k ), errors_min_percentage_numerical_decoded_absolute( k ), errors_min_percentage_numerical_decoded_relative( k ), errors_max_numerical_decoded_absolute( k ), errors_max_numerical_decoded_relative( k ), errors_max_percentage_numerical_decoded_absolute( k ), errors_max_percentage_numerical_decoded_relative( k ) );
    
    
    %% Compute the Transmission Subnetwork Numerical Stability Information.
    
    % Define the property retrieval settings.
    as_matrix_flag = true;
    
    % Define the stability analysis timestep seed.
    dt0 = 1e-6;                                                                                                                                                             % [s] Numerical Stability Time Step.
    
    % Retrieve the properties necessary to compute the numerical stability parameters for an absolute and relative transmission subnetwork.
    [ Cms_absolute, Gms_absolute, Rs_absolute, gs_absolute, dEs_absolute, Ias_absolute ] = network_absolute.get_numerical_stability_parameters( network_absolute.neuron_manager, network_absolute.synapse_manager, as_matrix_flag, undetected_option );
    [ Cms_relative, Gms_relative, Rs_relative, gs_relative, dEs_relative, Ias_relative ] = network_relative.get_numerical_stability_parameters( network_relative.neuron_manager, network_relative.synapse_manager, as_matrix_flag, undetected_option );
    
    % Compute the realtive transmission steady state output.
    [ ~, As_absolute, dts_absolute, condition_numbers_absolute ] = network_absolute.achieved_transmission_RK4_stability_analysis( Us_desired_absolute( :, 1 ), Cms_absolute, Gms_absolute, Rs_absolute, Ias_absolute, gs_absolute, dEs_absolute, dt0, network_absolute.neuron_manager, network_absolute.synapse_manager, network_absolute.applied_current_manager, undetected_option, network_absolute.network_utilities );
    [ ~, As_relative, dts_relative, condition_numbers_relative ] = network_relative.achieved_transmission_RK4_stability_analysis( Us_desired_relative( :, 1 ), Cms_relative, Gms_relative, Rs_relative, Ias_relative, gs_relative, dEs_relative, dt0, network_relative.neuron_manager, network_relative.synapse_manager, network_relative.applied_current_manager, undetected_option, network_relative.network_utilities );
    
    % Retrieve the maximum RK4 step size.
    [ dts_max_absolute( k ), indexes_dt_absolute ] = max( dts_absolute );
    [ dts_max_relative( k ), indexes_dt_relative ] = max( dts_relative );
    
    % Retrieve the maximum condition number.
    [ condition_numbers_max_absolute( k ), indexes_condition_number_absolute ] = max( condition_numbers_absolute );
    [ condition_numbers_max_relative( k ), indexes_condition_number_relative ] = max( condition_numbers_relative );
    
    
    %% Print the Numerical Stability Information.
    
    % % Print out the stability information.
    % network_absolute.numerical_method_utilities.print_numerical_stability_info( As_absolute, dts_absolute, network_dt, condition_numbers_absolute );
    % network_relative.numerical_method_utilities.print_numerical_stability_info( As_relative, dts_relative, network_dt, condition_numbers_relative );
    
    
end

%% Define Plotting Parameters.

% Define a scaling factor.
scale = 1e3;

% Define the line colors.
color1 = [ 0.0000, 0.4470, 0.7410, 1.0000 ];
color2 = [ 0.8500, 0.3250, 0.0980, 1.0000 ];

% Retrieve the numerical input.
xs_numerical_input = xs_numerical_absolute( :, 1 );
Us_numerical_input = Us_numerical_absolute( :, 1 );

% Define the input grid.
[ Cs, Xs_input ] = meshgrid( cs, xs_numerical_input );
[ ~, Us_input ] = meshgrid( cs, Us_numerical_input );


%% Plot the Encoded Steady State Behavior.

% Plot the encoded steady state behavior.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Absolute Encoded Steady State Response' ); hold on, grid on, rotate3d on, view( 45, 30 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, U1 [mV]' ), zlabel( 'Encoded Output, U2 [mV]' ), title( 'Transmission: Absolute Encoded Steady State Error' )
surf( Cs, scale*Us_input, scale*Us_desired_absolute_output, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, scale*Us_input, scale*Us_theoretical_absolute_output, 'Edgecolor', 'None', 'Facecolor', 'g', 'Facealpha', 0.5 )
surf( Cs, scale*Us_input, scale*Us_numerical_absolute_output, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Desired', 'Achieved (Theory)', 'Achieved (Numerical)' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_absolute_encoded_ss_response_gain' ] ) 

fig = figure( 'Color', 'w', 'Name', 'Transmission: Relative Encoded Steady State Response' ); hold on, grid on, rotate3d on, view( 45, 30 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, U1 [mV]' ), zlabel( 'Encoded Output, U2 [mV]' ), title( 'Transmission: Relative Encoded Steady State Error' )
surf( Cs, scale*Us_input, scale*Us_desired_relative_output, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, scale*Us_input, scale*Us_theoretical_relative_output, 'Edgecolor', 'None', 'Facecolor', 'g', 'Facealpha', 0.5 )
surf( Cs, scale*Us_input, scale*Us_numerical_relative_output, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Desired', 'Achieved (Theory)', 'Achieved (Numerical)' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_relative_encoded_ss_response_gain' ] ) 


%% Plot the Decoded Steady State Behavior.

% Plot the absolute decoded steady state behavior.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Absolute Decoded Steady State Response' ); hold on, grid on, rotate3d on, view( 45, 30 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Output, x2 [-]' ), title( 'Transmission: Absolute Decoded Steady State Error' )
surf( Cs, Xs_input, Xs_desired_absolute_output, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Xs_input, Xs_theoretical_absolute_output, 'Edgecolor', 'None', 'Facecolor', 'g', 'Facealpha', 0.5 )
surf( Cs, Xs_input, Xs_numerical_absolute_output, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Desired', 'Achieved (Theory)', 'Achieved (Numerical)' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_absolute_decoded_ss_response_gain' ] ) 

% Plot the relative decoded steady state behavior.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Absolute Decoded Steady State Response' ); hold on, grid on, rotate3d on, view( 45, 30 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Output, x2 [-]' ), title( 'Transmission: Relative Decoded Steady State Error' )
surf( Cs, Xs_input, Xs_desired_relative_output, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Xs_input, Xs_theoretical_relative_output, 'Edgecolor', 'None', 'Facecolor', 'g', 'Facealpha', 0.5 )
surf( Cs, Xs_input, Xs_numerical_relative_output, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Desired', 'Achieved (Theory)', 'Achieved (Numerical)' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_relative_decoded_ss_response_gain' ] ) 


%% Plot the Encoded Error vs Gain.

% Plot the encoded error vs gain.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error' ); 
subplot( 2, 1, 1 ), hold on, grid on, rotate3d on, view( 25, 10 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, U1 [mv]' ), zlabel( 'Encoded Error, E [mV]' ), title( 'Transmission: Theoretical Encoded Error' )
surf( Cs, scale*Us_input, scale*errors_theoretical_encoded_absolute, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, scale*Us_input, scale*errors_theoretical_encoded_relative, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )

subplot( 2, 1, 2 ), hold on, grid on, rotate3d on, view( 25, 10 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, U1 [mv]' ), zlabel( 'Encoded Error, E [mV]' ), title( 'Transmission: Numerical Encoded Error' )
surf( Cs, scale*Us_input, scale*errors_numerical_encoded_absolute, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, scale*Us_input, scale*errors_numerical_encoded_relative, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_gain' ] ) 

% Plot the encoded error vs gain summary.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Summary' );
subplot( 2, 2, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error, E [mV]' ), title( 'Transmission: Encoded Theoretical Absolute Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ scale*errors_min_theoretical_encoded_absolute; flipud( scale*errors_max_theoretical_encoded_absolute ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, scale*errors_rmse_theoretical_encoded_absolute, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, scale*errors_min_theoretical_encoded_absolute, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, scale*errors_max_theoretical_encoded_absolute, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 2, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error, E [mV]' ), title( 'Transmission: Encoded Theoretical Relative Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ scale*errors_min_theoretical_encoded_relative; flipud( scale*errors_max_theoretical_encoded_relative ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, scale*errors_rmse_theoretical_encoded_relative, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, scale*errors_min_theoretical_encoded_relative, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, scale*errors_max_theoretical_encoded_relative, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 2, 3 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error, E [mV]' ), title( 'Transmission: Encoded Numerical Absolute Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ scale*errors_min_numerical_encoded_absolute; flipud( scale*errors_max_numerical_encoded_absolute ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, scale*errors_rmse_numerical_encoded_absolute, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, scale*errors_min_numerical_encoded_absolute, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, scale*errors_max_numerical_encoded_absolute, '--', 'Color', color2, 'Linewidth', 1 )

subplot( 2, 2, 4 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error, E [mV]' ), title( 'Transmission: Encoded Numerical Relative Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ scale*errors_min_numerical_encoded_relative; flipud( scale*errors_max_numerical_encoded_relative ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, scale*errors_rmse_numerical_encoded_relative, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, scale*errors_min_numerical_encoded_relative, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, scale*errors_max_numerical_encoded_relative, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_gain_summary' ] ) 


%% Plot the Decoded Error vs Gain.

% Plot the decoded error vs gain.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error' ); 
subplot( 2, 1, 1 ), hold on, grid on, rotate3d on, view( 25, 10 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Error, E [-]' ), title( 'Transmission: Theoretical Decoded Error' )
surf( Cs, Xs_input, errors_theoretical_decoded_absolute, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Xs_input, errors_theoretical_decoded_relative, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )

subplot( 2, 1, 2 ), hold on, grid on, rotate3d on, view( 25, 10 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Error, E [-]' ), title( 'Transmission: Numerical Decoded Error' )
surf( Cs, Xs_input, errors_numerical_decoded_absolute, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Xs_input, errors_numerical_decoded_relative, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_gain' ] ) 

% Plot the decoded error vs gain summary.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Summary' );
subplot( 2, 2, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error, E [-]' ), title( 'Transmission: Decoded Theoretical Absolute Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_theoretical_decoded_absolute; flipud( errors_max_theoretical_decoded_absolute ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_rmse_theoretical_decoded_absolute, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_theoretical_decoded_absolute, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_theoretical_decoded_absolute, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 2, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error, E [-]' ), title( 'Transmission: Decoded Theoretical Relative Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_theoretical_decoded_relative; flipud( errors_max_theoretical_decoded_relative ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_rmse_theoretical_decoded_relative, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_theoretical_decoded_relative, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_theoretical_decoded_relative, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 2, 3 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error, E [-]' ), title( 'Transmission: Decoded Numerical Absolute Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_numerical_decoded_absolute; flipud( errors_max_numerical_decoded_absolute ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_rmse_numerical_decoded_absolute, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_numerical_decoded_absolute, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_numerical_decoded_absolute, '--', 'Color', color2, 'Linewidth', 1 )

subplot( 2, 2, 4 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error, E [-]' ), title( 'Transmission: Decoded Numerical Relative Error Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_numerical_decoded_relative; flipud( errors_max_numerical_decoded_relative ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_rmse_numerical_decoded_relative, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_numerical_decoded_relative, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_numerical_decoded_relative, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_gain_summary' ] ) 


%% Plot the Encoded Error Difference vs Gain.

% Plot the encoded error difference.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Difference' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, U1 [mV]' ), zlabel( 'Encoded Error Difference, dE [mV]' ), title( 'Transmission: Encoded Error Difference' )
surf( Cs, Us_input, errors_diff_theoretical_encoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Us_input, errors_diff_numerical_encoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_difference_gain' ] ) 

% Plot the encoded error difference percentage.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Difference Percentage' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, U1 [mV]' ), zlabel( 'Encoded Error Difference Percentage, dE [%]' ), title( 'Transmission: Encoded Error Difference Percentage' )
surf( Cs, Us_input, errors_percent_diff_theoretical_encoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Us_input, errors_percent_diff_numerical_encoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_difference_percentage_gain' ] ) 

% Plot the encoded error difference.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Difference Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Difference, dE [mV]' ), title( 'Transmission: Encoded Theoretical Error Difference  Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_diff_theoretical_encoded; flipud( errors_max_diff_theoretical_encoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_diff_theoretical_encoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_diff_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_diff_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Difference Summary, dE [mV]' ), title( 'Transmission: Encoded Numerical Error Difference Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_diff_numerical_encoded; flipud( errors_max_diff_numerical_encoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_diff_numerical_encoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_diff_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_diff_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_difference_gain_summary' ] ) 

% Plot the encoded error percentage difference.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Percentage Difference Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Percentage Difference, E [%]' ), title( 'Transmission: Encoded Theoretical Error Percentage Difference Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_diff_theoretical_encoded; flipud( errors_max_percent_diff_theoretical_encoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_diff_theoretical_encoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_percent_diff_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_percent_diff_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Percentage Difference, E [%]' ), title( 'Transmission: Encoded Numerical Error Percentage Difference Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_diff_numerical_encoded; flipud( errors_max_percent_diff_numerical_encoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_diff_numerical_encoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_percent_diff_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_percent_diff_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_difference_percentage_gain_summary' ] ) 


%% Plot Decoded Error Difference vs Gain.

% Plot the decoded error difference.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Difference' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Error Difference, E [-]' ), title( 'Transmission: Decoded Error Difference' )
surf( Cs, Xs_input, errors_diff_theoretical_decoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Xs_input, errors_diff_numerical_decoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_difference_gain' ] ) 

% Plot the decoded error difference percentage.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Difference Percentage' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Error Difference Percentage, E [%]' ), title( 'Transmission: Decoded Error Difference Percentage' )
surf( Cs, Xs_input, errors_percent_diff_theoretical_decoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Xs_input, errors_percent_diff_numerical_decoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_difference_percentage_gain' ] ) 

% % Plot the decoded mse difference.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded MSE Difference' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded MSE Difference, mse [-]' ), title( 'Transmission: Decoded MSE Difference' )
% plot( cs, errors_mse_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_mse_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_mse_difference_gain' ] ) 
% 
% % Plot the decoded mse difference percentage.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded MSE Difference Percentage' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded MSE Difference Percentage, mse [%]' ), title( 'Transmission: Decoded MSE Difference Percentage' )
% plot( cs, errors_mse_percent_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_mse_percent_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_mse_percent_difference_gain' ] ) 
% 
% % Plot the decoded std difference.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded STD Difference' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded STD Difference, std [-]' ), title( 'Transmission: Decoded STD Difference' )
% plot( cs, errors_std_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_std_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_std_difference_gain' ] ) 
% 
% % Plot the decoded std percentage difference.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded STD Difference Percentage' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded STD Difference Percentage, std [%]' ), title( 'Transmission: Decoded STD Difference Percentage' )
% plot( cs, errors_std_percent_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_std_percent_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_std_percent_difference_gain' ] ) 
% 
% % Plot the decoded minimum error difference.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Minimum Error Difference' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Minimum Error Difference, Emin [-]' ), title( 'Transmission: Decoded Minimum Error Difference' )
% plot( cs, errors_min_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_min_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_min_error_difference_gain' ] ) 
% 
% % Plot the decoded minimum error difference percentage.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Minimum Error Difference Percentage' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Minimum Error Difference Percentage, Emin [%]' ), title( 'Transmission: Decoded Minimum Error Difference Percentage' )
% plot( cs, errors_min_percent_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_min_percent_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_min_percent_error_difference_gain' ] ) 
% 
% % Plot the decoded maximum error difference.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Maximum Error Difference' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Maximum Error Difference, Emin [-]' ), title( 'Transmission: Decoded Maximum Error Difference' )
% plot( cs, errors_max_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_max_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_max_error_difference_gain' ] ) 
% 
% % Plot the decoded maximum error difference percentage.
% fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Maximum Error Difference Percentage' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Maximum Error Difference Percentage, Emax [%]' ), title( 'Transmission: Decoded Maximum Error Difference Percentage' )
% plot( cs, errors_max_percent_diff_theoretical_decoded, '-.', 'Color', color1, 'Linewidth', 3 )
% plot( cs, errors_max_percent_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 3 )
% legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
% saveas( fig, [ save_directory, '\', 'transmission_decoded_max_percent_error_difference_gain' ] ) 

% Plot the decoded error difference.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Difference Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Difference, E [-]' ), title( 'Transmission: Decoded Theoretical Error Difference Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_diff_theoretical_decoded; flipud( errors_max_diff_theoretical_decoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_diff_theoretical_decoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_diff_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_diff_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Difference, E [-]' ), title( 'Transmission: Decoded Numerical Error Difference Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_diff_numerical_decoded; flipud( errors_max_diff_numerical_decoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_diff_numerical_decoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_error_difference_gain_summary' ] ) 

% Plot the decoded error percentage difference.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Percentage Difference Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Percentage Difference, E [%]' ), title( 'Transmission: Decoded Theoretical Error Percentage Difference Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_diff_theoretical_decoded; flipud( errors_max_percent_diff_theoretical_decoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_diff_theoretical_decoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_percent_diff_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_percent_diff_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Percentage Difference, E [%]' ), title( 'Transmission: Decoded Numerical Error Percentage Difference Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_diff_numerical_decoded; flipud( errors_max_percent_diff_numerical_decoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_diff_numerical_decoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_percent_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_percent_diff_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )


%% Plot the Encoded Error Improvement vs Gain.

% Plot the encoded error improvement.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Improvement' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, x1 [-]' ), zlabel( 'Encoded Error Improvement, E [-]' ), title( 'Transmission: Encoded Error Improvement' )
surf( Cs, Us_input, errors_improv_theoretical_encoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Us_input, errors_improv_numerical_encoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_improvement_gain' ] ) 

% Plot the encoded error improvement percentage.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Improvement Percentage' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Input, x1 [-]' ), zlabel( 'Encoded Error Improvement Percentage, E [%]' ), title( 'Transmission: Encoded Error Improvement Percentage' )
surf( Cs, Us_input, errors_percent_improv_theoretical_encoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Us_input, errors_percent_improv_numerical_encoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_improvement_percentage_gain' ] ) 

% Plot the encoded error improvement.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Improvement Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Improvement, E [-]' ), title( 'Transmission: Encoded Theoretical Error Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_improv_theoretical_encoded; flipud( errors_max_improv_theoretical_encoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_improv_theoretical_encoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_improv_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_improv_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Improvement, E [-]' ), title( 'Transmission: Encoded Numerical Error Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_improv_numerical_encoded; flipud( errors_max_improv_numerical_encoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_improv_numerical_encoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_improv_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_improv_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_difference_gain_summary' ] ) 

% Plot the encoded error percentage improvement.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Encoded Error Percentage Improvement Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Percentage Improvement, E [%]' ), title( 'Transmission: Encoded Theoretical Error Percentage Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_improv_theoretical_encoded; flipud( errors_max_percent_improv_theoretical_encoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_improv_theoretical_encoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_percent_improv_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_percent_improv_theoretical_encoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Encoded Error Percentage Improvement, E [%]' ), title( 'Transmission: Encoded Numerical Error Percentage Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_improv_numerical_encoded; flipud( errors_max_percent_improv_numerical_encoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_improv_numerical_encoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_percent_improv_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_percent_improv_numerical_encoded, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_encoded_error_improvement_percentage_gain_summary' ] ) 


%% Plot the Decoded Error Improvement vs Gain.

% Plot the decoded error improvement.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Improvement' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Error Improvement, E [-]' ), title( 'Transmission: Decoded Error Improvement' )
surf( Cs, Us_input, errors_improv_theoretical_decoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Us_input, errors_improv_numerical_decoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_improvement_gain' ] ) 

% Plot the encoded error improvement percentage.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Improvement Percentage' ); hold on, grid on, rotate3d on, view( 45, 20 ), xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Input, x1 [-]' ), zlabel( 'Decoded Error Improvement Percentage, E [%]' ), title( 'Transmission: Decoded Error Improvement Percentage' )
surf( Cs, Us_input, errors_percent_improv_theoretical_decoded, 'Edgecolor', 'None', 'Facecolor', 'b', 'Facealpha', 0.5 )
surf( Cs, Us_input, errors_percent_improv_numerical_decoded, 'Edgecolor', 'None', 'Facecolor', 'r', 'Facealpha', 0.5 )
legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_improvement_percentage_gain' ] ) 

% Plot the encoded error improvement.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Improvement Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Improvement, E [-]' ), title( 'Transmission: Decoded Theoretical Error Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_improv_theoretical_decoded; flipud( errors_max_improv_theoretical_decoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_improv_theoretical_decoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_improv_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_improv_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Improvement, E [-]' ), title( 'Transmission: Decoded Numerical Error Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_improv_numerical_decoded; flipud( errors_max_improv_numerical_decoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_improv_numerical_decoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_improv_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_improv_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_difference_gain_summary' ] ) 

% Plot the encoded error percentage improvement.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Decoded Error Percentage Improvement Summary' );
subplot( 2, 1, 1 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Percentage Improvement, E [%]' ), title( 'Transmission: Decoded Theoretical Error Percentage Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_improv_theoretical_decoded; flipud( errors_max_percent_improv_theoretical_decoded ) ], color1( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_improv_theoretical_decoded, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, errors_min_percent_improv_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )
plot( cs, errors_max_percent_improv_theoretical_decoded, '--', 'Color', color1, 'Linewidth', 1 )

subplot( 2, 1, 2 ), hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Decoded Error Percentage Improvement, E [%]' ), title( 'Transmission: Decoded Numerical Error Percentage Improvement Summary' )
patch( [ cs'; flipud( cs' ) ], [ errors_min_percent_improv_numerical_decoded; flipud( errors_max_percent_improv_numerical_decoded ) ], color2( 1:end - 1 ), 'FaceAlpha', 0.5, 'EdgeColor', 'None' )
plot( cs, errors_mse_percent_improv_numerical_decoded, '-', 'Color', color2, 'Linewidth', 3 )
plot( cs, errors_min_percent_improv_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )
plot( cs, errors_max_percent_improv_numerical_decoded, '--', 'Color', color2, 'Linewidth', 1 )
saveas( fig, [ save_directory, '\', 'transmission_decoded_error_improvement_percentage_gain_summary' ] ) 


%% Plot the Numerical Stability Information vs Gain.

% Plot the maximum RK4 step size.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Maximum RK4 Step Size' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Max RK4 Step Size, dt [s]' ), title( 'Transmission: Maximum RK4 Step Size' )
plot( cs, dts_max_absolute, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, dts_max_relative, '-', 'Color', color2, 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_max_rk4_step_size_gain' ] ) 

% Plot the maximum condition number.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Maximum Condition Number' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Max Condition Number [-]' ), title( 'Transmission: Maximum Condition Number' )
plot( cs, condition_numbers_max_absolute, '-', 'Color', color1, 'Linewidth', 3 )
plot( cs, condition_numbers_max_relative, '-', 'Color', color2, 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_max_condition_number_gain' ] ) 


%% Plot Network Parameters vs Gain.

% Plot the sodium channel conductance vs gain.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Sodium Channel Conductance' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Sodium Channel Conductance, Gna [S]' ), title( 'Transmission: Sodium Channel Conductance' )
plot( cs, Gnas_absolute( :, 1 ), '-.', 'Color', color1, 'Linewidth', 3 )
plot( cs, Gnas_absolute( :, 1 ), '--', 'Color', color1, 'Linewidth', 3 )
plot( cs, Gnas_relative( :, 1 ), '-.', 'Color', color2, 'Linewidth', 3 )
plot( cs, Gnas_relative( :, 1 ), '--', 'Color', color2, 'Linewidth', 3 )
legend( { 'Absolute 1', 'Absolute 2', 'Relative 1', 'Relative 2' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_sodium_channel_conductance_gain' ] )

% Plot the maximum membrane voltage vs gain.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Maximum Membrane Voltage' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Maximum Membrane Voltage, R [V]' ), title( 'Transmission: Maximum Membrane Voltage' )
plot( cs, Rs2_absolute, '-.', 'Color', color1, 'Linewidth', 3 )
plot( cs, Rs2_relative, '-.', 'Color', color2, 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_maximum_membrane_voltage_gain' ] )

% Plot the synaptic reversal potential vs gain.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Synaptic Reversal Potential' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Synaptic Reversal Potential, dE [V]' ), title( 'Transmission: Synaptic Reversal Potential' )
plot( cs, dEs21_absolute, '-.', 'Color', color1, 'Linewidth', 3 )
plot( cs, dEs21_relative, '-.', 'Color', color2, 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_synaptic_reversal_potential_gain' ] )

% Plot the maximum synaptic conductance vs gain.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Maximum Synaptic Conductance' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Maximum Synaptic Conductance, gs [S]' ), title( 'Transmission: Maximum Synaptic Conductance' )
plot( cs, gs21_absolute, '-.', 'Color', color1, 'Linewidth', 3 )
plot( cs, gs21_relative, '-.', 'Color', color2, 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_max_synaptic_conductance_gain' ] )

% Plot the steady state applied current vs gain.
fig = figure( 'Color', 'w', 'Name', 'Transmission: Applied Current' ); hold on, grid on, xlabel( 'Gain, c [-]' ), ylabel( 'Applied Current, Ia [A]' ), title( 'Transmission: Applied Current' )
plot( cs, Ias2_absolute, '-.', 'Color', color1, 'Linewidth', 3 )
plot( cs, Ias2_relative, '-.', 'Color', color2, 'Linewidth', 3 )
legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
saveas( fig, [ save_directory, '\', 'transmission_applied_currents_gain' ] )

