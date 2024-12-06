classdef plotting_utilities_class
    
    % This class contains properties and methods related to plotting utilities.
    
    
    %% PLOTTING UTILITIES PROPERTIES
    
    % Define the class properties.
    properties
        
        
        
    end
    
    
    %% PLOTTING UTILITIES METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = plotting_utilities_class(  )
            
            
        end
        
        
        %% Plotting Functions.
        
        % Implement a function to compute the number of subplot rows and columns necessary to store a certain number of plots.
        function [ nrows, ncols ] = get_subplot_rows_columns( ~, n, bPreferColumns )
            
            %If the value is not supplied by the user, default to preferring to add rows.
            if nargin < 3, bPreferColumns = false; end
            
            %Determine whether n is an integer.
            if n ~= round( n )                                                            % If n is not an integer...
                
                % Round n to the nearest integer.
                n = round( n );                                                           
                
                % Throw a warning about rounding n.
                warning( 'n must be an integer.  Rounding n to the nearest integer.' )    
            
            end
            
            % Compute the square root of the integer of interest.
            nsr = sqrt( n );
            
            % Determine how many rows and columns to use in the subplot.
            if nsr == round( nsr )                    % If n is a perfect square...
            
                [ nrows, ncols ] = deal( nsr );         % Set the number of rows and columns to be the square root.
            
            else
                
                % Compute all divisors of n.
                dn = divisors( n );
                
                % Set the number of rows and columns to be the central divisors.
                nrows = dn( length( dn ) / 2 + 1 );
                ncols = dn( length( dn ) / 2 );
                
            end
            
            % Determine whether to give prefernce to columns or rows.
            if bPreferColumns                           % If we want to prefer columns...
                
                % Flip the row and column assignments, since we prefer rows by default.
                [ nrows, ncols ] = deal( ncols, nrows );
                
            end
            
        end
        
        
        % Implement a function to plot the steady state response of a subnetwork for a specific encoding scheme and gain.
        function fig = plot_steady_state_response( ~, xs, ys_desired, ys_theoretical, ys_numerical, scale, subnetwork_name, encoding_scheme, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 14, save_directory = './'; end
            if nargin < 13, save_flag = true; end
            if nargin < 12, unit = 'mV'; end
            if nargin < 11, output_variable_string = 'U2'; end
            if nargin < 10, input_variable_string = 'U1'; end
            if nargin < 9, encoded_string = 'Encoded'; end
            if nargin < 8, encoding_scheme = 'Absolute'; end
            if nargin < 7, subnetwork_name = 'Transmission'; end
            if nargin < 6, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( '%s %s: %s Steady State Response', encoding_scheme, subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Output, %s [%s]', encoded_string, output_variable_string, unit );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the desired, theoretical, and numerical steady state responses.
            plot( scale*xs, scale*ys_desired, '-', 'Linewidth', 3 )
            plot( scale*xs, scale*ys_theoretical, '-.', 'Linewidth', 3 )
            plot( scale*xs, scale*ys_numerical, '--', 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Desired', 'Achieved (Theory)', 'Achieved (Numerical)' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_%s_steady_state_response.png', lower( encoding_scheme ), lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the steady state response of a subnetwork for a specific encoding scheme and gain.
        function fig = plot_steady_state_response_comparison( ~, xs_absolute, ys_desired_absolute, ys_theoretical_absolute, ys_numerical_absolute, color_absolute, xs_relative, ys_desired_relative, ys_theoretical_relative, ys_numerical_relative, color_relative, scale, subnetwork_name, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 19, save_directory = './'; end
            if nargin < 18, save_flag = true; end
            if nargin < 17, unit = 'mV'; end
            if nargin < 16, output_variable_string = 'U2'; end
            if nargin < 15, input_variable_string = 'U1'; end
            if nargin < 14, encoded_string = 'Encoded'; end
            if nargin < 13, subnetwork_name = 'Transmission'; end
            if nargin < 12, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( 'Absolute vs Relative %s: %s Steady State Response', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Output, %s [%s]', encoded_string, output_variable_string, unit );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the absolute desired, theoretical, and numerical steady state responses.
            plot( scale*xs_absolute, scale*ys_desired_absolute, '-', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_absolute, scale*ys_theoretical_absolute, '-.', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_absolute, scale*ys_numerical_absolute, '--', 'Color', color_absolute, 'Linewidth', 3 )
            
            % Plot the relative desired, theoretical, and numerical steady state responses.
            plot( scale*xs_relative, scale*ys_desired_relative, '-', 'Color', color_relative, 'Linewidth', 3 )
            plot( scale*xs_relative, scale*ys_theoretical_relative, '-.', 'Color', color_relative, 'Linewidth', 3 )
            plot( scale*xs_relative, scale*ys_numerical_relative, '--', 'Color', color_relative, 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Absolute Desired', 'Absolute Achieved (Theory)', 'Absolute Achieved (Numerical)', 'Relative Desired', 'Relative Achieved (Theory)', 'Relative Achieved (Numerical)' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_steady_state_response.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
    end
    
    
end
    