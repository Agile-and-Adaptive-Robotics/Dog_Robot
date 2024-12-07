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
        
        
        % Implement a function to plot the steady state error of a subnetwork for a specific gain.
        function fig = plot_steady_state_error_comparison( ~, xs_absolute, errors_theoretical_absolute, errors_numerical_absolute, color_absolute, xs_relative, errors_theoretical_relative, errors_numerical_relative, color_relative, scale, subnetwork_name, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 17, save_directory = './'; end
            if nargin < 16, save_flag = true; end
            if nargin < 15, unit = 'mV'; end
            if nargin < 14, output_variable_string = 'dU'; end
            if nargin < 13, input_variable_string = 'U1'; end
            if nargin < 12, encoded_string = 'Encoded'; end
            if nargin < 11, subnetwork_name = 'Transmission'; end
            if nargin < 10, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( 'Absolute vs Relative %s: %s Steady State Error', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Error, %s [%s]', encoded_string, output_variable_string, unit );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the absolute and relative theoretical and numerical errors.
            plot( scale*xs_absolute, scale*errors_theoretical_absolute, '-.', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_absolute, scale*errors_numerical_absolute, '--', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_relative, scale*errors_theoretical_relative, '-.', 'Color', color_relative, 'Linewidth', 3 )
            plot( scale*xs_relative, scale*errors_numerical_relative, '--', 'Color', color_relative, 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Absolute Theoretical', 'Absolute Numerical', 'Relative Theoretical', 'Relative Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_steady_state_error.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the steady state error percentage of a subnetwork for a specific gain.
        function fig = plot_steady_state_error_percentage_comparison( ~, xs_absolute, error_percentages_theoretical_absolute, error_percentages_numerical_absolute, color_absolute, xs_relative, error_percentages_theoretical_relative, error_percentages_numerical_relative, color_relative, scale, subnetwork_name, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 17, save_directory = './'; end
            if nargin < 16, save_flag = true; end
            if nargin < 15, unit = 'mV'; end
            if nargin < 14, output_variable_string = 'dU'; end
            if nargin < 13, input_variable_string = 'U1'; end
            if nargin < 12, encoded_string = 'Encoded'; end
            if nargin < 11, subnetwork_name = 'Transmission'; end
            if nargin < 10, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( 'Absolute vs Relative %s: %s Steady State Error Percentage', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Error, %s [%%]', encoded_string, output_variable_string );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the absolute and relative theoretical and numerical errors.
            plot( scale*xs_absolute, error_percentages_theoretical_absolute, '-.', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_absolute, error_percentages_numerical_absolute, '--', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_relative, error_percentages_theoretical_relative, '-.', 'Color', color_relative, 'Linewidth', 3 )
            plot( scale*xs_relative, error_percentages_numerical_relative, '--', 'Color', color_relative, 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Absolute Theoretical', 'Absolute Numerical', 'Relative Theoretical', 'Relative Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_steady_state_error_percentage.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the steady state error difference of a subnetwork for a specific gain.
        function fig = plot_steady_state_error_difference( ~, xs_theoretical, error_difference_theoretical, xs_numerical, error_difference_numerical, scale, subnetwork_name, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 13, save_directory = './'; end
            if nargin < 12, save_flag = true; end
            if nargin < 11, unit = 'mV'; end
            if nargin < 10, output_variable_string = 'dU'; end
            if nargin < 9, input_variable_string = 'U1'; end
            if nargin < 8, encoded_string = 'Encoded'; end
            if nargin < 7, subnetwork_name = 'Transmission'; end
            if nargin < 6, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( '%s: %s Steady State Error Difference', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Error Difference, %s [%s]', encoded_string, output_variable_string, unit );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the absolute and relative theoretical and numerical errors.
            plot( scale*xs_theoretical, scale*error_difference_theoretical, '-.', 'Linewidth', 3 )
            plot( scale*xs_numerical, scale*error_difference_numerical, '--', 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal')
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_steady_state_error_difference.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the steady state error difference of a subnetwork for a specific gain.
        function fig = plot_steady_state_error_percentage_difference( ~, xs_theoretical, error_percentages_difference_theoretical, xs_numerical, error_percentages_difference_numerical, scale, subnetwork_name, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 13, save_directory = './'; end
            if nargin < 12, save_flag = true; end
            if nargin < 11, unit = 'mV'; end
            if nargin < 10, output_variable_string = 'dU'; end
            if nargin < 9, input_variable_string = 'U1'; end
            if nargin < 8, encoded_string = 'Encoded'; end
            if nargin < 7, subnetwork_name = 'Transmission'; end
            if nargin < 6, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( '%s: %s Steady State Error Percentage Difference', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Error Percentage Difference, %s [%%]', encoded_string, output_variable_string );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the absolute and relative theoretical and numerical errors.
            plot( scale*xs_theoretical, error_percentages_difference_theoretical, '-.', 'Linewidth', 3 )
            plot( scale*xs_numerical, error_percentages_difference_numerical, '--', 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal')
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_steady_state_error_percentage_difference.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the steady state error improvement of a subnetwork for a specific gain.
        function fig = plot_steady_state_error_improvement( ~, xs_theoretical, error_improvement_theoretical, xs_numerical, error_improvement_numerical, scale, subnetwork_name, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 13, save_directory = './'; end
            if nargin < 12, save_flag = true; end
            if nargin < 11, unit = 'mV'; end
            if nargin < 10, output_variable_string = 'dU'; end
            if nargin < 9, input_variable_string = 'U1'; end
            if nargin < 8, encoded_string = 'Encoded'; end
            if nargin < 7, subnetwork_name = 'Transmission'; end
            if nargin < 6, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( '%s: %s Steady State Error Improvement', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Error Improvement, %s [%s]', encoded_string, output_variable_string, unit );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the absolute and relative theoretical and numerical errors.
            plot( scale*xs_theoretical, scale*error_improvement_theoretical, '-.', 'Linewidth', 3 )
            plot( scale*xs_numerical, scale*error_improvement_numerical, '--', 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal')
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_steady_state_error_improvement.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the steady state error improvement of a subnetwork for a specific gain.
        function fig = plot_steady_state_error_percentage_improvement( ~, xs_theoretical, error_percentages_improvement_theoretical, xs_numerical, error_percentages_improvement_numerical, scale, subnetwork_name, encoded_string, input_variable_string, output_variable_string, unit, save_flag, save_directory )
            
            % Set the default input arguments.
            if nargin < 13, save_directory = './'; end
            if nargin < 12, save_flag = true; end
            if nargin < 11, unit = 'mV'; end
            if nargin < 10, output_variable_string = 'dU'; end
            if nargin < 9, input_variable_string = 'U1'; end
            if nargin < 8, encoded_string = 'Encoded'; end
            if nargin < 7, subnetwork_name = 'Transmission'; end
            if nargin < 6, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( '%s: %s Steady State Error Percentage Improvement', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( '%s Error Percentage Improvement, %s [%%]', encoded_string, output_variable_string );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the absolute and relative theoretical and numerical errors.
            plot( scale*xs_theoretical, error_percentages_improvement_theoretical, '-.', 'Linewidth', 3 )
            plot( scale*xs_numerical, error_percentages_improvement_numerical, '--', 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Theoretical', 'Numerical' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal')
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_%s_steady_state_error_percentage_improvement.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the maximum rk4 step size for a specific gain.
        function fig = plot_rk4_maximum_timestep( ~, xs_absolute, dts_absolute, color_absolute, xs_relative, dts_relative, color_relative, scale, subnetwork_name, encoded_string, input_variable_string, unit, save_flag, save_directory )
        
            % Set the default input arguments.
            if nargin < 14, save_directory = './'; end
            if nargin < 13, save_flag = true; end
            if nargin < 12, unit = 'mV'; end
            if nargin < 11, input_variable_string = 'U1'; end
            if nargin < 10, encoded_string = 'Encoded'; end
            if nargin < 9, subnetwork_name = 'Transmission'; end
            if nargin < 8, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( '%s: %s RK4 Maximum Timestep', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( 'RK4 Maximum Timestep, dt [s]' );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the desired, theoretical, and numerical steady state responses.
            plot( scale*xs_absolute, dts_absolute, '-', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_relative, dts_relative, '-', 'Color', color_relative, 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_rk4_maximum_timestep_%s.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        % Implement a function to plot the condition number for a specific gain.
        function fig = plot_condition_numbers( ~, xs_absolute, condition_numbers_absolute, color_absolute, xs_relative, condition_numbers_relative, color_relative, scale, subnetwork_name, encoded_string, input_variable_string, unit, save_flag, save_directory )
        
            % Set the default input arguments.
            if nargin < 14, save_directory = './'; end
            if nargin < 13, save_flag = true; end
            if nargin < 12, unit = 'mV'; end
            if nargin < 11, input_variable_string = 'U1'; end
            if nargin < 10, encoded_string = 'Encoded'; end
            if nargin < 9, subnetwork_name = 'Transmission'; end
            if nargin < 8, scale = 1; end
            
            % Compute the figure labels.
            title_string = sprintf( '%s: %s Condition Numbers', subnetwork_name, encoded_string );
            xlabel_string = sprintf( '%s Input, %s [%s]', encoded_string, input_variable_string, unit );
            ylabel_string = sprintf( 'Condition Numbers [-]' );

            % Create the figure.
            fig = figure( 'Color', 'w', 'Name', title_string ); hold on, grid on, xlabel( xlabel_string ), ylabel( ylabel_string ), title( title_string )
            
            % Plot the desired, theoretical, and numerical steady state responses.
            plot( scale*xs_absolute, condition_numbers_absolute, '-', 'Color', color_absolute, 'Linewidth', 3 )
            plot( scale*xs_relative, condition_numbers_relative, '-', 'Color', color_relative, 'Linewidth', 3 )
            
            % Add a legend to the figure.
            legend( { 'Absolute', 'Relative' }, 'Location', 'Bestoutside', 'Orientation', 'Horizontal' )
            
            % Determine whether to save the figure.
            if save_flag                            % If we want to save the figure...
                    
                % Define the file name.
                file_name = sprintf( '%s_condition_number_%s.png', lower( subnetwork_name ), lower( encoded_string ) );
                
                % Save the figure.
                saveas( fig, [ save_directory, '\', file_name ] ) 
            
            end
            
        end
        
        
        
    end
    
    
end
    