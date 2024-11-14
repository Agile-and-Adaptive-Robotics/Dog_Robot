classdef applied_voltage_manager_class
    
    % This class contains properties and methods related to managing applied voltages.
    
    
    %% APPLIED VOLTAGE MANAGER PROPERTIES
    
    % Define general class properties.
    properties
        
        applied_voltages
        num_applied_voltages
        
        array_utilities
        data_loader_utilities
        
    end
    
    
    % Define private, constant class properties.
    properties ( Access = private, Constant = true )

        
        
    end
    
    
    %% APPLIED VOLTAGE MANAGER METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = applied_voltage_manager_class( applied_voltages )
            
            % Create an instance of the data loader utilities class.
            self.data_loader_utilities = data_loader_utilities_class(  );
            
            % Create an instance of the array manager class.
            self.array_utilities = array_utilities_class(  );
            
            % Set the default properties.
            if nargin < 1, self.applied_voltages = [  ]; else, self.applied_voltages = applied_voltages; end
            
            % Compute the number of applied voltages.
            self.num_applied_voltages = length( self.applied_voltages );
            
        end
        
        
        %% General Get & Set Applied Voltage Property Functions.
        
        % Implement a function to retrieve the properties of specific applied voltages.
        function xs = get_applied_voltage_property( self, applied_voltage_IDs, applied_voltage_property, as_matrix_flag, applied_voltages, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 5, applied_voltages = self.applied_voltages; end
            if nargin < 4, as_matrix_flag = self.as_matrix_flag_DEFAULT; end
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Determine how many applied voltages to which we are going to apply the given method.
            num_properties_to_get = length( applied_voltage_IDs );
            
            % Preallocate a variable to store the applied voltage properties.
            xs = cell( 1, num_properties_to_get );
            
            % Retrieve the given applied voltage property for each applied voltage.
            for k = 1:num_properties_to_get                                                                         % Iterate through each of the properties to get...
                
                % Retrieve the index associated with this applied voltage ID.
                applied_voltage_index = self.get_applied_voltage_index( applied_voltage_IDs( k ), applied_voltages, undetected_option );
                
                % Define the eval string.
                eval_str = sprintf( 'xs{ k } = applied_voltages( %0.0f ).%s;', applied_voltage_index, applied_voltage_property );
                
                % Evaluate the given applied voltage property.
                eval( eval_str );
                
            end
            
            % Determine whether to convert the network properties to a matrix.
            if as_matrix_flag                                    % If we want the applied voltage properties as a matrix instead of a cell...
                
                % Convert the applied current properties from a cell to a matrix.
                xs = cell2mat( xs );
                
            end
            
        end
        
        
        % Implement a function to set the properties of specific applied voltages.
        function [ applied_voltages, self ] = set_applied_voltage_property( self, applied_voltage_IDs, applied_voltage_property_values, applied_voltage_property, applied_voltages, set_flag )
            
            % Set the default input arguments.
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end
            if nargin < 5, applied_voltages = self.applied_voltages; end
            
            % Compute the number of applied voltages.
            n_applied_voltages = length( applied_voltages );
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Validate the applied voltage property values.
            if ~isa( applied_voltage_property_values, 'cell' )                    % If the applied voltage property values are not a cell array...
                
                % Convert the applied voltage property values to a cell array.
                applied_voltage_property_values = num2cell( applied_voltage_property_values );
                
            end
            
            % Retreive the number of applied voltage IDs.
            num_applied_voltage_IDs = length( applied_voltage_IDs );
            %             num_applied_voltage_IDs = size( applied_voltage_IDs, 2 );
            
            % Retrieve the number of applied voltage property values.
            num_applied_voltage_property_values = length( applied_voltage_property_values );
            
            % Ensure that the provided neuron property values have the same length as the provided applied voltage IDs.
            if ( num_applied_voltage_IDs ~= num_applied_voltage_property_values )                                     % If the number of provided applied voltage IDs does not match the number of provided property values...
                
                % Determine whether to agument the property values.
                if num_applied_voltage_property_values == 1                                                  % If there is only one provided property value...
                    
                    % Agument the property value length to match the ID length.
                    %                     applied_voltage_property_values = applied_voltage_property_values*ones( 1, num_applied_voltage_IDs );
                    applied_voltage_property_values = repmat( applied_voltage_property_values, [ 1, num_applied_voltage_IDs ] );
                    
                else                                                                                % Otherwise...
                    
                    % Throw an error.
                    error( 'The number of provided applied voltage propety values must match the number of provided applied voltage IDs, unless a single applied voltage property value is provided.' )
                    
                end
                
            end
            
            
            % Set the properties of each applied voltage.
            for k = 1:n_applied_voltages                   % Iterate through each applied voltage...
                
                % Determine the index of the applied voltage property value that we want to apply to this applied voltage (if we want to set a property of this applied voltage).
                index = find( applied_voltages( k ).ID == applied_voltage_IDs, 1 );
                
                % Determine whether to set a property of this applied voltage.
                if ~isempty( index )                         % If a matching applied voltage ID was detected...
                    
                    % Create an evaluation string that sets the desired applied voltage property.
                    eval_string = sprintf( 'applied_voltages( %0.0f ).%s = applied_voltage_property_values{ %0.0f };', k, applied_voltage_property, index );
                    
                    % Evaluate the evaluation string.
                    eval( eval_string );
                    
                end
                
            end
            
            % Determine whether to update the applied voltages manager object.
            if set_flag, self.applied_voltages = applied_voltages; end
            
        end
        
        
        %% Specific Get & Set Functions.
        
        % Implement a function to retrieve all of the neuron IDs.
        function applied_voltage_IDs = get_all_applied_voltage_IDs( self, applied_voltages )
            
            % Set the default input arguments.
            if nargin < 2, applied_voltages = self.applied_voltages; end
            
            % Compute the number of applied voltages.
            n_applied_voltages = length( applied_voltages );
            
            % Preallocate a variable to store the applied voltage IDs.
            applied_voltage_IDs = zeros( 1, n_applied_voltages );
            
            % Retrieve the ID associated with each applied voltage.
            for k = 1:n_applied_voltages                             % Iterate through each of the applied voltages...
                
                % Retrieve this applied voltage ID.
                applied_voltage_IDs( k ) = applied_voltages( k ).ID;
                
            end
            
        end
        
        
        % Implement a function to retrieve the number of time steps of the specified applied voltages.
        function num_timesteps = get_num_timesteps( self, applied_voltage_IDs, applied_voltages, filter_disabled_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 5, process_option = self.process_option_DEFAULT; end
            if nargin < 4, filter_disabled_flag = self.filter_distabled_flag_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Determine how to compute the number of timesteps.
            if all( applied_voltage_IDs == -1 )                         % If all of the applied voltage IDs are invalid...
                
                % Set the number of timesteps to zero.
                num_timesteps = 0;
                
            else                                                        % Otherwise...
                
                % Remove any invalid applied voltage IDs.
                applied_voltage_IDs( applied_voltage_IDs == -1 ) = [  ];
                
                % Remove IDs associated with disabled applied voltages (if desired).
                if filter_disabled_flag, applied_voltage_IDs = self.remove_disabled_applied_voltage_IDs( applied_voltage_IDs, applied_voltages, undetected_option ); end
                
                % Retrieve the number of timesteps associated with each applied voltage.
                num_timesteps = self.get_applied_voltage_property( applied_voltage_IDs, 'num_timesteps', as_matrix_flag, applied_voltages, undetected_option );
                
                % Determine how to process the number of timesteps.
                if strcmpi( process_option, 'average' )                     % If we want the average time step...
                    
                    % Set the number of timesteps to be the average number of timesteps.
                    num_timesteps = mean( num_timesteps );
                    
                elseif strcmpi( process_option, 'max' )                    % If we want the maximum time step...
                    
                    % Set the number of timesteps to be the largest number of timesteps.
                    num_timesteps = max( num_timesteps );
                    
                elseif strcmpi( process_option, 'min' )                    % If we want the minimum time step...
                    
                    % Set the number of timesteps to be the smallest number of timesteps.
                    num_timesteps = min( num_timesteps );
                    
                elseif strcmpi( process_option, 'none' )                   % If we have selected no process options...
                    
                    % Do nothing.
                    
                else                                                        % Otherwise...
                    
                    % Throw an error.
                    error( 'Process option %s not recognized.', process_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to retrieve the step size of the specified applied voltages.
        function dt = get_dts( self, applied_voltage_IDs, applied_voltages, filter_disabled_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 5, process_option = self.process_option_DEFAULT; end
            if nargin < 4, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Determine how to compute the step size.
            if all( applied_voltage_IDs == -1 )                         % If all of the applied voltage IDs are invalid...
                
                % Set the step size to zero.
                dt = 1e-3;
                
            else                                                        % Otherwise...
                
                % Remove any invalid applied voltage IDs.
                applied_voltage_IDs( applied_voltage_IDs == -1 ) = [  ];
                
                % Remove IDs associated with disabled applied voltages (if desired).
                if filter_disabled_flag, applied_voltage_IDs = self.remove_disabled_applied_voltage_IDs( applied_voltage_IDs, applied_voltages, undetected_option ); end
                
                % Retrieve the step size associated with each applied voltage.
                dt = self.get_applied_voltage_property( applied_voltage_IDs, 'dt', as_matrix_flag, applied_voltages, undetected_option );
                
                % Determine whether the step size needs to be set.
                if isempty( dt ), dt = 1e-3; end
                
                % Determine how to process the step size.
                if strcmpi( process_option, 'average' )                     % If we want the average step size...
                    
                    % Set the step size to be the average step size.
                    dt = mean( dt );
                    
                elseif strcmpi( process_option, 'max' )                    % If we want the maximum step size...
                    
                    % Set the step size to be the largest step size.
                    dt = max( dt );
                    
                elseif strcmpi( process_option, 'min' )                    % If we want the minimum step size...
                    
                    % Set the step size to be the smallest step size.
                    dt = min( dt );
                    
                elseif strcmpi( process_option, 'none' )                   % If we have selected no process options...
                    
                    % Do nothing.
                    
                else                                                        % Otherwise...
                    
                    % Throw an error.
                    error( 'Process option %s not recognized.', process_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to retrieve the final time of the specified applied voltages.
        function tf = get_tfs( self, applied_voltage_IDs, applied_voltages, filter_disabled_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 5, process_option = self.process_option_DEFAULT; end
            if nargin < 4, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Determine how to compute the final time.
            if all( applied_voltage_IDs == -1 )             % If all of the applied voltage IDs are invalid...
                
                % Set the final time to zero.
                tf = 0;
                
            else                                            % Otherwise...
                
                % Remove any invalid applied voltage IDs.
                applied_voltage_IDs( applied_voltage_IDs == -1 ) = [  ];
                
                % Remove IDs associated with disabled applied voltages (if desired).
                if filter_disabled_flag, applied_voltage_IDs = self.remove_disabled_applied_voltage_IDs( applied_voltage_IDs, applied_voltages, undetected_option ); end
                
                % Retrieve the final time associated with each applied voltage.
                tf = self.get_applied_voltage_property( applied_voltage_IDs, 'tf', true, applied_voltages, undetected_option );
                
                % Determine how to process the final time.
                if strcmpi( process_option, 'average' )                     % If we want the average final time...
                    
                    % Set the step size to be the average final time.
                    tf = mean( tf );
                    
                elseif strcmpi( process_option, 'max' )                    % If we want the maximum final time...
                    
                    % Set the step size to be the largest final time.
                    tf = max( tf );
                    
                elseif strcmpi( process_option, 'min' )                    % If we want the minimum final time...
                    
                    % Set the step size to be the smallest final time.
                    tf = min( tf );
                    
                elseif strcmpi( process_option, 'none' )                   % If we have selected no process options...
                    
                    % Do nothing.
                    
                else                                                        % Otherwise...
                    
                    % Throw an error.
                    error( 'Process option %s not recognized.', process_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to retrieve the applied voltages.
        function [ ts, Vas ] = get_Vas( self, applied_voltage_IDs, dt, tf, applied_voltages, filter_disabled_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 6, filter_disabled_flag = self.filtered_disabled_flag_DEFAULT; end
            if nargin < 5, applied_voltages = self.applied_voltages; end
            if nargin < 4, tf = [  ]; end
            if nargin < 3, dt = [  ]; end
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Determine how many applied voltages to get.
            num_applied_voltages_to_get = length( applied_voltage_IDs );
            
            % Determine whether we need to set the final time.
            if isempty( tf )                                                                % If the final time is empty...
                
                % Compute the maximum final time among the given applied voltages.
                tf = self.get_tfs( applied_voltage_IDs, applied_voltages, filter_disabled_flag, 'max', undetected_option );
                
            end
            
            % Determine whether we need to set the step size.
            if isempty( dt )                                                                % If the step size is empty...
                
                % Compute the minimum step size among the given applied voltages.
                dt = self.get_dts( applied_voltage_IDs, applied_voltages, filter_disabled_flag, 'min', undetected_option );

            end
            
            % Compute the number of time steps.
            num_timesteps = floor( round( tf/dt, 8 ) ) + 1;
            
            % % Preallocate a variable to store the applied voltage properties.
            % ts = cell( num_timesteps, num_applied_voltages_to_get );
            % Vas = cell( num_timesteps, num_applied_voltages_to_get );
            
            % Preallocate a variable to store the applied voltage properties.
            ts = zeros( num_timesteps, num_applied_voltages_to_get );
            Vas = zeros( num_timesteps, num_applied_voltages_to_get );

            % Retrieve the given neuron property for each applied voltage.
            for k = 1:num_applied_voltages_to_get                           % Iterate through each of the voltages to retrieve...
                
                % Retrieve the index associated with this applied voltage ID.
                applied_voltage_index = self.get_applied_voltage_index( applied_voltage_IDs( k ), applied_voltages, undetected_option );
                
                % Determine how to retrieve this applied voltage.
                if ( applied_voltage_index >= 0 ) && ( applied_voltages( applied_voltage_index ).enabled_flag )                                                      % If the applied voltage ID is greater than or equal to zero...
                    
                    % Retrieve the applied voltages.
                    [ ts( :, k ), Vas( :, k ) ] = applied_voltages( applied_voltage_index ).sample_Vas( dt, tf, applied_voltages( applied_voltage_index ).ts, applied_voltages( applied_voltage_index ).Vas );
                    
                    % % Retrieve the applied voltages.
                    % [ ts{ k }, Vas{ k } ] = applied_voltages( applied_voltage_index ).sample_Vas( dt, tf, applied_voltages( applied_voltage_index ).ts, applied_voltages( applied_voltage_index ).Vas );
                    
                elseif ( applied_voltage_index == -1 ) || ( ~applied_voltages( applied_voltage_index ).enabled_flag )                                                % If the applied voltage ID is negative one...
                    
                    % % Set the applied voltage to zero.
                    % ts( :, k ) = num2cell( zeros( num_timesteps, 1 ) );
                    % Vas( :, k ) = num2cell( zeros( num_timesteps, 1 ) );
                    
                    % Set the applied voltage to NaN.
                    ts( :, k ) = NaN( num_timesteps, 1 );
                    Vas( :, k ) = NaN( num_timesteps, 1 );

                else                                                                                    % Otherwise...
                    
                    % Throw an error.
                    error( 'Applied voltage ID %0.2f not recognized.', applied_voltage_IDs( k ) )
                    
                end
                
            end
            
        end
        
        
        %% Applied Voltage Index & ID Functions
        
        % Implement a function to retrieve the index associated with a given applied_voltage ID.
        function applied_voltage_index = get_applied_voltage_index( self, applied_voltage_ID, applied_voltages, undetected_option )
            
            % Set the default input argument.
            if nargin < 4, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Compute the number of applied voltages.
            n_applied_voltages = length( applied_voltages );
            
            % Set a flag variable to indicate whether a matching applied_voltage index has been found.
            match_found_flag = false;
            
            % Initialize the applied_voltage index.
            applied_voltage_index = 0;
            
            while ( applied_voltage_index < n_applied_voltages ) && ( ~match_found_flag )
                
                % Advance the applied_voltage index.
                applied_voltage_index = applied_voltage_index + 1;
                
                % Check whether this applied_voltage index is a match.
                if applied_voltages( applied_voltage_index ).ID == applied_voltage_ID                       % If this applied_voltage has the correct applied_voltage ID...
                    
                    % Set the match found flag to true.
                    match_found_flag = true;
                    
                end
                
            end
            
            % Determine whether to adjust the applied voltage index.
            if ~match_found_flag                                                       % If a match was not found...
                
                % Determine how to handle when a match is not found.
                if strcmpi( undetected_option, 'error' )                            % If the undetected option is set to 'error'...
                    
                    % Throw an error.
                    error( 'No applied voltage with ID %0.0f.', applied_voltage_ID )
                    
                elseif strcmpi( undetected_option, 'warning' )                     % If the undetected option is set to 'warning'...
                    
                    % Throw a warning.
                    warning( 'No applied voltage with ID %0.0f.', applied_voltage_ID )
                    
                    % Set the applied voltage index to negative one.
                    applied_voltage_index = -1;
                    
                elseif strcmpi( undetected_option, 'ignore' )                       % If the undetected option is set to 'ignore'...
                    
                    % Set the applied voltage index to negative one.
                    applied_voltage_index = -1;
                    
                else                                                                % Otherwise...
                    
                    % Throw an error.
                    error( 'Undetected option %s not recognized.', undetected_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to validate applied_voltage IDs.
        function applied_voltage_IDs = validate_applied_voltage_IDs( self, applied_voltage_IDs, applied_voltages )
            
            % Set the default input arguments.
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Compute the number of applied voltages.
            n_applied_voltages = length( applied_voltages );
            
            % Determine whether we want get the desired applied_voltage property from all of the applied_voltages.
            if isa( applied_voltage_IDs, 'char' )                                                      % If the applied_voltage IDs variable is a character array instead of an integer srray...
                
                % Determine whether this is a valid character array.
                if  strcmpi( applied_voltage_IDs, 'all' )                  % If the character array is either 'all' or 'All'...
                    
                    % Preallocate an array to store the applied_voltage IDs.
                    applied_voltage_IDs = zeros( 1, n_applied_voltages );
                    
                    % Retrieve the applied_voltage ID associated with each applied_voltage.
                    for k = 1:n_applied_voltages                   % Iterate through each applied_voltage...
                        
                        % Store the applied_voltage ID associated with the voltage applied_voltage.
                        applied_voltage_IDs( k ) = applied_voltages( k ).ID;
                        
                    end
                    
                else                                                                        % Otherwise...
                    
                    % Throw an error.
                    error( 'Applied voltage ID must be either an array of valid applied voltage IDs or one of the strings: ''all'' or ''All''.' )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to check if the existing applied voltage IDs are unique.
        function [ unique_flag, match_logicals ] = unique_existing_applied_voltage_IDs( self, applied_voltages )
            
            % Set the default input arguments.
            if nargin < 2, applied_voltages = self.applied_voltages; end
            
            % Compute the number of applied voltages.
            n_applied_voltages = length( applied_voltages );
            
            % Retrieve all of the existing applied voltage IDs.
            applied_voltage_IDs = self.get_all_applied_voltage_IDs( applied_voltages );
            
            % Determine whether all entries are unique.
            if length( unique( applied_voltage_IDs ) ) == n_applied_voltages                    % If all of the applied voltage IDs are unique...
                
                % Set the unique flag to true.
                unique_flag = true;
                
                % Set the logicals array to true.
                match_logicals = false( 1, n_applied_voltages );
                
            else                                                                     % Otherwise...
                
                % Set the unique flag to false.
                unique_flag = false;
                
                % Set the logicals array to true.
                match_logicals = false( 1, n_applied_voltages );
                
                % Determine which applied voltages have duplicate IDs.
                for k1 = 1:n_applied_voltages                          % Iterate through each applied voltage...
                    
                    % Initialize the loop variable.
                    k2 = 0;
                    
                    % Determine whether there is another applied voltage with the same ID.
                    while ( k2 < n_applied_voltages ) && ( ~match_logicals( k1 ) ) && ( k1 ~= ( k2 + 1 ) )                    % While we haven't checked all of the applied voltages and we haven't found a match...
                        
                        % Advance the loop variable.
                        k2 = k2 + 1;
                        
                        % Determine whether this applied voltage is a match.
                        if applied_voltages( k2 ).ID == applied_voltage_IDs( k1 )                              % If this applied voltage ID is a match...
                            
                            % Set this match logical to true.
                            match_logicals( k1 ) = true;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % Implement a function to check if a proposed applied voltage ID is unique.
        function [ unique_flag, match_logicals, match_indexes ] = unique_applied_voltage_ID( self, applied_voltage_ID, applied_voltages, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Retrieve all of the existing applied voltage IDs.
            applied_voltage_IDs = self.get_all_applied_voltage_IDs( applied_voltages );
            
            % Determine whether the given applied voltage ID is one of the existing applied voltage IDs ( if so, provide the matching logicals and indexes ).
            [ match_found_flag, match_logicals, match_indexes ] = array_utilities.is_value_in_array( applied_voltage_ID, applied_voltage_IDs );
            
            % Define the uniqueness flag.
            unique_flag = ~match_found_flag;
            
        end
        
        
        % Implement a function to generate a unique applied voltage ID.
        function applied_voltage_ID = generate_unique_applied_voltage_ID( self, applied_voltages, array_utilities )
            
            % Set the default input arguments.
            if nargin < 3, array_utilities = self.array_utilities; end
            if nargin < 2, applied_voltages = self.applied_voltages; end
            
            % Retrieve the existing applied voltage IDs.
            existing_applied_voltage_IDs = self.get_all_applied_voltage_IDs( applied_voltages );
            
            % Generate a unique applied voltage ID.
            applied_voltage_ID = array_utilities.get_lowest_natural_number( existing_applied_voltage_IDs );
            
        end
        
        
        % Implement a function to generate multiple unique applied voltage IDs.
        function applied_voltage_IDs = generate_unique_applied_voltage_IDs( self, num_IDs, applied_voltages, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Retrieve the existing applied voltage IDs.
            existing_applied_voltage_IDs = self.get_all_applied_voltage_IDs( applied_voltages );
            
            % Preallocate an array to store the newly generated applied voltage IDs.
            applied_voltage_IDs = zeros( 1, num_IDs );
            
            % Generate each of the new IDs.
            for k = 1:num_IDs                           % Iterate through each of the new IDs...
                
                % Generate a unique applied voltage ID.
                applied_voltage_IDs( k ) = array_utilities.get_lowest_natural_number( [ existing_applied_voltage_IDs, applied_voltage_IDs( 1:( k - 1 ) ) ] );
                
            end
            
        end
        
        
        % Implement a function to check whether a proposed applied voltage ID is a unique natural.
        function unique_flag_natural = unique_natural_applied_voltage_ID( self, applied_voltage_ID, applied_voltages, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Initialize the unique natural to false.
            unique_flag_natural = false;
            
            % Determine whether this applied voltage ID is unique.
            unique_flag = self.unique_applied_voltage_ID( applied_voltage_ID, applied_voltages, array_utilities );
            
            % Determine whether this applied voltage ID is a unique natural.
            if unique_flag && ( applied_voltage_ID > 0 ) && ( round( applied_voltage_ID ) == applied_voltage_ID )                     % If this applied voltage ID is a unique natural...
                
                % Set the unique natural flag to true.
                unique_flag_natural = true;
                
            end
            
        end
        
        
        % Implement a function to remove disabled applied voltage IDs.
        function applied_voltage_IDs = remove_disabled_applied_voltage_IDs( self, applied_voltage_IDs, applied_voltages, undetected_option )
            
            % Set the default input arguments.
            if nargin < 4, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Retrieve the number of applied voltage IDs.
            num_applied_voltage_IDs = length( applied_voltage_IDs );
            
            % Create an array to store the indexes to remove.
            remove_indexes = zeros( 1, num_applied_voltage_IDs );
            
            % Remove any IDs associated with disabled applied voltages.
            for k = 1:num_applied_voltage_IDs                       % Iterate through each of the applied voltage IDs...
                
                % Retrieve the index associated with this applied voltage ID.
                applied_voltage_index = self.get_applied_voltage_index( applied_voltage_IDs( k ), applied_voltages, undetected_option );
                
                % Determine whether to remove this applied voltage ID.
                if ( applied_voltage_index == -1 ) || ( ~applied_voltages( applied_voltage_index ).enabled_flag )                         % If this applied voltage index is invalid or this applied voltage is disabled...
                    
                    % Store the indexes to remove.
                    remove_indexes( k ) = k;
                    
                end
                
            end
            
            % Remove any extra zeros from the remove indexes.
            remove_indexes( remove_indexes == 0 ) = [  ];
            
            % Remove the applied voltage IDs.
            applied_voltage_IDs( remove_indexes ) = [  ];
            
        end
        
        
        %% Neuron ID Functions.
        
        % Implement a function to return the applied voltage ID associated with a given neuron ID.
        function applied_voltage_ID = to_neuron_ID2applied_voltage_ID( self, to_neuron_ID, applied_voltages, undetected_option )
            
            % NOTE: This function assumes that only one applied voltage applies to each neuron.
            
            % Set the default input argument.
            if nargin < 4, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Compute the number of applied voltages.
            n_applied_voltages = length( applied_voltages );
            
            % Initialize the applied voltage detected flag.
            b_applied_voltage_detected = false;
            
            % Initialize the loop counter.
            k = 0;
            
            % Search for the applied voltage(s) that connect the specified neurons.
            while ( ~b_applied_voltage_detected ) && ( k < n_applied_voltages )              % While a matching applied voltage has not yet been detected and we haven't looked through all of the applied voltages...
                
                % Advance the loop counter.
                k = k + 1;
                
                % Determine whether this applied voltage connects the specified neurons.
                if ( applied_voltages( k ).neuron_ID == to_neuron_ID )
                    
                    % Set the applied voltage detected flag to true.
                    b_applied_voltage_detected = true;
                    
                end
                
            end
            
            % Determine whether a matching applied voltage was detected.
            if b_applied_voltage_detected                                   % If we found a matching applied voltage....
                
                % Retrieve the ID of the matching applied voltage.
                applied_voltage_ID = applied_voltages( k ).ID;
                
            else                                                    % Otherwise...
                
                % Determine how to handle the situation where we can not find a applied voltage that connects the selected neurons.
                if strcmpi( undetected_option, 'error' )                                    % If the error option is selected...
                    
                    % Throw an error.
                    error( 'No applied voltage found that stimulates neuron %0.0f.', to_neuron_ID )
                    
                elseif strcmpi( undetected_option, 'warning' )                              % If the warning option is selected...
                    
                    % Throw a warning.
                    warning( 'No applied voltage found that stimulates neuron %0.0f.', to_neuron_ID )
                    
                    % Set the synapse ID to be negative one.
                    applied_voltage_ID = -1;
                    
                elseif strcmpi( undetected_option, 'ignore' )                               % If the ignore option is selected...
                    
                    % Set the applied voltage ID to be negative one.
                    applied_voltage_ID = -1;
                    
                else                                                                        % Otherwise...
                    
                    % Throw an error.
                    error( 'undetected_option %s unrecognized.', undetected_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function return the applied voltage IDs associated with given neuron IDs.
        function applied_voltage_IDs = neuron_IDs2applied_voltage_IDs( self, to_neuron_IDs, applied_voltages, undetected_option )
            
            % Set the default input argument.
            if nargin < 4, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_currents_DEFAULT; end
            
            % Retrieve the number of applied voltages to find.
            num_applied_voltages_to_find = length( to_neuron_IDs );
            
            % Preallocate an array to store the applied voltage IDs.
            applied_voltage_IDs = zeros( 1, num_applied_voltages_to_find );
            
            % Search for each applied voltage ID.
            for k = 1:num_applied_voltages_to_find                              % Iterate through each set of neurons for which we are searching for a connecting applied voltage...
                
                % Retrieve the ID of the applied voltage that connects to this neuron.
                applied_voltage_IDs( k ) = self.to_neuron_ID2applied_voltage_ID( to_neuron_IDs( k ), applied_voltages, undetected_option );
                
            end
            
        end
        
        
        % Implement a function to return the applied voltages associated with given neuron IDs.
        function [ ts, Vas ] = to_neuron_IDs2Vas( self, to_neuron_IDs, dt, tf, applied_voltages, filter_disabled_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 6, filter_disabled_flag = self.filtered_disabled_flag_DEFAULT; end
            if nargin < 5, applied_voltages = self.applied_voltages; end
            if nargin < 4, tf = [  ]; end
            if nargin < 3, dt = [  ]; end
            
            % Retrieve the applied voltage IDs.
            applied_voltage_IDs = self.neuron_IDs2applied_voltage_IDs( to_neuron_IDs, applied_voltages, undetected_option );
            
            % Retrieve the applied voltages.
            [ ts, Vas ] = self.get_Vas( applied_voltage_IDs, dt, tf, applied_voltages, filter_disabled_flag, undetected_option );
            
        end
        
        
        %% Enable & Disable Functions.
        
        % Implement a function to enable an applied voltage.
        function [ enabled_flag, applied_voltages, self ] = enable_applied_voltage( self, applied_voltage_ID, applied_voltages, set_flag, undetected_option )
        
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Retrieve the index associated with this applied voltage.
            applied_voltage_index = self.get_applied_current_index( applied_voltage_ID, applied_voltages, undetected_option );
            
            % Enable this applied voltage.
            [ enabled_flag, applied_voltages( applied_voltage_index ) ] = applied_voltages( applied_voltage_index ).enable( true );
            
            % Determine whether to update the applied voltage manager object.
            if set_flag, self.applied_voltages = applied_voltages; end
            
        end
        
        
        % Implement a function to enable applied voltages.
        function [ enabled_flags, applied_voltages, self ] = enable_applied_voltages( self, applied_voltage_IDs, applied_voltages, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_voltages = self.applied_voltages; end                   	% [class] Array of Applied Voltage Class Objects.
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Determine the number of applied voltages to enable.
            num_applied_voltages_to_enable = length( applied_voltage_IDs );
            
            % Preallocate an array to store the enabled flags.
            enabled_flags = false( 1, num_applied_voltages_to_enable );
            
            % Enable all of the specified applied voltages.
            for k = 1:num_applied_voltages_to_enable                      % Iterate through all of the specified applied voltages...
                
                % Enable this applied voltage.
                [ enabled_flags( k ), applied_voltages, self ] = self.enable_applied_voltage( applied_voltage_IDs( k ), applied_voltages, set_flag, undetected_option );
                
            end
            
        end
        
        
        % Implement a function to disable an applied voltage.
        function [ enabled_flag, applied_voltages, self ] = disable_applied_voltage( self, applied_voltage_ID, applied_voltages, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_voltages = self.applied_voltages; end                 	% [class] Array of Neuron Class Objects.
            
            % Retrieve the index associated with this applied voltage.
            applied_current_index = self.get_applied_voltage_index( applied_voltage_ID, applied_voltages, undetected_option );
            
            % Disable this applied voltage.
            [ enabled_flag, applied_voltages( applied_current_index ) ] = applied_voltages( applied_current_index ).disable( true );
            
            % Determine whether to update the applied voltage manager object.
            if set_flag, self.applied_voltages = applied_voltages; end
            
        end
        
        
        % Implement a function to disable applied voltages.
        function [ enabled_flags, applied_voltages, self ] = disable_applied_voltages( self, applied_voltage_IDs, applied_voltages, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undected_option_DEFAULT; end            % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_voltages = self.applied_voltages; end                  	% [class] Array of Neuron Class Objects.
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Determine the number of applied voltages to disable.
            num_applied_voltages_to_enable = length( applied_voltage_IDs );
            
            % Preallocate an array to store the enabled flags.
            enabled_flags = false( 1, num_applied_voltages_to_enable );
            
            % Disable all of the specified applied voltages.
            for k = 1:num_applied_voltages_to_enable                      % Iterate through all of the specified applied voltages...
                
                % Disable this applied voltage.
                [ enabled_flags( k ), applied_voltages, self ] = self.disable_applied_voltage( applied_voltage_ID, applied_voltages, set_flag, undetected_option );
                
            end
            
        end
        
        
        % Implement a function to toggle an applied voltage's enabled flag.
        function [ enabled_flag, applied_voltages, self ] = toggle_enabled_applied_voltage( self, applied_voltage_ID, applied_voltages, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_voltages = self.applied_voltages; end                   	% [class] Array of Neuron Class Objects.
            
            % Retrieve the index associated with this applied voltage.
            applied_voltage_index = self.get_applied_current_index( applied_voltage_ID, applied_voltages, undetected_option );
            
            % Toggle whether this applied current is enabled.
            [ enabled_flag, applied_voltages( applied_voltage_index ) ] = applied_voltages( applied_voltage_index ).toggle_enabled( applied_voltages( applied_voltage_index ).enabled_flag, true );
            
            % Determine whether to update the applied current manager object.
            if set_flag, self.applied_voltages = applied_voltages; end
            
        end
        
        
        % Implement a function to toggle applied voltage enable state.
        function [ enabled_flags, applied_voltages, self ] = toggle_enabled_applied_voltages( self, applied_voltage_IDs, applied_voltages, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_voltages = self.applied_voltages; end                    % [class] Array of Applied Voltage Class Objects.
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Determine the number of applied voltages to disable.
            num_applied_voltages_to_enable = length( applied_voltage_IDs );
            
            % Preallocate an array to store the enabled flags.
            enabled_flags = false( 1, num_applied_voltages_to_enable );
            
            % Disable all of the specified applied voltages.
            for k = 1:num_applied_voltages_to_enable                      % Iterate through all of the specified applied voltages...
                
                % Toggle this applied voltage.
                [ enabled_flags( k ), applied_voltages, self ] = self.toggle_enabled_applied_voltage( applied_voltage_IDs( k ), applied_voltages, set_flag, undetected_option );
                
            end
            
        end
        
        
        %% Validation Functions.
        
        % Ensure that each neuron has only one applied voltage.
        function one_to_one_flag = one_to_one_applied_voltages( self, applied_voltages, array_utilities )
            
            % Set the default input arguments.
            if nargin < 3, array_utilities = self.array_utilites; end
            if nargin < 2, applied_voltages = self.applied_voltages; end
            
            % Compute the number of applied voltages.
            n_applied_voltages = length( applied_voltages );
            
            % Set the one-to-one flag.
            one_to_one_flag = true;
            
            % Initialize a counter variable.
            k = 0;
            
            % Preallocate arrays to store the neuron IDs.
            to_neuron_IDs = zeros( 1, n_applied_voltages );
            enabled_flags = false( 1, n_applied_voltages );
            
            % Determine whether there is only one synapse between each neuron.
            while ( one_to_one_flag ) && ( k < n_applied_voltages )                             % While we haven't found an applied voltage repetition and we haven't checked all of the applied voltages...
                
                % Advance the loop counter.
                k = k + 1;
                
                % Store these from neuron and to neuron IDs.
                to_neuron_IDs( k ) = applied_voltages( k ).neuron_ID;
                
                % Determine whether we need to check this synapse for repetition.
                if k ~= 1                               % If this is not the first iteration...
                    
                    % Determine whether this neuron ID is unique.
                    [ neuron_ID_match, neuron_ID_match_logicals ] = array_utilities.is_value_in_array( to_neuron_IDs( k ), to_neuron_IDs( 1:( k  - 1 ) ) );
                    
                    % Determine whether this applied voltage is a duplicate.
                    if neuron_ID_match && enabled_flags( k ) && any( neuron_ID_match_logicals & enabled_flags( 1:( k  - 1 ) ) )                           % If this neuron ID is a match, this applied voltage is enabled, and the matching applied voltage is enabled...
                        
                        % Set the one-to-one flag to false (this applied voltage is a duplicate).
                        one_to_one_flag = false;
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % Implement a function to validate the compatibiltiy of applied voltage properties.
        function valid_flag = validate_applied_voltage_properties( self, n_applied_voltages, IDs, names, to_neuron_IDs, ts, Vas, enabled_flags, applied_voltages, array_utilities )
            
            % Set the default input arguments.
            if nargin < 10, array_utilities = self.array_utilities; end
            if nargin < 9, applied_voltages = self.applied_voltages; end
            if nargin < 8, enabled_flags = self.enabled_flag_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 7, Vas = self.Vas_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_IDs_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 4, names = repmat( { '' }, 1, n_applied_voltages ); end
            if nargin < 3, IDs = self.generate_unique_applied_current_IDs( n_applied_voltages, applied_voltages, array_utilities ); end
            
            % Determine whether to convert the names property to a cell.
            if ~iscell( names ), names = { names }; end
            
            % Determine whether the applied current properties are relevant.
            valid_flag = ( n_applied_voltages == length( IDs ) ) && ( n_applied_voltages == length( names ) &&  n_applied_voltages == length( to_neuron_IDs ) ) && ( n_applied_voltages == length( ts ) ) && ( n_applied_voltages == length( Vas ) ) && ( n_applied_voltages == length( enabled_flags ) );
            
        end
        
        
        %% Process Functions.
        
        % Implement a function to process applied voltage creation inputs.
        function [ n_applied_voltages, IDs, names, to_neuron_IDs, ts, Vas, enabled_flags ] = process_applied_voltage_creation_inputs( self, n_applied_voltages, IDs, names, to_neuron_IDs, ts, Vas, enabled_flags, applied_voltages, array_utilities )
        
            % Set the default input arguments.
            if nargin < 10, array_utilities = self.array_utilities; end
            if nargin < 9, applied_voltages = self.applied_voltages; end
            if nargin < 8, enabled_flags = self.enabled_flag_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 7, Vas = self.Vas_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_IDs_DEFAULT*ones( 1, n_applied_voltages ); end
            if nargin < 4, names = repmat( { '' }, 1, n_applied_voltages ); end
            if nargin < 3, IDs = self.generate_unique_applied_current_IDs( n_applied_voltages, applied_voltages, array_utilities ); end

            % Convert the applied voltage parameters from cells to arrays as appropriate.
            enabled_flags = array_utilities.cell2array( enabled_flags );
            Vas = array_utilities.cell2array( Vas );
            ts = array_utilities.cell2array( ts );
            to_neuron_IDs = array_utilities.cell2array( to_neuron_IDs );
            names = array_utilities.cell2array( names );
            IDs = array_utilities.cell2array( IDs );
            n_applied_voltages = array_utilities.cell2array( n_applied_voltages );
            
            % Ensure that the applied voltage properties match the required number of applied voltages.
            assert( self.validate_applied_voltage_properties( n_applied_voltages, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities ), 'Provided applied voltage properties  must be of consistent size.' )
            
        end

        
        % Implement a function to process the applied voltage creation outputs.
        function [ IDs, applied_voltages ] = process_applied_voltage_creation_outputs( self, IDs, applied_voltages, as_cell_flag, array_utilities )
            
           	% Set the default input arguments.
            if nargin < 5, array_utilities = self.array_utilities; end                      % [class] Array Utilities Class.
            if nargin < 4, as_cell_flag = self.as_cell_flag_DEFAULT; end                  	% [T/F] As Cell Flag (Determines whether neurons are returned in an array or a cell.)
            if nargin < 3, applied_voltages = self.applied_voltages; end                  	% [class] Array of Applied Voltage Class Objects.
            
            % Determine whether to embed the new applied voltage IDs and objects in cells.
            if as_cell_flag                                                                 % If we want to embed the new applied voltage IDs and objects into cells...
                
                % Determine whether to embed the applied voltage IDs into a cell.
                if ~iscell( IDs )                                                           % If the IDs are not already a cell...
                
                    % Embed applied voltage IDs into a cell.
                    IDs = { IDs };
                
                end
                
                % Determine whether to embed the applied voltage objects into a cell.
                if ~iscell( applied_voltages )                                                       % If the applied voltages are not already a cell...
                
                    % Embed applied voltage objects into a cell.
                    applied_voltages = { applied_voltages };
                    
                end
                
            else                                                                            % Otherwise...
                
                % Determine whether to embed the applied voltage IDs into an array.
                if iscell( IDs )                                                            % If the applied voltage IDs are a cell...
                
                    % Convert the applied voltage IDs cell to a regular array.
                    IDs = array_utilities.cell2array( IDs );
                    
                end
                
                % Determine whether to embed the applied voltage objects into an array.
                if iscell( applied_voltages )                                                        % If the applied voltage objects are a cell...
                
                    % Convert the applied voltage objects cell to a regular array.
                    applied_voltages = array_utilities.cell2array( applied_voltages );
                    
                end
                
            end
            
        end
        
        
        % Implement a function to update the applied voltage manager.
        function [ applied_voltages, self ] = update_applied_voltage_manager( self, applied_voltages, applied_voltage_manager, set_flag )
        
            % Set the default input arguments.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end            % [T/F] Set Flag (Determines whether output self object is updated.)
            
            % Determine whether to update the applied voltage manager object.
            if set_flag                                                  	% If we want to update the applied voltage manager object...
                
                % Update the applied voltage manager object.
                self = applied_voltage_manager;
            
            else                                                            % Otherwise...
                
                % Reset the applied voltages object.
                applied_voltages = self.applied_voltages;
            
            end
            
        end
        
            
        %% Applied Voltage Creation Functions.
        
        % Implement a function to create a new applied voltage.
        function [ ID_new, applied_voltage_new, applied_voltages, self ] = create_applied_voltage( self, ID, name, to_neuron_ID, ts, Vas, enabled_flag, applied_voltages, set_flag, as_cell_flag, array_utilities )
            
            % Set the default input arguments.
            if nargin < 11, array_utilities = self.array_utilities; end
            if nargin < 10, as_cell_flag = self.as_cell_flag_DEFAULT; end
            if nargin < 9, set_flag = self.set_flag_DEFAULT; end
            if nargin < 8, applied_voltages = self.applied_voltages; end 
            if nargin < 7, enabled_flag = true; end
            if nargin < 6, Vas = { [ ] }; end
            if nargin < 5, ts = 0; end
            if nargin < 4, to_neuron_ID = 0; end
            if nargin < 3, name = ''; end
            if nargin < 2, ID = self.generate_unique_applied_voltage_ID( applied_voltages, array_utilities ); end
            
            % Process the applied voltage creation properties.
            [ ~, ID, name, to_neuron_ID, ts, Vas, enabled_flag ] = self.process_applied_voltage_creation_inputs( 1, ID, name, to_neuron_ID, ts, Vas, enabled_flag, applied_voltages, array_utilities );
            
            % Ensure that this neuron ID is a unique natural.
            assert( self.unique_natural_applied_voltage_ID( ID, applied_voltages, array_utilities ), 'Proposed applied voltage ID %0.2f is not a unique natural number.', ID )
            
            % Create an instance of the applied voltage manager.
            applied_voltage_manager = self;
            
            % Create an instance of the applied voltage class.
            applied_voltage_new = applied_voltage_class( ID, name, to_neuron_ID, ts, Vas, enabled_flag, array_utilities );
            
            % Retrieve the new applied voltage ID.
            ID_new = applied_voltage_new.ID;
            
            % Determine whether to embed the new applied voltage ID and object in cells.
            [ ID_new, applied_voltage_new ] = self.process_applied_voltage_creation_outputs( ID_new, applied_voltage_new, as_cell_flag, array_utilities );
            
            % Append this applied voltage to the array of existing applied voltages.
            applied_voltages = [ applied_voltages, applied_voltage_new ];
            
            % Update the applied voltage manager to reflect the update applied currents object.
            applied_voltage_manager.applied_voltages = applied_voltages;
            applied_voltage_manager.num_applied_voltages = length( applied_voltages );
            
            % Determine whether to update the applied voltage manager object.
            [ applied_voltages, self ] = self.update_applied_current_manager( applied_voltages, applied_voltage_manager, set_flag );
            
        end
        
        
        % Implement a function to create multiple applied voltages.
        function [ IDs_new, applied_voltages_new, applied_voltages, self ] = create_applied_voltages( self, n_applied_voltages_to_create, IDs, names, to_neuron_IDs, ts, Vas, enabled_flags, applied_voltages, set_flag, as_cell_flag, array_utilities )
            
            % Set the default input arguments.
            if nargin < 7, enabled_flags = self.enabled_flag*ones( 1, n_applied_voltages_to_create ); end
            if nargin < 6, Vas = self.Vas_DEFAULT*ones( 1, n_applied_voltages_to_create ); end
            if nargin < 5, ts = self.ts_DEFAULT*ones( 1, n_applied_voltages_to_create ); end
            if nargin < 4, to_neuron_IDs = self.to_neuron_ID_DEFAULT*ones( 1, n_applied_voltages_to_create ); end
            if nargin < 3, names = repmat( { '' }, 1, n_applied_voltages_to_create ); end
            if nargin < 2, IDs = self.generate_unique_applied_current_IDs( n_applied_voltages_to_create, applied_voltages, array_utilities ); end
            
            % Process the applied voltage creation inputs.
            [ n_applied_voltages, IDs, names, to_neuron_IDs, ts, Vas, enabled_flags ] = self.process_applied_current_creation_inputs( n_applied_voltages_to_create, IDs, names, to_neuron_IDs, ts, Vas, enabled_flags, applied_voltages, array_utilities );
            
            % Preallocate an array to store the new applied voltages.
            applied_voltages_new = repmat( applied_voltage_class(  ), [ 1, n_applied_voltages_to_create ] );
            
            % Preallocate an array to store the new applied voltage IDs.
            IDs_new = zeros( 1, n_applied_voltages );
            
            % Create an instance of the applied voltage manager that can be updated.
            applied_voltage_manager = self;
            
            % Create each of the specified applied voltages.
            for k = 1:n_applied_voltages                                                                                           % Iterate through each of the applied voltages we want to create...
                
                % Create this applied voltage.
                [ IDs_new( k ), applied_voltages_new( k ), applied_voltages, applied_voltage_manager ] = applied_voltage_manager.create_applied_voltage( IDs( k ), names { k }, to_neuron_IDs( k ), ts( :, k ), Vas( :, k ), enabled_flags( k ), applied_voltages, true, false, array_utilities );
                
            end
            
            % Determine whether to embed the new applied voltage ID and object in cells.
            [ IDs_new, applied_voltages_new ] = self.process_applied_voltage_creation_outputs( IDs_new, applied_voltages_new, as_cell_flag, array_utilities );
            
            % Determine whether to update the applied voltage manager object.
            [ applied_voltages, self ] = self.update_applied_voltage_manager( applied_voltages, applied_voltage_manager, set_flag );
            
        end
        
        
        % Implement a function to delete an applied voltage.
        function [ applied_voltages, self ] = delete_applied_voltage( self, applied_voltage_ID, applied_voltages, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end
            if nargin < 3, applied_voltages = self.applied_voltages; end
            
            % Create an instance of the applied voltage manager.
            applied_voltage_manager = self;
            
            % Retrieve the index associated with this applied voltage.
            applied_voltage_index = self.get_applied_voltage_index( applied_voltage_ID, applied_voltages, undetected_option );
            
            % Remove this applied voltage from the array of applied voltages.
            applied_voltages( applied_voltage_index ) = [  ];
            
            % Update the applied voltage manager to reflect these changes.
            applied_voltage_manager.applied_voltages = applied_voltages;
            applied_voltage_manager.num_applied_voltages = length( applied_voltages );
            
           % Determine whether to update the applied voltages and applied voltage manager objects.
           [ applied_voltages, self ] = self.update_applied_voltage_manager( applied_voltages, applied_voltage_manager, set_flag );

        end
        
        
        % Implement a function to delete multiple applied voltages.
        function [ applied_voltages, self ] = delete_applied_voltages( self, applied_voltage_IDs, applied_voltages, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_voltages = self.applied_voltages; end                   	% [class] Array of Applied Voltage Class Objects.
            if nargin < 2, applied_voltage_IDs = 'all'; end
            
            % Validate the applied voltage IDs.
            applied_voltage_IDs = self.validate_applied_voltage_IDs( applied_voltage_IDs, applied_voltages );
            
            % Retrieve the number of applied voltages to delete.
            num_applied_voltages_to_delete = length( applied_voltage_IDs );
            
            % Delete each of the specified applied voltages.
            for k = 1:num_applied_voltages_to_delete                      % Iterate through each of the applied voltages we want to delete...
                
                % Delete this applied voltage.
                [ applied_voltages, self ] = self.delete_applied_voltage( applied_voltage_IDs( k ), applied_voltages, set_flag, undetected_option );
                
            end
            
        end
        
        
        %% Save & Load Functions
        
        % Implement a function to save applied voltage manager data as a matlab object.
        function save( self, directory, file_name )
            
            % Set the default input arguments.
            if nargin < 3, file_name = 'Applied_Voltage_Manager.mat'; end
            if nargin < 2, directory = '.'; end
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Save the neuron data.
            save( full_path, self )
            
        end
        
        
        % Implement a function to load applied voltage manager data as a matlab object.
        function [ data, self ] = load( self, directory, file_name, set_flag )
            
            % Set the default input arguments.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end
            if nargin < 3, file_name = 'Applied_Voltage_Manager.mat'; end
            if nargin < 2, directory = '.'; end
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Load the data.
            data = load( full_path );
            
            % Determine whether to update the applied voltage manager object.
            if set_flag, self = data; end
            
        end
        
        
        % Implement a function to load applied voltage data from an xlsx file.
        function [ applied_voltages, self ] = load_xlsx( self, file_name, directory, append_flag, verbose_flag, applied_voltages, set_flag, data_loader_utilities )
            
            % Set the default input arguments.
            if nargin < 8, data_loader_utilities = self.data_loader_utilities; end          % [class] Data Load Utilities Class.
            if nargin < 7, set_flag  = self.set_flag; end                                   % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 6, applied_voltages = self.applied_voltages; end                                    % [class] Array of Applied Voltages Class Objects.
            if nargin < 5, verbose_flag = true; end
            if nargin < 4, append_flag = false; end
            if nargin < 3, directory = '.'; end
            if nargin < 2, file_name = 'Applied_Voltage_Data.xlsx'; end
            
            % Determine whether to print status messages.
            if verbose_flag, fprintf( 'LOADING APPLIED VOLTAGE DATA. Please Wait...\n' ), end
            
            % Start a timer.
            tic
            
            % Load the applied voltage data.
            [ applied_voltage_IDs, applied_voltage_names, applied_voltage_to_neuron_IDs, applied_voltage_ts, applied_voltage_Ias ] = data_loader_utilities.load_applied_voltage_data( file_name, directory );
            
            % Define the number of applied voltages.
            num_applied_voltages_to_load = length( applied_voltage_IDs );
            
            % Preallocate an array of applied voltages.
            applied_voltages_to_load = repmat( applied_voltage_class(  ), 1, num_applied_voltages_to_load );
            
            % Create each applied voltage object.
            for k = 1:num_applied_voltages_to_load               % Iterate through each of the applied voltages...
                
                % Create this applied voltage.
                applied_voltages_to_load( k ) = applied_voltage_class( applied_voltage_IDs( k ), applied_voltage_names{ k }, applied_voltage_to_neuron_IDs( k ), applied_voltage_ts( :, k ), applied_voltage_Ias( :, k ) );
                
            end
            
            % Determine whether to append the applied voltages we just loaded.
            if append_flag                         % If we want to append the applied voltages we just loaded...
                
                % Append the applied voltages we just loaded to the array of existing applied voltages.
                applied_voltages = [ applied_voltages, applied_voltages_to_load ];
                
                % Update the number of applied voltages.
                n_applied_voltages = length( applied_voltages );
                
            else                                % Otherwise...
                
                % Replace the existing applied voltages with the applied currents we just loaded.
                applied_voltages = applied_voltages_to_load;
                
                % Update the number of applied voltages.
                n_applied_voltages = length( applied_voltages );
                
            end
            
            % Determine whether to update the applied voltage manager properties.
            if set_flag                                             	% If we want to update the applied voltage manager properties...
                
                % Update the applied voltages property.
                self.applied_voltages = applied_voltages;
                
                % Update the number of applied voltages.
                self.num_applied_voltages = n_applied_voltages;
                
            end
            
            % Retrieve the elapsed time.
            elapsed_time = toc;
            
            % Determine whether to print status messages.
            if verbose_flag, fprintf( 'LOADING APPLIED VOTLAGE DATA. Please Wait... Done. %0.3f [s] \n\n', elapsed_time ), end
            
        end
        
        
    end
end

