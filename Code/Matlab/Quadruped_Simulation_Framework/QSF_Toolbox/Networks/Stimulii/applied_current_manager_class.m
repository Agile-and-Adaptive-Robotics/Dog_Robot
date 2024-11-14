classdef applied_current_manager_class
    
    % This class contains properties and methods related to managing applied currents.
    
    
    %% APPLIED CURRENT MANAGER PROPERTIES
    
    % Define general class properties.
    properties
        
        applied_currents
        num_applied_currents
        
        array_utilities
        data_loader_utilities
        applied_current_utilities
        
    end
    
    
    % Define private, constant class properties.
    properties ( Access = private, Constant = true )
        
        % Define the neuron parameters.
        R_DEFAULT = 20e-3;                                                                                  	% [V] Activation Domain.
        Gm_DEFAULT = 1e-6;                                                                                   	% [S] Membrane Conductance.
        to_neuron_ID_DEFAULT = -1;                                                                            	% [#] Neuron ID.
        
        % Define subnetwork neuron quantities.
        n_mcpg_applied_currents_DEFAULT = 1;                                                                	% [#] Number of Multistate CPG Applied Currents.
        n_centering_applied_currents_DEFAULT = 1;                                                               % [#] Number of Centering Applied Currents.
        n_dc_applied_currents_DEFAULT = 1;                                                                   	% [#] Number of Double Centering Applied Currents.
        n_inversion_applied_currents_DEFAULT = 1;                                                           	% [#] Number of Inversion Applied Currents.
        n_multiplication_applied_currents_DEFAULT = 1;                                                      	% [#] Number of Multiplication Applied Currents.
        n_integration_applied_currents_DEFAULT = 2;                                                          	% [#] Number of Integration Applied Currents.
        n_vbi_applied_currents = 2;                                                                             % [#] Number of Voltage Based Integration Applied Currents.
        n_svbi_applied_currents_DEFAULT = 3;                                                                    % [#] Number of Split Voltage Based Integration Applied Currents.
        
        % Define the default applied current properties.
        ts_DEFAULT = 0;                                                                                         % [s] Applied Current Times.
        Ias_DEFAULT = 0;                                                                                        % [A] Applied Current Magnitudes.

        % Define the simulation parameters.
        dt_DEFAULT = 1e-3;                                                                                  	% [s] Simulation Time Step.
        tf_DEFAULT = 1;                                                                                        	% [s] Simulation Duration.
        
        % Define the option defaults.
        undetected_option_DEFAULT = 'error';                                                                  	% [str] Undetected Option (Either 'error', 'warning', or 'ignore'.) (Determines how to handle situations where applied current IDs are provided that do not match an existing IDs.)
        process_option_DEFAULT = 'none';                                                                      	% [str] Process Option (Either 'max', 'min', 'mean', or 'none'.) (Determines the type of post processing that should be done to applied current properties when they are retrieved.)
        
        % Define the flag defaults.
        filter_disabled_flag_DEFAULT = true;                                                                  	% [T/F] Filter Disabled Flag. (Determines whether to filter out disabled applied currents.)
        
    end
    
    
    %% APPLIED CURRENT MANAGER METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = applied_current_manager_class( applied_currents, data_loader_utilities, array_utilities )
            
            % Set the default input arguments.
            if nargin < 3, array_utilities = array_utilities_class(  ); end                         % [class] Array Utilities Class.
            if nargin < 2, data_loader_utilities = data_loader_utilities_class(  ); end             % [class] Data Loader Utilities Class.
            if nargin < 1, applied_currents = [  ]; end                                             % [class] Array of Applied Current Class Objects.
            
            % Store the utility class properties.
            self.array_utilities = array_utilities;                                                 
            self.data_loader_utilities = data_loader_utilities;
            
            % Store the applied currents property.
            self.applied_currents = applied_currents;
            
            % Compute the number of applied currents.
            self.num_applied_currents = length( applied_currents );
            
        end
        
        
        %% Applied Current Name Functions.
            
        % Implement a function to generate names for applied currents.
        function [ names, applied_currents, self ] = generate_names( self, applied_current_IDs, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end
            if nargin < 3, applied_currents = self.applied_currents; end
            if nargin < 2, applied_current_IDs = self.get_all_applied_current_IDs( applied_currents ); end
            
            % Determine how to generate the applied current names.
            if isempty( applied_currents )              	% If there are no existing applied currents...

                % Convert the applied current IDs to applied current names.
                names = self.applied_current_utilities.IDs2names( applied_current_IDs );
                
            else                                            % If there are existing applied currents...

                % Determine the number of applied currents.
                n_applied_currents = length( applied_currents );

                % Preallocate a cell to store the applied current names.
                names = cell( 1, n_applied_currents );

                % Generate names for each of the applied currents.
                for k = 1:n_applied_currents                         % Iterate through each of the applied currents...

                    % Retrieve the index associated with this applied current.
                    applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );

                   % Generate a name for this applied current.
                   [ names{ k }, applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).generate_name( applied_current_IDs( k ), true, applied_currents( applied_current_index ).applied_current_utilities );

                end

                % Determine whether to update the applied current manager object.
                if set_flag, self.applied_currents = applied_currents; end

            end
            
        end
        
        
        %% General Get & Set Applied Current Property Functions
        
        % Implement a function to retrieve the properties of specific applied currents.
        function xs = get_applied_current_property( self, applied_current_IDs, applied_current_property, as_matrix_flag, applied_currents, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 5, applied_currents = self.applied_currents; end                   	% [class] Array of Neuron Class Objects.
            if nargin < 4, as_matrix_flag = self.as_matrix_flag_DEFAULT; end             	% [T/F] As Matrix Flag. (Determines whether to return the neuron property as a matrix or as a cell.)
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_properties_to_get = length( applied_current_IDs );
            
            % Preallocate a variable to store the applied current properties.
            xs = cell( 1, num_properties_to_get );
            
            % Retrieve the given applied current property for each applied current.
            for k = 1:num_properties_to_get                                                 % Iterate through each of the properties to get...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Define the eval string.
                eval_str = sprintf( 'xs{ k } = applied_currents( %0.0f ).%s;', applied_current_index, applied_current_property );
                
                % Evaluate the given applied current property.
                eval( eval_str );
                
            end
            
            % Determine whether to convert the network properties to a matrix.
            if as_matrix_flag                                                               % If we want the applied current properties as a matrix instead of a cell...
                
                % Convert the applied current properties from a cell to a matrix.
                xs = cell2mat( xs );
                
            end
            
        end
        
        
        % Implement a function to set the properties of specific applied currents.
        function [ applied_currents, self ] = set_applied_current_property( self, applied_current_IDs, applied_current_property_values, applied_current_property, applied_currents, set_flag )
            
            % Set the default input arguments.
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                         	% [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                    % [class] Array of Applied Current Class Objects.
            
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Validate the applied current property values.
            if ~isa( applied_current_property_values, 'cell' )                              % If the applied current property values are not a cell array...
                
                % Convert the applied current property values to a cell array.
                applied_current_property_values = num2cell( applied_current_property_values );
                
            end
            
            % Retreive the number of applied current IDs.
            num_applied_current_IDs = length( applied_current_IDs );
            %             num_applied_current_IDs = size( applied_current_IDs, 2 );
            
            % Retrieve the number of applied current property values.
            num_applied_current_property_values = length( applied_current_property_values );
            
            % Ensure that the provided neuron property values have the same length as the provided applied current IDs.
            if ( num_applied_current_IDs ~= num_applied_current_property_values )           % If the number of provided applied current IDs does not match the number of provided property values...
                
                % Determine whether to agument the property values.
                if num_applied_current_property_values == 1                                 % If there is only one provided property value...
                    
                    % Augment the property value length to match the ID length.
                    %                     applied_current_property_values = applied_current_property_values*ones( 1, num_applied_current_IDs );
                    applied_current_property_values = repmat( applied_current_property_values, [ 1, num_applied_current_IDs ] );
                    
                else                                                                    	% Otherwise...
                    
                    % Throw an error.
                    error( 'The number of provided applied current propety values must match the number of provided applied current IDs, unless a single applied current property value is provided.' )
                    
                end
                
            end
            
            % Set the properties of each applied current.
            for k = 1:n_applied_currents                                                    % Iterate through each applied current...
                
                % Determine the index of the applied current property value that we want to apply to this applied current (if we want to set a property of this applied current).
                index = find( applied_currents( k ).ID == applied_current_IDs, 1 );
                
                % Determine whether to set a property of this applied current.
                if ~isempty( index )                                                        % If a matching applied current ID was detected...
                    
                    % Create an evaluation string that sets the desired applied current property.
                    eval_string = sprintf( 'applied_currents( %0.0f ).%s = applied_current_property_values{ %0.0f };', k, applied_current_property, index );
                    
                    % Evaluate the evaluation string.
                    eval( eval_string );
                    
                end
                
            end
            
            % Determine whether to update the applied currents manager object.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        %% Specific Get & Set Functions.
        
        % Implement a function to retrieve all of the neuron IDs.
        function applied_current_IDs = get_all_applied_current_IDs( self, applied_currents )
            
            % Set the default input arguments.
            if nargin < 2, applied_currents = self.applied_currents; end            % [class] Array of Applied Current Class Objects.
            
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Preallocate a variable to store the applied current IDs.
            applied_current_IDs = zeros( 1, n_applied_currents );
            
            % Retrieve the ID associated with each applied current.
            for k = 1:n_applied_currents                                            % Iterate through each of the applied currents...
                
                % Retrieve this applied current ID.
                applied_current_IDs( k ) = applied_currents( k ).ID;
                
            end
            
        end
        
        
        % Implement a function to retrieve the number of time steps of the specified applied currents.
        function num_timesteps = get_num_timesteps( self, applied_current_IDs, applied_currents, filter_disabled_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option. (Must be either 'error', 'warning', or 'ignore'.)
            if nargin < 5, process_option = self.process_option_DEFAULT; end                        % [str] Process Option. (Must be either 'max', 'min', 'mean', or 'none'.)
            if nargin < 4, filter_disabled_flag = self.filter_distabled_flag_DEFAULT; end           % [T/F] Filter Disabled Flag. (Determines whether disabled applied currents are considered.)
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Determine how to compute the number of timesteps.
            if all( applied_current_IDs == -1 )                                                     % If all of the applied current IDs are invalid...
                
                % Set the number of timesteps to zero.
                num_timesteps = 0;
                
            else                                                                                    % Otherwise...
                
                % Remove any invalid applied current IDs.
                applied_current_IDs( applied_current_IDs == -1 ) = [  ];
                
                % Remove IDs associated with disabled applied currents (if desired).
                if filter_disabled_flag, applied_current_IDs = self.remove_disabled_applied_current_IDs( applied_current_IDs, applied_currents, undetected_option ); end
                
                % Retrieve the number of timesteps associated with each applied current.
                num_timesteps = self.get_applied_current_property( applied_current_IDs, 'num_timesteps', true, applied_currents, undetected_option );
                
                % Determine how to process the number of timesteps.
                if strcmpi( process_option, 'average' )                                             % If we want the average time step...
                    
                    % Set the number of timesteps to be the average number of timesteps.
                    num_timesteps = mean( num_timesteps );
                    
                elseif strcmpi( process_option, 'max' )                                             % If we want the maximum time step...
                    
                    % Set the number of timesteps to be the largest number of timesteps.
                    num_timesteps = max( num_timesteps );
                    
                elseif strcmpi( process_option, 'min' )                                             % If we want the minimum time step...
                    
                    % Set the number of timesteps to be the smallest number of timesteps.
                    num_timesteps = min( num_timesteps );
                    
                elseif strcmpi( process_option, 'none' )                                            % If we have selected no process options...
                    
                    % Do nothing.
                    
                else                                                                                % Otherwise...
                    
                    % Throw an error.
                    error( 'Process option %s not recognized.', process_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to retrieve the step size of the specified applied currents.
        function dt = get_dts( self, applied_current_IDs, applied_currents, filter_disabled_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option. (Must be either 'error', 'warning', or 'ignore'.)
            if nargin < 5, process_option = self.process_option_DEFAULT; end                        % [str] Process Option. (Must be either 'max', 'min', 'mean', or 'none'.)
            if nargin < 4, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end            % [T/F] Filter Disabled Flag. (Determines whether to considered disabled applied currents.)
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Determine how to compute the step size.
            if all( applied_current_IDs == -1 )                                                     % If all of the applied current IDs are invalid...
                
                % Set the step size to zero.
                dt = 1e-3;
                
            else                                                                                    % Otherwise...
                
                % Remove any invalid applied current IDs.
                applied_current_IDs( applied_current_IDs == -1 ) = [  ];
                
                % Remove IDs associated with disabled applied currents (if desired).
                if filter_disabled_flag, applied_current_IDs = self.remove_disabled_applied_current_IDs( applied_current_IDs, applied_currents, undetected_option ); end
                
                % Retrieve the step size associated with each applied current.
                dt = self.get_applied_current_property( applied_current_IDs, 'dt', true, applied_currents, undetected_option );
                
                % Determine whether the step size needs to be set.
                if isempty( dt ), dt = 1e-3; end
                
                % Determine how to process the step size.
                if strcmpi( process_option, 'average' )                                             % If we want the average step size...
                    
                    % Set the step size to be the average step size.
                    dt = mean( dt );
                    
                elseif strcmpi( process_option, 'max' )                                             % If we want the maximum step size...
                    
                    % Set the step size to be the largest step size.
                    dt = max( dt );
                    
                elseif strcmpi( process_option, 'min' )                                             % If we want the minimum step size...
                    
                    % Set the step size to be the smallest step size.
                    dt = min( dt );
                    
                elseif strcmpi( process_option, 'none' )                                            % If we have selected no process options...
                    
                    % Do nothing.
                    
                else                                                                                % Otherwise...
                    
                    % Throw an error.
                    error( 'Process option %s not recognized.', process_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to retrieve the final time of the specified applied currents.
        function tf = get_tfs( self, applied_current_IDs, applied_currents, filter_disabled_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 6, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option. (Must be either 'error', 'warning', or 'ignore'.)
            if nargin < 5, process_option = self.process_option_DEFAULT; end                        % [str] Process Option. (Must be either 'max', 'min', 'mean', or 'none'.)
            if nargin < 4, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end            % [T/F] Filter Disabled Flag. (Determines whether to considered disabled applied currents.)    
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Determine how to compute the final time.
            if all( applied_current_IDs == -1 )                                                     % If all of the applied current IDs are invalid...
                
                % Set the final time to zero.
                tf = 0;
                
            else                                                                                    % Otherwise...
                
                % Remove any invalid applied current IDs.
                applied_current_IDs( applied_current_IDs == -1 ) = [  ];
                
                % Remove IDs associated with disabled applied currents (if desired).
                if filter_disabled_flag, applied_current_IDs = self.remove_disabled_applied_current_IDs( applied_current_IDs, applied_currents, undetected_option ); end
                
                % Retrieve the final time associated with each applied current.
                tf = self.get_applied_current_property( applied_current_IDs, 'tf', true, applied_currents, undetected_option );
                
                % Determine how to process the final time.
                if strcmpi( process_option, 'average' )                                             % If we want the average final time...
                    
                    % Set the step size to be the average final time.
                    tf = mean( tf );
                    
                elseif strcmpi( process_option, 'max' )                                             % If we want the maximum final time...
                    
                    % Set the step size to be the largest final time.
                    tf = max( tf );
                    
                elseif strcmpi( process_option, 'min' )                                             % If we want the minimum final time...
                    
                    % Set the step size to be the smallest final time.
                    tf = min( tf );
                    
                elseif strcmpi( process_option, 'none' )                                            % If we have selected no process options...
                    
                    % Do nothing.
                    
                else                                                                                % Otherwise...
                    
                    % Throw an error.
                    error( 'Process option %s not recognized.', process_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to retrieve the applied currents.
        function [ ts, Ias ] = get_Ias( self, applied_current_IDs, dt, tf, applied_currents, filter_disabled_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option. (Must be either 'error', 'warning', or 'ignore'.)
            if nargin < 6, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end            % [T/F] Filter Disabled Flag. (Determines whether to considered disabled applied currents.)  
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, tf = [  ]; end                                                           % [s] Final Simulation Time.
            if nargin < 3, dt = [  ]; end                                                           % [s] Simulation Timestep.
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Determine how many applied currents to get.
            num_applied_currents_to_get = length( applied_current_IDs );
            
            % Determine whether we need to set the final time.
            if isempty( tf )                                                                % If the final time is empty...
                
                % Compute the maximum final time among the given applied currents.
                tf = self.get_tfs( applied_current_IDs, applied_currents, filter_disabled_flag, 'max', undetected_option );
                
            end
            
            % Determine whether we need to set the step size.
            if isempty( dt )                                                                % If the step size is empty...
                
                % Compute the minimum step size among the given applied currents.
                dt = self.get_dts( applied_current_IDs, applied_currents, filter_disabled_flag, 'min', undetected_option );
                
            end
            
            % Compute the number of time steps.
            num_timesteps = floor( round( tf/dt, 8 ) ) + 1;
            
            % Preallocate a variable to store the applied current properties.
            ts = zeros( num_timesteps, num_applied_currents_to_get );
            Ias = zeros( num_timesteps, num_applied_currents_to_get );
            
            % Retrieve the given neuron property for each applied current.
            for k = 1:num_applied_currents_to_get                           % Iterate through each of the currents to retrieve...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );            % Undetected option hard set to 'ignore' before revisions.
                
                % Determine how to retrieve this applied current.
                if ( applied_current_index >= 0 ) && ( applied_currents( applied_current_index ).enabled_flag )                                                      % If the applied current ID is greater than or equal to zero...
                    
                    % Retrieve the applied currents.
                    [ ts( :, k ), Ias( :, k ) ] = applied_currents( applied_current_index ).sample_Ias( dt, tf, applied_currents( applied_current_index ).ts, applied_currents( applied_current_index ).Ias );
                    
                elseif ( applied_current_index == -1 ) || ( ~applied_currents( applied_current_index ).enabled_flag )                                                % If the applied current ID is negative one...
                    
                    % Set the applied current to zero.
                    ts = zeros( num_timesteps, 1 );
                    Ias( :, k ) = zeros( num_timesteps, 1 );
                    
                else                                                                                    % Otherwise...
                    
                    % Throw an error.
                    error( 'Applied current ID %0.2f not recognized.', applied_current_IDs( k ) )
                    
                end
                
            end
            
        end
        
        
        %% Applied Current Index & ID Functions.
        
        % Implement a function to retrieve the index associated with a given applied_current ID.
        function applied_current_index = get_applied_current_index( self, applied_current_ID, applied_currents, undetected_option )
            
            % Set the default input argument.
            if nargin < 4, undetected_option = 'error'; end                                         % [str] Undetected Option. (Must be either 'error', 'warning', or 'ignore'.)
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Set a flag variable to indicate whether a matching applied_current index has been found.
            match_found_flag = false;
            
            % Initialize the applied_current index.
            applied_current_index = 0;
            
            while ( applied_current_index < n_applied_currents ) && ( ~match_found_flag )                          % While we have not yet checked all of the applied currenst and have not yet found an ID match...
                
                % Advance the applied_current index.
                applied_current_index = applied_current_index + 1;
                
                % Check whether this applied_current index is a match.
                if applied_currents( applied_current_index ).ID == applied_current_ID                       % If this applied_current has the correct applied_current ID...
                    
                    % Set the match found flag to true.
                    match_found_flag = true;
                    
                end
                
            end
            
            % Determine whether to adjust the applied current index.
            if ~match_found_flag                                                       % If a match was not found...
                
                % Determine how to handle when a match is not found.
                if strcmpi( undetected_option, 'error' )                            % If the undetected option is set to 'error'...
                    
                    % Throw an error.
                    error( 'No applied current with ID %0.0f.', applied_current_ID )
                    
                elseif strcmpi( undetected_option, 'warning' )                     % If the undetected option is set to 'warning'...
                    
                    % Throw a warning.
                    warning( 'No applied current with ID %0.0f.', applied_current_ID )
                    
                    % Set the applied current index to negative one.
                    applied_current_index = -1;
                    
                elseif strcmpi( undetected_option, 'ignore' )                       % If the undetected option is set to 'ignore'...
                    
                    % Set the applied current index to negative one.
                    applied_current_index = -1;
                    
                else                                                                % Otherwise...
                    
                    % Throw an error.
                    error( 'Undetected option %s not recognized.', undetected_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to validate applied_current IDs.
        function applied_current_IDs = validate_applied_current_IDs( self, applied_current_IDs, applied_currents )
            
            % Set the default input arguments.
            if nargin < 3, applied_currents = self.applied_currents; end            % [class] Array of Applied Current Class Objects.
            
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Determine whether we want get the desired applied_current property from all of the applied_currents.
            if isa( applied_current_IDs, 'char' )                                 	% If the applied_current IDs variable is a character array instead of an integer srray...
                
                % Determine whether this is a valid character array.
                if  strcmpi( applied_current_IDs, 'all' )                           % If the character array is either 'all' or 'All'...
                    
                    % Preallocate an array to store the applied_current IDs.
                    applied_current_IDs = zeros( 1, n_applied_currents );
                    
                    % Retrieve the applied_current ID associated with each applied_current.
                    for k = 1:n_applied_currents                                  	% Iterate through each applied_current...
                        
                        % Store the applied_current ID associated with the current applied_current.
                        applied_current_IDs( k ) = applied_currents( k ).ID;
                        
                    end
                    
                else                                                             	% Otherwise...
                    
                    % Throw an error.
                    error( 'Applied current ID must be either an array of valid applied_current IDs or one of the strings: ''all'' or ''All''.' )
                    
                end
                
            end
            
        end
        
        
        % Implement a function to check if the existing applied current IDs are unique.
        function [ unique_flag, match_logicals ] = unique_existing_applied_current_IDs( self, applied_currents )
            
            % Set the default input arguments.
            if nargin < 2, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Retrieve all of the existing applied current IDs.
            applied_current_IDs = self.get_all_applied_current_IDs( applied_currents );
            
            % Determine whether all entries are unique.
            if length( unique( applied_current_IDs ) ) == n_applied_currents                    % If all of the applied current IDs are unique...
                
                % Set the unique flag to true.
                unique_flag = true;
                
                % Set the logicals array to true.
                match_logicals = false( 1, n_applied_currents );
                
            else                                                                     % Otherwise...
                
                % Set the unique flag to false.
                unique_flag = false;
                
                % Set the logicals array to true.
                match_logicals = false( 1, n_applied_currents );
                
                % Determine which applied currents have duplicate IDs.
                for k1 = 1:n_applied_currents                          % Iterate through each applied current...
                    
                    % Initialize the loop variable.
                    k2 = 0;
                    
                    % Determine whether there is another applied current with the same ID.
                    while ( k2 < n_applied_currents ) && ( ~match_logicals( k1 ) ) && ( k1 ~= ( k2 + 1 ) )                    % While we haven't checked all of the applied currents and we haven't found a match...
                        
                        % Advance the loop variable.
                        k2 = k2 + 1;
                        
                        % Determine whether this applied current is a match.
                        if applied_currents( k2 ).ID == applied_current_IDs( k1 )                              % If this applied current ID is a match...
                            
                            % Set this match logical to true.
                            match_logicals( k1 ) = true;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % Implement a function to check if a proposed applied current ID is unique.
        function [ unique_flag, match_logicals, match_indexes ] = unique_applied_current_ID( self, applied_current_ID, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end          % [class] Array Utilities Class.
            if nargin < 3, applied_currents = self.applied_currents; end        % [class] Array of Applied Current Class Objects.
            
            % Retrieve all of the existing applied current IDs.
            applied_current_IDs = self.get_all_applied_current_IDs( applied_currents );
            
            % Determine whether the given applied current ID is one of the existing applied current IDs ( if so, provide the matching logicals and indexes ).
            [ match_found_flag, match_logicals, match_indexes ] = array_utilities.is_value_in_array( applied_current_ID, applied_current_IDs );
            
            % Define the uniqueness flag.
            unique_flag = ~match_found_flag;
            
        end
        
        
        % Implement a function to generate a unique applied current ID.
        function applied_current_ID = generate_unique_applied_current_ID( self, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 3, array_utilities = self.array_utilities; end         	% [class] Array Utilities Class.
            if nargin < 2, applied_currents = self.applied_currents; end        % [class] Array of Applied Current Class Objects.
            
            % Retrieve the existing applied current IDs.
            existing_applied_current_IDs = self.get_all_applied_current_IDs( applied_currents );
            
            % Generate a unique applied current ID.
            applied_current_ID = array_utilities.get_lowest_natural_number( existing_applied_current_IDs );
            
        end
        
        
        % Implement a function to generate multiple unique applied current IDs.
        function applied_current_IDs = generate_unique_applied_current_IDs( self, num_IDs, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end          % [class] Array Utilities Class.
            if nargin < 3, applied_currents = self.applied_currents; end      	% [class] Array of Applied Current Class Objects.
            
            % Retrieve the existing applied current IDs.
            existing_applied_current_IDs = self.get_all_applied_current_IDs( applied_currents );
            
            % Preallocate an array to store the newly generated applied current IDs.
            applied_current_IDs = zeros( 1, num_IDs );
            
            % Generate each of the new IDs.
            for k = 1:num_IDs                                                   % Iterate through each of the new IDs...
                
                % Generate a unique applied current ID.
                applied_current_IDs( k ) = array_utilities.get_lowest_natural_number( [ existing_applied_current_IDs, applied_current_IDs( 1:( k - 1 ) ) ] );
                
            end
            
        end
        
        
        % Implement a function to check whether a proposed applied current ID is a unique natural.
        function unique_flag_natural = unique_natural_applied_current_ID( self, applied_current_ID, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                                                      % [class] Array Utilities Class.
            if nargin < 3, applied_currents = self.applied_currents; end                                                    % [class] Array of Applied Current Class Objects.
            
            % Initialize the unique natural to false.
            unique_flag_natural = false;
            
            % Determine whether this applied current ID is unique.
            unique_flag = self.unique_applied_current_ID( applied_current_ID, applied_currents, array_utilities );
            
            % Determine whether this applied current ID is a unique natural.
            if unique_flag && ( applied_current_ID > 0 ) && ( round( applied_current_ID ) == applied_current_ID )           % If this applied current ID is a unique natural...
                
                % Set the unique natural flag to true.
                unique_flag_natural = true;
                
            end
            
        end
        
        
        % Implement a function to remove disabled applied current IDs.
        function applied_current_IDs = remove_disabled_applied_current_IDs( self, applied_current_IDs, applied_currents, undetected_option )
            
            % Set the default input arguments.
            if nargin < 4, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 3, applied_currents = self.applied_currents; end                   	% [class] Array of Applied Current Class Objects.
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Retrieve the number of applied current IDs.
            num_applied_current_IDs = length( applied_current_IDs );
            
            % Create an array to store the indexes to remove.
            remove_indexes = zeros( 1, num_applied_current_IDs );
            
            % Remove any IDs associated with disabled applied currents.
            for k = 1:num_applied_current_IDs                       % Iterate through each of the applied current IDs...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );             % Undetected option was hard set to 'ignore' before revisions.
                
                % Determine whether to remove this applied current ID.
                if ( applied_current_index == -1 ) || ( ~applied_currents( applied_current_index ).enabled_flag )                         % If this applied current index is invalid or this applied current is disabled...
                    
                    % Store the indexes to remove.
                    remove_indexes( k ) = k;
                    
                end
                
            end
            
            % Remove any extra zeros from the remove indexes.
            remove_indexes( remove_indexes == 0 ) = [  ];
            
            % Remove the applied current IDs.
            applied_current_IDs( remove_indexes ) = [  ];
            
        end
        
        
        %% Neuron ID Functions.
        
        % Implement a function to return the applied current ID associated with a given neuron ID.
        function applied_current_ID = to_neuron_ID2applied_current_ID( self, to_neuron_ID, applied_currents, undetected_option )
            
            % NOTE: This function assumes that only one applied current applies to each neuron.
            
            % Set the default input argument.
            if nargin < 4, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 3, applied_currents = self.applied_currents; end                   	% [class] Array of Applied Current Class Objects.
            
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Initialize the applied current detected flag.
            b_applied_current_detected = false;
            
            % Initialize the loop counter.
            k = 0;
            
            % Search for the applied current(s) that connect the specified neurons.
            while ( ~b_applied_current_detected ) && ( k < n_applied_currents )              % While a matching applied current has not yet been detected and we haven't looked through all of the applied currents...
                
                % Advance the loop counter.
                k = k + 1;
                
                % Determine whether this applied current connects the specified neurons.
                if ( applied_currents( k ).to_neuron_ID == to_neuron_ID )
                    
                    % Set the applied current detected flag to true.
                    b_applied_current_detected = true;
                    
                end
                
            end
            
            % Determine whether a matching applied current was detected.
            if b_applied_current_detected                                   % If we found a matching applied current....
                
                % Retrieve the ID of the matching applied current.
                applied_current_ID = applied_currents( k ).ID;
                
            else                                                    % Otherwise...
                
                % Determine how to handle the situation where we can not find a applied current that connects the selected neurons.
                if strcmpi( undetected_option, 'error' )                                    % If the error option is selected...
                    
                    % Throw an error.
                    error( 'No applied current found that stimulates neuron %0.0f.', to_neuron_ID )
                    
                elseif strcmpi( undetected_option, 'warning' )                              % If the warning option is selected...
                    
                    % Throw a warning.
                    warning( 'No applied current found that stimulates neuron %0.0f.', to_neuron_ID )
                    
                    % Set the synapse ID to be negative one.
                    applied_current_ID = -1;
                    
                elseif strcmpi( undetected_option, 'ignore' )                               % If the ignore option is selected...
                    
                    % Set the applied current ID to be negative one.
                    applied_current_ID = -1;
                    
                else                                                                        % Otherwise...
                    
                    % Throw an error.
                    error( 'undetected_option %s unrecognized.', undetected_option )
                    
                end
                
            end
            
        end
        
        
        % Implement a function return the applied current IDs associated with given neuron IDs.
        function applied_current_IDs = to_neuron_IDs2applied_current_IDs( self, to_neuron_IDs, applied_currents, undetected_option )
            
            % Set the default input argument.
            if nargin < 4, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected option. (Must be either 'error', 'warning', or 'ignore'.)
            if nargin < 3, applied_currents = self.applied_currents; end                    % [class] Array of Applied Current Class Objects.
            
            % Retrieve the number of applied currents to find.
            num_applied_currents_to_find = length( to_neuron_IDs );
            
            % Preallocate an array to store the applied current IDs.
            applied_current_IDs = zeros( 1, num_applied_currents_to_find );
            
            % Search for each applied current ID.
            for k = 1:num_applied_currents_to_find                                          % Iterate through each set of neurons for which we are searching for a connecting applied current...
                
                % Retrieve the ID of the applied current that connects to this neuron.
                applied_current_IDs( k ) = self.to_neuron_ID2applied_current_ID( to_neuron_IDs( k ), applied_currents, undetected_option );
                
            end
            
        end
        
        
        % Implement a function to return the applied currents associated with given neuron IDs.
        function [ ts, Ias ] = to_neuron_IDs2Ias( self, to_neuron_IDs, dt, tf, applied_currents, filter_disabled_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 8, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 7, process_option = self.process_option_DEFAULT; end                        % [str] Process Option. (Must be either 'max', 'min', 'mean', or 'none'.)
            if nargin < 6, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end            % [T/F] Filter Disabled Flag. (Determines whether to considered disabled applied currents.)  
            if nargin < 5, applied_currents = self.applied_currents_DEFAULT; end                    % [class] Array of Applied Current Class Objects.
            if nargin < 4, tf = [  ]; end                                                           % [s] Final Simulation Time.
            if nargin < 3, dt = [  ]; end                                                           % [s] Simulation Timestep.
            
            % Retrieve the applied current IDs.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( to_neuron_IDs, applied_currents, undetected_option );
            
            % Retrieve the applied currents.
            [ ts, Ias ] = self.get_Ias( applied_current_IDs, dt, tf, applied_currents, filter_disabled_flag, undetected_option );
        end
        
        
        %% Enable & Disable Functions.
        
        % Implement a function to enable an applied current.
        function [ enabled_flag, applied_currents, self ] = enable_applied_current( self, applied_current_ID, applied_currents, set_flag, undetected_option )
        
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Retrieve the index associated with this applied current.
            applied_current_index = self.get_applied_current_index( applied_current_ID, applied_currents, undetected_option );
            
            % Enable this applied current.
            [ enabled_flag, applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).enable( true );
            
            % Determine whether to update the applied current manager object.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to enable applied currents.
        function [ enabled_flags, applied_currents, self ] = enable_applied_currents( self, applied_current_IDs, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_currents = self.applied_currents; end                   	% [class] Array of Applied Current Class Objects.
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Determine the number of applied currents to enable.
            num_applied_currents_to_enable = length( applied_current_IDs );
            
            % Preallocate an array to store the enabled flags.
            enabled_flags = false( 1, num_applied_currents_to_enable );
            
            % Enable all of the specified applied currents.
            for k = 1:num_applied_currents_to_enable                      % Iterate through all of the specified applied currents...
                
                % Enable this applied current.
                [ enabled_flags( k ), applied_currents, self ] = self.enable_applied_current( applied_current_IDs( k ), applied_currents, set_flag, undetected_option );
                
            end
            
        end
        
        
        % Implement a function to disable an applied current.
        function [ enabled_flag, applied_currents, self ] = disable_applied_current( self, applied_current_ID, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_currents = self.applied_currents; end                 	% [class] Array of Neuron Class Objects.
            
            % Retrieve the index associated with this applied current.
            applied_current_index = self.get_applied_current_index( applied_current_ID, applied_currents, undetected_option );
            
            % Disable this applied current.
            [ enabled_flag, applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).disable( true );
            
            % Determine whether to update the applied current manager object.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to disable applied currents.
        function [ enabled_flags, applied_currents, self ] = disable_applied_currents( self, applied_current_IDs, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undected_option_DEFAULT; end            % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_currents = self.applied_currents; end                   	% [class] Array of Neuron Class Objects.
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_currents_IDs( applied_current_IDs, applied_currents );
            
            % Determine the number of applied currents to disable.
            num_applied_currents_IDs = length( applied_current_IDs );
            
            % Preallocate an array to store the enabled flags.
            enabled_flags = false( 1, num_applied_currents_IDs );
            
            % Disable all of the specified applied currents.
            for k = 1:num_applied_current_IDs                                               % Iterate through all of the specified applied currents...
                
                % Disable this applied current.
                [ enabled_flags( k ), applied_currents, self ] = self.disable_neuron( applied_current_IDs( k ), applied_currents, set_flag, undetected_option );
                
            end
            
        end
        
        
        % Implement a function to toggle an applied current's enabled flag.
        function [ enabled_flag, applied_currents, self ] = toggle_enabled_applied_current( self, applied_current_ID, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_currents = self.applied_currents; end                   	% [class] Array of Neuron Class Objects.
            
            % Retrieve the index associated with this applied current.
            applied_current_index = self.get_applied_current_index( applied_current_ID, applied_currents, undetected_option );
            
            % Toggle whether this applied current is enabled.
            [ enabled_flag, applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).toggle_enabled( applied_currents( applied_current_index ).enabled_flag, true );
            
            % Determine whether to update the applied current manager object.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to toggle multiple applied current enable states.
        function [ enabled_flags, applied_currents, self ] = toggle_enabled_applied_currents( self, applied_current_IDs, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_currents = self.applied_currents; end                    % [class] Array of Applied Current Class Objects.
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Determine the number of applied_currents to disable.
            num_applied_current_IDs = length( applied_current_IDs );
            
            % Preallocate an array to store the enabled flags.
            enabled_flags = false( 1, num_applied_current_IDs );
            
            % Disable all of the specified applied currents.
            for k = 1:num_applied_current_IDs                                                        % Iterate through all of the specified applied currents...
                
                % Toggle this applied current.
                [ enabled_flags( k ), applied_currents, self ] = self.toggle_enabled_applied_current( applied_current_IDs( k ), applied_currents, set_flag, undetected_option );
                
            end
            
        end
        
        
        %% Validation Functions.
        
        % Ensure that each neuron has only one applied current.
        function one_to_one_flag = one_to_one_applied_currents( self, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 3, array_utilities = self.array_utilities; end                                                      % [class] Array Utilities Class.
            if nargin < 2, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Set the one-to-one flag.
            one_to_one_flag = true;
            
            % Initialize a counter variable.
            k = 0;
            
            % Preallocate arrays to store the neuron IDs.
            to_neuron_IDs = zeros( 1, n_applied_currents );
            enabled_flags = false( 1, n_applied_currents );
            
            % Determine whether there is only one synapse between each neuron.
            while ( one_to_one_flag ) && ( k < n_applied_currents )                             % While we haven't found an applied current repetition and we haven't checked all of the applied currents...
                
                % Advance the loop counter.
                k = k + 1;
                
                % Store these from neuron and to neuron IDs.
                to_neuron_IDs( k ) = applied_currents( k ).to_neuron_ID;
                
                % Determine whether we need to check this synapse for repetition.
                if k ~= 1                               % If this is not the first iteration...
                    
                    % Determine whether this neuron ID is unique.
                    [ to_neuron_ID_match, to_neuron_ID_match_logicals ] = array_utilities.is_value_in_array( to_neuron_IDs( k ), to_neuron_IDs( 1:( k  - 1 ) ) );
                    
                    % Determine whether this applied current is a duplicate.
                    if to_neuron_ID_match && enabled_flags( k ) && any( to_neuron_ID_match_logicals & enabled_flags( 1:( k  - 1 ) ) )                           % If this neuron ID is a match, this applied current is enabled, and the matching applied current is enabled...
                        
                        % Set the one-to-one flag to false (this applied current is a duplicate).
                        one_to_one_flag = false;
                        
                    end
                    
                end
                
            end
            
        end
        
        
        %% Applied Current Magnitude Parameter Processing Functions.
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to process inversion Ias2 parameters.
        function parameters = process_inversion_Ias2_parameters( self, parameters, encoding_scheme )
        
            % Set the default input arguments.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm2 = self.Gm_DEFAULT;
                    R2 = self.R_DEFAULT;                           
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm2, R2 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm2 = self.Gm_DEFAULT;
                    R2 = self.R_DEFAULT;                           
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm2, R2 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end

        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to process reduced inversion Ias2 parameters.
        function parameters = process_reduced_inversion_Ias2_parameters( self, parameters, encoding_scheme )
        
            % Set the default input arguments.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm2 = self.Gm_DEFAULT;
                    R2 = self.R_DEFAULT;                           
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm2, R2 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm2 = self.Gm_DEFAULT;
                    R2 = self.R_DEFAULT;                           
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm2, R2 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end

        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to process multiplication Ias3 parameters.
        function parameters = process_multiplication_Ias3_parameters( self, parameters, encoding_scheme, applied_currents )
        
            % Set the default input arguments.
            if nargin < 4, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm3 = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R3 = self.R_DEFAULT*ones( 1, n_applied_currents );
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm3, R3 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm3 = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R3 = self.R_DEFAULT*ones( 1, n_applied_currents );                        
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm3, R3 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to process reduced multiplication Ias3 parameters.
        function parameters = process_reduced_multiplication_Ias3_parameters( self, parameters, encoding_scheme, applied_currents )
        
            % Set the default input arguments.
            if nargin < 4, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm3 = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R3 = self.R_DEFAULT*ones( 1, n_applied_currents );
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm3, R3 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm3 = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R3 = self.R_DEFAULT*ones( 1, n_applied_currents );                        
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm3, R3 };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------

        % Implement a function to process integration Ias parameters.
        function parameters = process_integration_Ias_parameters( self, parameters, encoding_scheme, applied_currents )
        
            % Set the default input arguments.
            if nargin < 4, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );                        
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end

        
        % Implement a function to process voltage based integration Ias parameters.
        function parameters = process_vbi_Ias_parameters( self, parameters, encoding_scheme, applied_currents )
        
            % Set the default input arguments.
            if nargin < 4, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );                        
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % Implement a function to process split voltage based integration Ias parameters.
        function parameters = process_svbi_Ias1_parameters( self, parameters, encoding_scheme, applied_currents )
        
            % Set the default input arguments.
            if nargin < 4, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );                        
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % Implement a function to process split voltage based integration Ias parameters.
        function parameters = process_svbi_Ias2_parameters( self, parameters, encoding_scheme, applied_currents )
        
            % Set the default input arguments.
            if nargin < 4, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 2, parameters = {  }; end
           
            % Compute the number of applied currents.
            n_applied_currents = length( applied_currents );
            
            % Determine how to create the parameters cell.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If this operation is using an absolute encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using an absolute encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If this operation uses a relative encoding scheme...
                
                % Determine how to create the parameters cell given that this operation is using a relative encoding scheme.
                if isempty( parameters )                                                % If no parameters were provided...
                    
                    % Set the default input and output voltage offsets.
                    Gm = self.Gm_DEFAULT*ones( 1, n_applied_currents );
                    R = self.R_DEFAULT*ones( 1, n_applied_currents );                        
                    
                    % Store the required parameters in a cell.
                    parameters = { Gm, R };
                    
                else                                                                    % Otherwise...
                    
                    % Determine whether the parameters cell has a valid number of entries.
                    if length( parameters ) ~= 2                                        % If there is anything other than three parameter entries...
                        
                        % Throw an error.
                        error( 'Invalid parameters detected.' )
                        
                    end
                    
                end
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        %% Compute Functions.
    
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the transmission subnetwork output applied currents.
        function Ias2 = compute_transmission_Ias2( self, encoding_scheme, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Define the dummy applied current properties.
            ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities );      % [#] Dummy Applied Current ID.
            name = 'Dummy Current';                                                                 % [str] Dummy Applied Curret Name.
            to_neuron_ID = 0;                                                                       % [#] Dummy Neuron ID to Which Current is Applied.
            ts = 0;                                                                                 % [s] Dummy Current Time Vector.
            Ias = 0;                                                                                % [A] Dummy Applied Current Magnitudes.
            enabled_flag = false;                                                                   % [F] Dummy Enabled Flags.
                        
            % Create a dummy applied current to calculate the applied current magnitude.
            [ ~, applied_current, ~, ~ ] = self.create_applied_current( ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, false, false, array_utilities );
            
            % Compute the magnitude for the output applied current.                        
            [ Ias2, ~ ] = applied_current.compute_transmission_Ias2( encoding_scheme, false, applied_current.applied_current_utilities );
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the addition subnetwork output applied currents.
        function Ias_n = compute_addition_Iasn( self, encoding_scheme, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Define the dummy applied current properties.
            ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities );      % [#] Dummy Applied Current ID.
            name = 'Dummy Current';                                                                 % [str] Dummy Applied Curret Name.
            to_neuron_ID = 0;                                                                       % [#] Dummy Neuron ID to Which Current is Applied.
            ts = 0;                                                                                 % [s] Dummy Current Time Vector.
            Ias = 0;                                                                                % [A] Dummy Applied Current Magnitudes.
            enabled_flag = false;                                                                   % [F] Dummy Enabled Flags.
                        
            % Create a dummy applied current to calculate the applied current magnitude.
            [ ~, applied_current, ~, ~ ] = self.create_applied_current( ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, false, false, array_utilities );
            
            % Compute the magnitude for the output applied current.                        
            [ Ias_n, ~ ] = applied_current.compute_addition_Iasn( encoding_scheme, false, applied_current.applied_current_utilities );
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the subtraction subnetwork output applied currents.
        function Ias_n = compute_subtraction_Iasn( self, encoding_scheme, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Define the dummy applied current properties.
            ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities );      % [#] Dummy Applied Current ID.
            name = 'Dummy Current';                                                                 % [str] Dummy Applied Curret Name.
            to_neuron_ID = 0;                                                                       % [#] Dummy Neuron ID to Which Current is Applied.
            ts = 0;                                                                                 % [s] Dummy Current Time Vector.
            Ias = 0;                                                                                % [A] Dummy Applied Current Magnitudes.
            enabled_flag = false;                                                                   % [F] Dummy Enabled Flags.
                        
            % Create a dummy applied current to calculate the applied current magnitude.
            [ ~, applied_current, ~, ~ ] = self.create_applied_current( ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, false, false, array_utilities );
            
            % Compute the magnitude for the output applied current.                        
            [ Ias_n, ~ ] = applied_current.compute_subtraction_Iasn( encoding_scheme, false, applied_current.applied_current_utilities );
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the inversion subnetwork output applied currents.
        function [ Ias2, applied_currents, self ] = compute_inversion_Ias2( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                	% [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                 % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.
            parameters = self.process_inversion_Ias2_parameters( parameters, encoding_scheme );
            
            % Retrieve the index associated with the output applied current.
            applied_current_index = self.get_applied_current_index( applied_current_IDs( end ), applied_currents, undetected_option );

            % Compute the magnitude for the output applied current.            
            [ Ias2, applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_inversion_Ias2( parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to compute the magnitude of the reduced inversion subnetwork output applied currents.
        function [ Ias2, applied_currents, self ] = compute_reduced_inversion_Ias2( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.
            parameters = self.process_reduced_inversion_Ias2_parameters( parameters, encoding_scheme );
            
            % Retrieve the index associated with the output applied current.
            applied_current_index = self.get_applied_current_index( applied_current_IDs( end ), applied_currents, undetected_option );

            % Compute the magnitude for the output applied current.
            [ Ias2, applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_reduced_inversion_Ias2( parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
                
        
        % ---------- Division Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the division subnetwork output applied currents.
        function Ias3 = compute_division_Ias3( self, encoding_scheme, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Define the dummy applied current properties.
            ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities );      % [#] Dummy Applied Current ID.
            name = 'Dummy Current';                                                                 % [str] Dummy Applied Curret Name.
            to_neuron_ID = 0;                                                                       % [#] Dummy Neuron ID to Which Current is Applied.
            ts = 0;                                                                                 % [s] Dummy Current Time Vector.
            Ias = 0;                                                                                % [A] Dummy Applied Current Magnitudes.
            enabled_flag = false;                                                                   % [F] Dummy Enabled Flags.
                        
            % Create a dummy applied current to calculate the applied current magnitude.
            [ ~, applied_current, ~, ~ ] = self.create_applied_current( ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, false, false, array_utilities );
            
            % Compute the magnitude for the output applied current.                        
            [ Ias3, ~ ] = applied_current.compute_division_Ias3( encoding_scheme, false, applied_current.applied_current_utilities );
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the reduced division subnetwork output applied currents.
        function Ias3 = compute_reduced_division_Ias3( self, encoding_scheme, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Define the dummy applied current properties.
            ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities );      % [#] Dummy Applied Current ID.
            name = 'Dummy Current';                                                                 % [str] Dummy Applied Curret Name.
            to_neuron_ID = 0;                                                                       % [#] Dummy Neuron ID to Which Current is Applied.
            ts = 0;                                                                                 % [s] Dummy Current Time Vector.
            Ias = 0;                                                                                % [A] Dummy Applied Current Magnitudes.
            enabled_flag = false;                                                                   % [F] Dummy Enabled Flags.
                        
            % Create a dummy applied current to calculate the applied current magnitude.
            [ ~, applied_current, ~, ~ ] = self.create_applied_current( ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, false, false, array_utilities );
            
            % Compute the magnitude for the output applied current.                        
            [ Ias3, ~ ] = applied_current.compute_reduced_division_Ias3( encoding_scheme, false, applied_current.applied_current_utilities );
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the division after inversion subnetwork output applied currents.
        function Ias3 = compute_dai_Ias3( self, encoding_scheme, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Define the dummy applied current properties.
            ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities );      % [#] Dummy Applied Current ID.
            name = 'Dummy Current';                                                                 % [str] Dummy Applied Curret Name.
            to_neuron_ID = 0;                                                                       % [#] Dummy Neuron ID to Which Current is Applied.
            ts = 0;                                                                                 % [s] Dummy Current Time Vector.
            Ias = 0;                                                                                % [A] Dummy Applied Current Magnitudes.
            enabled_flag = false;                                                                   % [F] Dummy Enabled Flags.
                        
            % Create a dummy applied current to calculate the applied current magnitude.
            [ ~, applied_current, ~, ~ ] = self.create_applied_current( ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, false, false, array_utilities );
            
            % Compute the magnitude for the output applied current.                        
            [ Ias3, ~ ] = applied_current.compute_dai_Ias3( encoding_scheme, false, applied_current.applied_current_utilities );
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of the reduced division after inversion subnetwork output applied currents.
        function Ias3 = compute_reduced_dai_Ias3( self, encoding_scheme, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Define the dummy applied current properties.
            ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities );      % [#] Dummy Applied Current ID.
            name = 'Dummy Current';                                                                 % [str] Dummy Applied Curret Name.
            to_neuron_ID = 0;                                                                       % [#] Dummy Neuron ID to Which Current is Applied.
            ts = 0;                                                                                 % [s] Dummy Current Time Vector.
            Ias = 0;                                                                                % [A] Dummy Applied Current Magnitudes.
            enabled_flag = false;                                                                   % [F] Dummy Enabled Flags.
                        
            % Create a dummy applied current to calculate the applied current magnitude.
            [ ~, applied_current, ~, ~ ] = self.create_applied_current( ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, false, false, array_utilities );
            
            % Compute the magnitude for the output applied current.                        
            [ Ias3, ~ ] = applied_current.compute_reduced_dai_Ias3( encoding_scheme, false, applied_current.applied_current_utilities );
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of multiplication subnetwork applied currents.
        function [ Ias3, applied_currents, self ] = compute_multiplication_Ias3( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.            
            parameters = self.process_multiplication_Ias3_parameters( parameters, encoding_scheme, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias3 = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the parameters associated with this applied current.
                these_parameters = { parameters{ 1 }{ k }, parameters{ 2 }{ k } };
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the magnitude for this applied current.
                [ Ias3( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_multiplication_Ias3( these_parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to compute the magnitude of reduced multiplication subnetwork applied currents.
        function [ Ias3, applied_currents, self ] = compute_reduced_multiplication_Ias3( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.            
            parameters = self.process_reduced_multiplication_Ias3_parameters( parameters, encoding_scheme, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias3 = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the parameters associated with this applied current.
                these_parameters = { parameters{ 1 }{ k }, parameters{ 2 }{ k } };
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the magnitude for this applied current.
                [ Ias3( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_reduced_multiplication_Ias3( these_parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % ---------- Centering Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of centering subnetwork applied currents.
        function [ Ias, applied_currents, self ] = compute_centering_Ias( self, applied_current_IDs, Gm, R, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                 	% [class] Array of Applied Current Class Objects.
            if nargin < 4, R = self.R_DEFAULT; end                                       	% [V] Activation Domain.
            if nargin < 3, Gm = self.Gm_DEFAULT; end                                        % [S] Membrane Conductance.
            if nargin < 2, applied_current_IDs = 'all'; end                                	% [-] Applied Current IDs.
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the magnitude for this applied current.
                [ Ias( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_centering_Ias( Gm, R, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of integration subnetwork applied currents.
        function [ Ias, applied_currents, self ] = compute_integration_Ias( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.
            parameters = self.process_integration_Ias_parameters( parameters, encoding_scheme, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the parameters associated with this applied current.
                these_parameters = { parameters{ 1 }{ k }, parameters{ 2 }{ k } };
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the magnitude for this applied current.
                [ Ias( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_integration_Ias( these_parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to compute the magnitude of voltage based integration subnetwork applied currents.
        function [ Ias, applied_currents, self ] = compute_vbi_Ias( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.
            parameters = self.process_vbi_Ias_parameters( parameters, encoding_scheme, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the parameters associated with this applied current.
                these_parameters = { parameters{ 1 }{ k }, parameters{ 2 }{ k } };
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the magnitude for this applied current.
                [ Ias( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_vbi_Ias( these_parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to compute the first magnitude of split voltage based integration subnetwork applied currents.
        function [ Ias, applied_currents, self ] = compute_svbi_Ias1( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.
            parameters = self.process_svbi_Ias1_parameters( parameters, encoding_scheme, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the parameters associated with this applied current.
                these_parameters = { parameters{ 1 }{ k }, parameters{ 2 }{ k } };
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the magnitude for this applied current.
                [ Ias( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_svbi_Ias1( these_parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to compute the second magnitude of split voltage based integration subnetwork applied currents.
        function [ Ias, applied_currents, self ] = compute_svbi_Ias2( self, applied_current_IDs, parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFUALT; end
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, parameters = {  }; end
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Process the parameters.
            parameters = self.process_svbi_Ias2_parameters( parameters, encoding_scheme, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the parameters associated with this applied current.
                these_parameters = { parameters{ 1 }{ k }, parameters{ 2 }{ k } };
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the magnitude for this applied current.
                [ Ias( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_svbi_Ias2( these_parameters, encoding_scheme, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to compute the time vector of multistate cpg applied currents.
        function [ ts, applied_currents, self ] = compute_mcpg_ts( self, applied_current_IDs, dt, tf, applied_currents, filter_disabled_flag, set_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 9, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 8, process_option = self.process_option_DEFAULT; end                                        % [str] Process Option. (Must be either 'max', 'min', 'mean', or 'none'.)
            if nargin < 7, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 6, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end                            % [T/F] Filter Disabled Flag. (Determines whether to considered disabled applied currents.)  
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, tf = self.tf_DEFAULT; end                                                                % [s] Simulation Duration
            if nargin < 3, dt = self.dt_DEFAULT; end                                                                % [s] Simulation Step Size
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Compute the number of timesteps associated with the given applied currents.
            num_timesteps = get_num_timesteps( applied_current_IDs( 1 ), applied_currents, filter_disabled_flag, process_option, undetected_option );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            ts = zeros( num_timesteps, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the time vector for this applied current.
                [ ts( :, k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_mcpg_ts( dt, tf, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to compute the magnitude of multistate cpg applied currents.
        function [ Ias, applied_currents, self ] = compute_mcpg_Ias( self, applied_current_IDs, dt, tf, applied_currents, filter_disabled_flag, set_flag, process_option, undetected_option )
            
            % Set the default input arguments.
            if nargin < 9, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 8, process_option = self.process_option_DEFAULT; end                                        % [str] Process Option. (Must be either 'max', 'min', 'mean', or 'none'.)
            if nargin < 7, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 6, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end                            % [T/F] Filter Disabled Flag. (Determines whether to considered disabled applied currents.)  
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, tf = self.tf_DEFAULT; end                                                                % [s] Simulation Duration
            if nargin < 3, dt = self.dt_DEFAULT; end                                                                % [s] Simulation Step Size
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Compute the number of timesteps associated with the given applied currents.
            num_timesteps = get_num_timesteps( applied_current_IDs( 1 ), applied_currents, filter_disabled_flag, process_option, undetected_option );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( num_timesteps, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the time vector for this applied current.
                [ Ias( :, k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_mcpg_Ias( dt, tf, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to compute the magnitude of driven multistate cpg applied currents.
        function [ Ias, applied_currents, self ] = compute_dmcpg_Ias( self, applied_current_IDs, Gm, R, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, R = self.R_DEFAULT; end                                                                  % [V] Activation Domain
            if nargin < 3, Gm = self.Gm_DEFAULT; end                                                                % [S] Membrane Conductance
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the time vector for this applied current.
                [ Ias( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_dmcpg_Ias( Gm, R, true, applied_currents( applied_current_index ).applied_current_utilities );
                
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        % Implement a function to compute the magnitude of the applied currents that connect the driven multistate cpg double centered lead lag and centered doube subtraction subnetworks.
        function [ Ias, applied_currents, self ] = compute_dmcpgdcll2cds_Ias( self, applied_current_IDs, Gm, R, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, R = self.R_DEFAULT; end                                                                  % [V] Activation Domain
            if nargin < 3, Gm = self.Gm_DEFAULT; end                                                                % [S] Membrane Conductance
            if nargin < 2, applied_current_IDs = 'all'; end                                                         % [-] Applied Current IDs
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Determine how many applied currents to which we are going to apply the given method.
            num_applied_currents_to_evaluate = length( applied_current_IDs );
            
            % Preallocate an array to store the time vectors associated with the applied currents.
            Ias = zeros( 1, num_applied_currents_to_evaluate );
            
            % Evaluate the given applied current method for each neuron.
            for k = 1:num_applied_currents_to_evaluate               % Iterate through each of the applied currents of interest...
                
                % Retrieve the index associated with this applied current ID.
                applied_current_index = self.get_applied_current_index( applied_current_IDs( k ), applied_currents, undetected_option );
                
                % Compute the time vector for this applied current.
                [ Ias( k ), applied_currents( applied_current_index ) ] = applied_currents( applied_current_index ).compute_dmcpgdcll2cds_Ias( Gm, R, true, applied_currents( applied_current_index ).applied_current_utilities );
               
            end
            
            % Determine whether to update the applied current manager.
            if set_flag, self.applied_currents = applied_currents; end
            
        end
        
        
        
        %% Validation Functions.
        
        % Implement a function to validate the compatibiltiy of applied current properties.
        function valid_flag = validate_applied_current_properties( self, n_applied_currents, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities )
            
            % Set the default input arguments.
            if nargin < 10, array_utilities = self.array_utilities; end                                                      % [class] Array Utilities Class.
            if nargin < 9, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 8, enabled_flags = true( 1, n_applied_currents ); end
            if nargin < 7, Ias = self.Ias_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_IDs_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 4, names = repmat( { '' }, 1, n_applied_currents ); end
            if nargin < 3, IDs = self.generate_unique_applied_current_IDs( n_applied_currents, applied_currents, array_utilities ); end
            
            % Determine whether to convert the names property to a cell.
            if ~iscell( names ), names = { names }; end
            
            % Determine whether the applied current properties are relevant.
            valid_flag = ( n_applied_currents == length( IDs ) ) && ( n_applied_currents == length( names ) &&  n_applied_currents == length( to_neuron_IDs ) ) && ( n_applied_currents == size( ts, 2 ) ) && ( n_applied_currents == size( Ias, 2 ) ) && ( n_applied_currents == length( enabled_flags ) );
            
        end
        
        
        %% Process Functions.
        
        % Implement a function to process applied current creation inputs.
        function [ n_applied_currents, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags ] = process_applied_current_creation_inputs( self, n_applied_currents, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 10, array_utilities = self.array_utilities; end                                                      % [class] Array Utilities Class.
            if nargin < 9, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 8, enabled_flags = true( 1, n_applied_currents ); end
            if nargin < 7, Ias = self.Ias_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_IDs_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 4, names = repmat( { '' }, 1, n_applied_currents ); end
            if nargin < 3, IDs = self.generate_unique_applied_current_IDs( n_applied_currents, applied_currents, array_utilities ); end
            
            % Convert the applied current parameters from cells to arrays as appropriate.
            enabled_flags = array_utilities.cell2array( enabled_flags );
            Ias = array_utilities.cell2array( Ias );
            ts = array_utilities.cell2array( ts );
            to_neuron_IDs = array_utilities.cell2array( to_neuron_IDs );
            names = array_utilities.cell2array( names );
            IDs = array_utilities.cell2array( IDs );
            n_applied_currents = array_utilities.cell2array( n_applied_currents );

            % Ensure that the applied current properties match the required number of applied currents.
            assert( self.validate_applied_current_properties( n_applied_currents, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities ), 'Provided applied current properties  must be of consistent size.' )
            
        end
        
        
        % Implement a function to process the applied current creation outputs.
        function [ IDs, applied_currents ] = process_applied_current_creation_outputs( self, IDs, applied_currents, as_cell_flag, array_utilities )
            
           	% Set the default input arguments.
            if nargin < 5, array_utilities = self.array_utilities; end                      % [class] Array Utilities Class.
            if nargin < 4, as_cell_flag = self.as_cell_flag_DEFAULT; end                  	% [T/F] As Cell Flag (Determines whether neurons are returned in an array or a cell.)
            if nargin < 3, applied_currents = self.applied_currents; end                  	% [class] Array of Applied Current Class Objects.
            
            % Determine whether to embed the new applied current IDs and objects in cells.
            if as_cell_flag                                                                 % If we want to embed the new applied current IDs and objects into cells...
                
                % Determine whether to embed the applied current IDs into a cell.
                if ~iscell( IDs )                                                           % If the IDs are not already a cell...
                
                    % Embed applied current IDs into a cell.
                    IDs = { IDs };
                
                end
                
                % Determine whether to embed the applied current objects into a cell.
                if ~iscell( applied_currents )                                                       % If the applied currents are not already a cell...
                
                    % Embed applied current objects into a cell.
                    applied_currents = { applied_currents };
                    
                end
                
            else                                                                            % Otherwise...
                
                % Determine whether to embed the applied current IDs into an array.
                if iscell( IDs )                                                            % If the applied current IDs are a cell...
                
                    % Convert the applied current IDs cell to a regular array.
                    IDs = array_utilities.cell2array( IDs );
                    
                end
                
                % Determine whether to embed the applied current objects into an array.
                if iscell( applied_currents )                                                        % If the applied current objects are a cell...
                
                    % Convert the applied current objects cell to a regular array.
                    applied_currents = array_utilities.cell2array( applied_currents );
                    
                end
                
            end
            
        end
        
        
        % Implement a function to update the applied current manager.
        function [ applied_currents, self ] = update_applied_current_manager( self, applied_currents, applied_current_manager, set_flag )
        
            % Set the default input arguments.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end            % [T/F] Set Flag (Determines whether output self object is updated.)
            
            % Determine whether to update the applied current manager object.
            if set_flag                                                  	% If we want to update the applied current manager object...
                
                % Update the applied current manager object.
                self = applied_current_manager;
            
            else                                                            % Otherwise...
                
                % Reset the applied currents object.
                applied_currents = self.applied_currents;
            
            end
            
        end
        
        
        %% Applied Current Parameter Packing Functions.
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to pack the parameters for an absolute inversion subnetwork.
        function inversion_parameters = pack_absolute_inversion_parameters( self, R2, Gm2 )
            
            % Set the default input arguments.
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, R2 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            inversion_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            inversion_parameters{ 1 } = R2;
            inversion_parameters{ 2 } = Gm2;            
            
        end
        
        
        % Implement a function to pack the parameters for a relative inversion subnetwork.
        function inversion_parameters = pack_relative_inversion_parameters( self, R2, Gm2 )
            
            % Set the default input arguments.
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, R2 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            inversion_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            inversion_parameters{ 1 } = R2;
            inversion_parameters{ 2 } = Gm2;            
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to pack the parameters for a reduced absolute inversion subnetwork.
        function reduced_inversion_parameters = pack_reduced_absolute_inversion_parameters( self, R2, Gm2 )
            
            % Set the default input arguments.
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, R2 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            reduced_inversion_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            reduced_inversion_parameters{ 1 } = R2;
            reduced_inversion_parameters{ 2 } = Gm2;            
            
        end
        
        
        % Implement a function to pack the parameters for a reduced relative inversion subnetwork.
        function reduced_inversion_parameters = pack_reduced_relative_inversion_parameters( self, R2, Gm2 )
            
            % Set the default input arguments.
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, R2 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            reduced_inversion_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            reduced_inversion_parameters{ 1 } = R2;
            reduced_inversion_parameters{ 2 } = Gm2;            
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to pack the parameters for an absolute multiplication subnetwork.
        function multiplication_parameters = pack_absolute_multiplication_parameters( self, R3, Gm3 )
            
            % Set the default input arguments.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            multiplication_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            multiplication_parameters{ 1 } = R3;
            multiplication_parameters{ 2 } = Gm3;            
            
        end
        
        
        % Implement a function to pack the parameters for a relative multiplication subnetwork.
        function multiplication_parameters = pack_relative_multiplication_parameters( self, R3, Gm3 )
            
            % Set the default input arguments.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            multiplication_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            multiplication_parameters{ 1 } = R3;
            multiplication_parameters{ 2 } = Gm3;            
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------
        
        % Implement a function to pack the parameters for a reduced absolute multiplication subnetwork.
        function reduced_multiplication_parameters = pack_reduced_absolute_multiplication_parameters( self, R3, Gm3 )
            
            % Set the default input arguments.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            reduced_multiplication_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            reduced_multiplication_parameters{ 1 } = R3;
            reduced_multiplication_parameters{ 2 } = Gm3;            
            
        end
        
        
        % Implement a function to pack the parameters for a reduced relative multiplication subnetwork.
        function reduced_multiplication_parameters = pack_reduced_relative_multiplication_parameters( self, R3, Gm3 )
            
            % Set the default input arguments.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Preallocate a cell array to store the parameters.
            reduced_multiplication_parameters = cell( 1, 2 );
            
            % Pack the parameters.
            reduced_multiplication_parameters{ 1 } = R3;
            reduced_multiplication_parameters{ 2 } = Gm3;            
            
        end
        
        
        %% Applied Current Unpacking Parameter Functions.
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to unpack absolute inversion subnetwork parameters.
        function [ R2, Gm2 ] = unpack_absolute_inversion_parameters( self, inversion_parameters )
            
            % Set the default input arguments.
            if nargin < 2, inversion_parameters = {  }; end                                                      	% [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( inversion_parameters )                                                                      % If the parameters are empty...
            
                % Set the parameters to default values.
                R2 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm2 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( inversion_parameters ) == 2                                                              % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R2 = inversion_parameters{ 1 };                                                                     % [V] Activation Domain.
                Gm2 = inversion_parameters{ 2 };                                                                    % [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack relative inversion subnetwork parameters.
        function [ R2, Gm2 ] = unpack_relative_inversion_parameters( self, inversion_parameters )
            
            % Set the default input arguments.
            if nargin < 2, inversion_parameters = {  }; end                                                      	% [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( inversion_parameters )                                                                      % If the parameters are empty...
            
                % Set the parameters to default values.
                R2 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm2 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( inversion_parameters ) == 2                                                              % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R2 = inversion_parameters{ 1 };                                                                     % [V] Activation Domain.
                Gm2 = inversion_parameters{ 2 };                                                                    % [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to unpack reduced absolute inversion subnetwork parameters.
        function [ R2, Gm2 ] = unpack_reduced_absolute_inversion_parameters( self, reduced_inversion_parameters )
            
            % Set the default input arguments.
            if nargin < 2, reduced_inversion_parameters = {  }; end                                                      	% [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( reduced_inversion_parameters )                                                                      % If the parameters are empty...
            
                % Set the parameters to default values.
                R2 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm2 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( reduced_inversion_parameters ) == 2                                                              % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R2 = reduced_inversion_parameters{ 1 };                                                                     % [V] Activation Domain.
                Gm2 = reduced_inversion_parameters{ 2 };                                                                    % [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack reduced relative inversion subnetwork parameters.
        function [ R2, Gm2 ] = unpack_reduced_relative_inversion_parameters( self, reduced_inversion_parameters )
            
            % Set the default input arguments.
            if nargin < 2, reduced_inversion_parameters = {  }; end                                              	% [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( reduced_inversion_parameters )                                                           	% If the parameters are empty...
            
                % Set the parameters to default values.
                R2 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm2 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( reduced_inversion_parameters ) == 2                                                     	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                R2 = reduced_inversion_parameters{ 1 };                                                            	% [V] Activation Domain.
                Gm2 = reduced_inversion_parameters{ 2 };                                                          	% [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to unpack absolute multiplication subnetwork parameters.
        function [ R3, Gm3 ] = unpack_absolute_multiplication_parameters( self, multiplication_parameters )
            
            % Set the default input arguments.
            if nargin < 2, multiplication_parameters = {  }; end                                                 	% [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( multiplication_parameters )                                                                	% If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( multiplication_parameters ) == 2                                                       	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = multiplication_parameters{ 1 };                                                              	% [V] Activation Domain.
                Gm3 = multiplication_parameters{ 2 };                                                             	% [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack relative multiplication subnetwork parameters.
         function [ R3, Gm3 ] = unpack_relative_multiplication_parameters( self, multiplication_parameters )
            
            % Set the default input arguments.
            if nargin < 2, multiplication_parameters = {  }; end                                                 	% [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( multiplication_parameters )                                                                	% If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( multiplication_parameters ) == 2                                                       	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = multiplication_parameters{ 1 };                                                              	% [V] Activation Domain.
                Gm3 = multiplication_parameters{ 2 };                                                             	% [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
           
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to unpack reduced absolute multiplication subnetwork parameters.
        function [ R3, Gm3 ] = unpack_reduced_absolute_multiplication_parameters( self, reduced_multiplication_parameters )
            
            % Set the default input arguments.
            if nargin < 2, reduced_multiplication_parameters = {  }; end                                            % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( reduced_multiplication_parameters )                                                       	% If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( reduced_multiplication_parameters ) == 2                                               	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = reduced_multiplication_parameters{ 1 };                                                       	% [V] Activation Domain.
                Gm3 = reduced_multiplication_parameters{ 2 };                                                   	% [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack reduced relative multiplication subnetwork parameters.
        function [ R3, Gm3 ] = unpack_reduced_relative_multiplication_parameters( self, reduced_multiplication_parameters )
            
            % Set the default input arguments.
            if nargin < 2, reduced_multiplication_parameters = {  }; end                                            % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( reduced_multiplication_parameters )                                                       	% If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                                                                % [V] Activation Domain.                                                                                            % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                                                              % [S] Membrane Conductance.                                                                                          % [S] Membrane Conductance.
                
            elseif length( reduced_multiplication_parameters ) == 2                                               	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = reduced_multiplication_parameters{ 1 };                                                       	% [V] Activation Domain.
                Gm3 = reduced_multiplication_parameters{ 2 };                                                   	% [S] Membrane Conductance.
            
            else                                                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
            
        
        %% Applied Current Creation Functions.
        
        % Implement a function to create a new applied current.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_applied_current( self, ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Set the default input arguments.
            if nargin < 11, array_utilities = self.array_utilities; end                                                     % [class] Array Utilities Class.
            if nargin < 10, as_cell_flag = self.as_cell_flag_DEFAULT; end                                                   % [T/F] As Cell Flag.
            if nargin < 9, set_flag = self.set_flag_DEFAULT; end                                                            % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 8, applied_currents = self.applied_currents; end                                                    % [class] Array of Applied Current Class Objects.
            if nargin < 7, enabled_flag = self.enabled_flag_DEFAULT; end                                                    % [T/F] Enabled Flag.
            if nargin < 6, Ias = self.Ias_DEFAULT; end                                                                      % [A] Applied Current Manitudes.
            if nargin < 5, ts = self.ts_DEFAULT; end                                                                        % [s] Applied Current Times.
            if nargin < 4, to_neuron_ID = self.to_neuron_ID_DEFAULT; end                                                    % [#] Neuron ID to which Current is Applied.
            if nargin < 3, name = ''; end                                                                                   % [str] Applied Current Name.
            if nargin < 2, ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end           % [#] Applied Current ID.
            
            % Process the applied current creation properties.
            [ ~, ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that this neuron ID is a unique natural.
            assert( self.unique_natural_applied_current_ID( ID, applied_currents, array_utilities ), 'Proposed applied current ID %0.2f is not a unique natural number.', ID )
            
            % Make an instance of the applied current manager.
            applied_current_manager = self;
            
            % Create an instance of the applied current class.
            applied_current_new = applied_current_class( ID, name, to_neuron_ID, ts, Ias, enabled_flag, array_utilities );
            
            % Retrieve the new applied current Id.
            ID_new = applied_current_new.ID;
            
            % Determine whether to embed the new applied current ID and object in cells.
            [ ID_new, applied_current_new ] = self.process_applied_current_creation_outputs( ID_new, applied_current_new, as_cell_flag, array_utilities );

            % Append this applied current to the array of existing applied currents.
            applied_currents = [ applied_currents, applied_current_new ];
            
            % Update the applied current manager to reflect the update applied currents object.
            applied_current_manager.applied_currents = applied_currents;
            applied_current_manager.num_applied_currents = length( applied_currents );
            
            % Determine whether to update the applied current manager object.
            [ applied_currents, self ] = self.update_applied_current_manager( applied_currents, applied_current_manager, set_flag );
            
        end
        

        % Implement a function to create multiple applied currents.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_applied_currents( self, n_applied_currents_to_create, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Set the default input arguments.
            if nargin < 7, enabled_flags = true( 1, n_applied_currents_to_create ); end
            if nargin < 6, Ias = self.Ias_DEFAULT*ones( 1, n_applied_currents_to_create ); end
            if nargin < 5, ts = self.ts_DEFAULT*ones( 1, n_applied_currents_to_create ); end
            if nargin < 4, to_neuron_IDs = self.to_neuron_ID_DEFAULT*ones( 1, n_applied_currents_to_create ); end
            if nargin < 3, names = repmat( { '' }, 1, n_applied_currents_to_create ); end
            if nargin < 2, IDs = self.generate_unique_applied_current_IDs( n_applied_currents_to_create, applied_currents, array_utilities ); end
            
            % Process the applied current creation inputs.
            [ n_applied_currents, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags ] = self.process_applied_current_creation_inputs( n_applied_currents_to_create, IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities );
            
            % Preallocate an array to store the new applied currents.
            applied_currents_new = repmat( applied_current_class(  ), [ 1, n_applied_currents_to_create ] );
            
            % Preallocate an array to store the new applied current IDs.
            IDs_new = zeros( 1, n_applied_currents );
            
            % Create an instance of the applied current manager that can be updated.
            applied_current_manager = self;
            
            % Create each of the specified applied currents.
            for k = 1:n_applied_currents                                                                                           % Iterate through each of the applied currents we want to create...
                
                % Create this applied current.
                [ IDs_new( k ), applied_currents_new( k ), applied_currents, applied_current_manager ] = applied_current_manager.create_applied_current( IDs( k ), names { k }, to_neuron_IDs( k ), ts( :, k ), Ias( :, k ), enabled_flags( k ), applied_currents, true, false, array_utilities );
                
            end
            
            % Determine whether to embed the new applied current ID and object in cells.
            [ IDs_new, applied_currents_new ] = self.process_applied_current_creation_outputs( IDs_new, applied_currents_new, as_cell_flag, array_utilities );
            
            % Determine whether to update the applied current manager object.
            [ applied_currents, self ] = self.update_applied_current_manager( applied_currents, applied_current_manager, set_flag );
            
        end
        
        
        % Implement a function to delete an applied current.
        function [ applied_currents, self ] = delete_applied_current( self, applied_current_ID, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Create an instance of the applied current manager.
            applied_current_manager = self;
            
            % Retrieve the index associated with this applied current.
            applied_current_index = self.get_applied_current_index( applied_current_ID, applied_currents, undetected_option );
            
            % Remove this applied current from the array of applied currents.
            applied_currents( applied_current_index ) = [  ];
            
            % Update the applied current manager to reflect these changes.
            applied_current_manager.applied_currents = applied_currents;
            applied_current_manager.num_applied_currents = length( applied_currents );
            
           % Determine whether to update the applied currents and applied current manager objects.
           [ applied_currents, self ] = self.update_applied_current_manager( applied_currents, applied_current_manager, set_flag );

        end
        
        
        % Implement a function to delete multiple applied currents.
        function [ applied_currents, self ] = delete_applied_currents( self, applied_current_IDs, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 5, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, applied_current_IDs = 'all'; end
            
            % Validate the applied current IDs.
            applied_current_IDs = self.validate_applied_current_IDs( applied_current_IDs, applied_currents );
            
            % Retrieve the number of applied currents to delete.
            num_applied_currents_to_delete = length( applied_current_IDs );
            
            % Delete each of the specified applied currents.
            for k = 1:num_applied_currents_to_delete                      % Iterate through each of the applied currents we want to delete...
                
                % Delete this applied current.
                [ applied_currents, self ] = self.delete_applied_current( applied_current_IDs( k ), applied_currents, set_flag, undetected_option );
                
            end
            
        end
        
        
        %% Subnetwork Applied Current Creation Functions
        
        %{
        % Implement a function to create the applied currents for an transmission subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_transmission_applied_currents( self, applied_currents, as_cell_flag )
        
            % Set the default input arguments.
            if nargin < 3, as_cell_flag = self.as_cell_flag_DEFAULT; end
            if nargin < 2, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Determine how to generate the applied current IDs and objects.
            if as_cell_flag                     % If we want the applied current IDs and objects to be cells...
                
                % Set the new applied current IDs and objects to be empty cells.
                IDs_new = {  };
                applied_currents_new = {  };
            
            else                                % Otherwise...
                
                % Set the new applied current IDs and obejcts to be empty arrays.
                IDs_new = [  ];
                applied_currents_new = [  ];
                
            end
        
        end
        
        
        % Implement a function to create the applied currents for an addition subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_addition_applied_currents( self, applied_currents, as_cell_flag )
        
            % Set the default input arguments.
            if nargin < 3, as_cell_flag = self.as_cell_flag_DEFAULT; end
            if nargin < 2, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Determine how to generate the applied current IDs and objects.
            if as_cell_flag                     % If we want the applied current IDs and objects to be cells...
                
                % Set the new applied current IDs and objects to be empty cells.
                IDs_new = {  };
                applied_currents_new = {  };
            
            else                                % Otherwise...
                
                % Set the new applied current IDs and obejcts to be empty arrays.
                IDs_new = [  ];
                applied_currents_new = [  ];
                
            end
        
        end
        
        
        % Implement a function to create the applied currents for a subtraction subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_subtraction_applied_currents( self, applied_currents, as_cell_flag )
        
            % Set the default input arguments.
            if nargin < 3, as_cell_flag = self.as_cell_flag_DEFAULT; end
            if nargin < 2, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Determine how to generate the applied current IDs and objects.
            if as_cell_flag                     % If we want the applied current IDs and objects to be cells...
                
                % Set the new applied current IDs and objects to be empty cells.
                IDs_new = {  };
                applied_currents_new = {  };
            
            else                                % Otherwise...
                
                % Set the new applied current IDs and obejcts to be empty arrays.
                IDs_new = [  ];
                applied_currents_new = [  ];
                
            end
        
        end
        
        %}
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to create the applied currents for an inversion subnetwork.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_inversion_applied_current( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
        
            % Set the number of neurons.
            n_neurons = self.num_inversion_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flags_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = ''; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_ID == self.to_neuron_ID_DEFAULT, to_neuron_ID = neuron_IDs( 2 ); end

            % Determine whether to compute the name.
            if isempty( name ), name = sprintf( 'Inversion Applied Current %0.0f', to_neuron_ID ); end
            
            % Create the subnetwork applied current.
            [ ID_new, applied_current_new, applied_currents, self ] = self.create_applied_current( applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
         
        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to create the applied currents for a reduced inversion subnetwork.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_reduced_inversion_applied_current( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
        
            % Set the number of neurons.
            n_neurons = self.num_inversion_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flags_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = ''; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_ID == self.to_neuron_ID_DEFAULT, to_neuron_ID = neuron_IDs( 2 ); end

            % Determine whether to compute the name.
            if isempty( name ), name = sprintf( 'Reduced Inversion Applied Current %0.0f', to_neuron_ID ); end
            
            % Create the subnetwork applied current.
            [ ID_new, applied_current_new, applied_currents, self ] = self.create_applied_current( applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
       %{ 
        % Implement a function to create the applied currents for a division subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_division_applied_currents( self, applied_currents, as_cell_flag )
            
            % Set the default input arguments.
            if nargin < 3, as_cell_flag = self.as_cell_flag_DEFAULT; end
            if nargin < 2, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            
            % Determine how to generate the applied current IDs and objects.
            if as_cell_flag                     % If we want the applied current IDs and objects to be cells...
                
                % Set the new applied current IDs and objects to be empty cells.
                IDs_new = {  };
                applied_currents_new = {  };
            
            else                                % Otherwise...
                
                % Set the new applied current IDs and obejcts to be empty arrays.
                IDs_new = [  ];
                applied_currents_new = [  ];
                
            end
            
        end
        %}
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to create the applied currents for a multiplication subnetwork.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_multiplication_applied_currents( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Set the number of neurons.
            n_neurons = self.num_multiplication_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flags_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = ''; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_ID == self.to_neuron_ID_DEFAULT, to_neuron_ID = neuron_IDs( 3 ); end

            % Determine whether to compute the name.
            if isempty( name ), name = sprintf( 'Multiplication Applied Current %0.0f', to_neuron_ID ); end
            
            % Create the subnetwork applied current.
            [ ID_new, applied_current_new, applied_currents, self ] = self.create_applied_current( applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
                
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to create the applied currents for a reduced multiplication subnetwork.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_reduced_multiplication_applied_currents( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Set the number of neurons.
            n_neurons = self.num_multiplication_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flags_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = ''; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_ID == self.to_neuron_ID_DEFAULT, to_neuron_ID = neuron_IDs( 3 ); end

            % Determine whether to compute the name.
            if isempty( name ), name = sprintf( 'Reduced Multiplication Applied Current %0.0f', to_neuron_ID ); end
            
            % Create the subnetwork applied current.
            [ ID_new, applied_current_new, applied_currents, self ] = self.create_applied_current( applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        % Implement a function to create the applied currents for an integration subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_integration_applied_currents( self, neuron_IDs, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Set the number of neurons.
            n_neurons = self.num_integration_neurons_DEFAULT;
            n_applied_currents = self.n_integration_applied_currents_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flags = self.enabled_flags_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 7, Ias = self.Ias_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, names = ''; end
            if nargin < 3, applied_current_IDs = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags ] = self.process_applied_current_creation_inputs( n_applied_currents, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_IDs == self.to_neuron_ID_DEFAULT, to_neuron_IDs = [ neuron_IDs( 1 ), neuron_IDs( 2 ) ]; end

            % Determine whether to compute the name.
            if all( [ names{ : } ] == '' )
                
                names = cells( 1, n_applied_currents );
                
                for k = 1:n_applied_currents
                
                    names{ k } = sprintf( 'Integration Applied Current %0.0f', to_neuron_IDs( k ) );
            
                end
            
            end
            
            % Create the subnetwork applied current.
            [ IDs_new, applied_currents_new, applied_currents, self ] = self.create_applied_current( applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
        % Implement a function to create the applied currents for a voltage based integration subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_vbi_applied_currents( self, neuron_IDs, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Set the number of neurons.
            n_neurons = self.num_integration_neurons_DEFAULT;
            n_applied_currents = self.n_integration_applied_currents_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flags = self.enabled_flags_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 7, Ias = self.Ias_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, names = ''; end
            if nargin < 3, applied_current_IDs = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags ] = self.process_applied_current_creation_inputs( n_applied_currents, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_IDs == self.to_neuron_ID_DEFAULT, to_neuron_IDs = [ neuron_IDs( 3 ), neuron_IDs( 4 ) ]; end

            % Determine whether to compute the name.     
            if isempty( names ), names = { sprintf( 'VBI Applied Current %0.0f', to_neuron_ID( 3 ) ), sprintf( 'VBI Applied Current %0.0f', to_neuron_ID( 4 ) ) }; end
             
            % Create the subnetwork applied current.
            [ IDs_new, applied_currents_new, applied_currents, self ] = self.create_applied_current( applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
        % Implement a function to create the applied currents for a split voltage based integration subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_svbi_applied_currents( self, neuron_IDs, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Set the number of neurons.
            n_neurons = self.num_integration_neurons_DEFAULT;
            n_applied_currents = self.n_integration_applied_currents_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flags = self.enabled_flags_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 7, Ias = self.Ias_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, names = ''; end
            if nargin < 3, applied_current_IDs = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags ] = self.process_applied_current_creation_inputs( n_applied_currents, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_IDs == self.to_neuron_ID_DEFAULT, to_neuron_IDs = [ neuron_IDs( 3 ), neuron_IDs( 4 ), neuron_IDs( 9 ) ]; end

            % Determine whether to compute the name.     
            if isempty( names ), names = { sprintf( 'SVBI Applied Current %0.0f', to_neuron_ID( 3 ) ), sprintf( 'SVBI Applied Current %0.0f', to_neuron_ID( 4 ) ), sprintf( 'SVBI Applied Current %0.0f', to_neuron_ID( 9 ) ) }; end
             
            % Create the subnetwork applied current.
            [ IDs_new, applied_currents_new, applied_currents, self ] = self.create_applied_current( applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
        % Implement a function to create the applied currents for a modulated split voltage based integration subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_msvbi_applied_currents( self, neuron_IDs, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Create the applied currents for a modulated split voltage based integration subnetwork.
            [ IDs_new, applied_currents_new, applied_currents, self ] = self.create_svbi_applied_currents( neuron_IDs, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
        % Implement a function to create the applied currents for a modulated split voltage based integration subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_mssvbi_applied_currents( self, neuron_IDs, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Create the modulated split voltage based integration applied currents.
            [ IDs_new, applied_currents_new, applied_currents, self ] = self.create_msvbi_applied_currents( neuron_IDs( 5:end ), applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities );

        end
        
        
        % ---------- Centering Subnetwork Functions ----------
        
        % Implement a function to create the applied currents for a centering subnetwork.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_centering_applied_currents( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
       
            % Set the number of neurons.
            n_neurons = self.num_centering_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flags_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = ''; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_ID == self.to_neuron_ID_DEFAULT, to_neuron_ID = neuron_IDs( 2 ); end

            % Determine whether to compute the name.
            if isempty( name ), name = sprintf( 'Centering Applied Current %0.0f', to_neuron_ID ); end
            
            % Create the subnetwork applied current.
            [ ID_new, applied_current_new, applied_currents, self ] = self.create_applied_current( applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
            
        
        % Implement a function to create the applied currents for a double centering subnetwork.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_double_centering_applied_currents( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
       
            % Set the number of neurons.
            n_neurons = self.num_double_centering_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flags_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = ''; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Determine whether to compute the to neuron ID.
            if to_neuron_ID == self.to_neuron_ID_DEFAULT, to_neuron_ID = neuron_IDs( 2 ); end

            % Determine whether to compute the name.
            if isempty( name ), name = sprintf( 'Double Centering Applied Current %0.0f', to_neuron_ID ); end
            
            % Create the subnetwork applied current.
            [ ID_new, applied_current_new, applied_currents, self ] = self.create_applied_current( applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
        % Implement a function to create the applied currents for a centered double subtraction subnetwork.
        function [ ID_new, applied_current_new, applied_currents, self ] = create_cds_applied_currents( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
        
            % Define the number of neurons from the various subnetworks.
            n_ds_neurons = self.num_double_subtraction_neurons_DEFAULT;                                                     % [#] Number of DS Neurons.
            n_dc_neurons = self.num_double_centering_neurons_DEFAULT;                                                       % [#] Number of DC Neurons.
            n_neurons = n_ds_neurons + n_dc_neurons;                                                                        % [#] Number of Neurons.
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flags_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = ''; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the applied current creation inputs.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Ensure that the neuron properties match the require number of neurons.
            assert( n_neurons == length( neuron_IDs ), 'Provided neuron properties must be of consistent size.' )
            
            % Define the neuron indexes.
            i_start_dc_neurons = i_end_ds_neurons + 1;
            i_end_dc_neurons = i_end_ds_neurons + n_dc_neurons;

            % Create the double centering applied currents.
            [ ID_new, applied_current_new, applied_currents, self ] = self.create_double_centering_applied_currents( neuron_IDs( i_start_dc_neurons:i_end_dc_neurons ), applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        

        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to create the applied currents for a multistate CPG subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_mcpg_applied_currents( self, neuron_IDs, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Compute the number of neurons and the number of applied currents.
            if nargin < 2, n_neurons = self.num_mcpg_neurons_DEFAULT; else, n_neurons = length( neuron_IDs ); end
            
            % Set the default input arguments.
            if nargin < 8, enabled_flag = self.enabled_flag_DEFAULT; end
            if nargin < 7, Ias = self.Ias_DEFAULT; end
            if nargin < 6, ts = self.ts_DEFAULT; end
            if nargin < 5, to_neuron_ID = self.to_neuron_ID_DEFAULT; end
            if nargin < 4, name = { '' }; end
            if nargin < 3, applied_current_ID = self.generate_unique_applied_current_ID( applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the input information.
            [ ~, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag ] = self.process_applied_current_creation_inputs( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, array_utilities );
            
            % Determine whether to use default to neuron IDs.
            if to_neuron_ID == self.to_neuron_ID_DEFAULT, to_neuron_ID = neuron_IDs( 1 ); end
            
            % Determine whether to use default names.
            if isempty( name ), name = sprintf( 'MCPG Applied Current %0.0f -> Neuron %0.0f', applied_current_ID, to_neuron_ID ); end

            % Create the subnetwork applied currents.
            [ IDs_new, applied_currents_new, applied_currents, self ] = self.create_applied_currents( 1, applied_current_ID, name, to_neuron_ID, ts, Ias, enabled_flag, applied_currents, set_flag, as_cell_flag, array_utilities );
            
        end
        
        
        % Implement a function to create the applied currents for a driven multistate CPG subnetwork.
        function [ IDs_new, applied_currents_new, applied_currents, self ] = create_dmcpg_applied_currents( self, neuron_IDs, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, set_flag, as_cell_flag, array_utilities )
            
            % Compute the number of neurons and the number of applied currents.
            if nargin < 2, n_neurons = self.num_dmcpg_neurons_DEFAULT; else, n_neurons = length( neuron_IDs ); end
            n_applied_currents = self.num_dmcpg_applied_currents_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 8, enabled_flags = self.enabled_flag_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 7, Ias = self.Ias_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 6, ts = self.ts_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 5, to_neuron_IDs = self.to_neuron_ID_DEFAULT*ones( 1, n_applied_currents ); end
            if nargin < 4, names = repmat( { '' }, 1, n_applied_currents ); end
            if nargin < 3, applied_current_IDs = self.generate_unique_applied_current_IDs( n_applied_currents, applied_currents, array_utilities ); end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the input information.
            [ ~, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags ] = self.process_applied_current_creation_inputs( n_applied_currents, applied_current_IDs, names, to_neuron_IDs, ts, Ias, enabled_flags, applied_currents, array_utilities );
            
            % Preallocate a cell array to store the new applied current IDs and objects.
            IDs_new = cell( 1, 2 );
            applied_currents_new = cell( 1, 2 );
            
            % Define the indexes for the mcpg subnetwork.
            i_start_mcpg_neurons = 1; i_end_mcpg_neurons = n_neurons - 1;
            i_start_mcpg_applied_currents = 1; i_end_mcpg_applied_currents = 1;
            
            % Create the applied currents for a multistate cpg subnetwork.
            [ IDs_new{ 1 }, applied_currents_new{ 1 }, applied_currents, applied_current_manager ] = self.create_mcpg_applied_currents( neuron_IDs( i_start_mcpg_neurons:i_end_mcpg_neurons ), applied_current_IDs( i_start_mcpg_applied_currents:i_end_mcpg_applied_currents ), names{ i_start_mcpg_applied_currents:i_end_mcpg_applied_currents }, to_neuron_IDs( i_start_mcpg_applied_currents:i_end_mcpg_applied_currents ), ts( :, i_start_mcpg_applied_currents:i_end_mcpg_applied_currents ), Ias( :, i_start_mcpg_applied_currents:i_end_mcpg_applied_currents ), enabled_flags( i_start_mcpg_applied_currents:i_end_mcpg_applied_currents ), applied_currents, true, false, array_utilities );
            
            % Determine whether to update the default or provided to neuron ID.
            if all( to_neuron_IDs( end ) == -1 )
            
                % Set the to neuron ID.
                to_neuron_IDs( end ) = neuron_IDs( end );
                
            end
            
            % Determine whether to update the default or provided CPG drive neuron name.
            if isempty( names{ end } )                                                                    % If the final name is empty...
                
                % Set the drive neuron name.
                names{ end } = sprintf( 'DMCPG Applied Current %0.0f -> Neuron %0.0f', applied_current_IDs( end ), to_neuron_IDs( end ) );

            end 
            
            % Create the an additional applied current for the drive neuron.
            [ IDs_new{ 2 }, applied_currents_new{ 2 }, applied_currents, applied_current_manager ] = applied_current_manager.create_applied_current( applied_current_IDs( end ), names{ end }, to_neuron_IDs( end ), ts( :, end ), Ias( :, end ), enabled_flags( end ), applied_currents, true, false, array_utilities );
            
            % Determine how to format the neuron IDs and objects.
            [ IDs_new, applied_currents_new ] = self.process_applied_current_creation_outputs( IDs_new, applied_currents_new, as_cell_flag, array_utilities );
            
            % Update the applied current manager and applied current objects as appropriate.
            [ applied_currents, self ] = self.update_applied_current_manager( applied_currents, applied_current_manager, set_flag );
            
        end
        
        
        %{
%         % Implement a function to create the applied currents for a driven multistate CPG split lead lag subnetwork.
%         function [ self, applied_current_IDs_cell ] = create_dmcpg_sll_applied_currents( self, to_neuron_IDs_cell )
%             
%             % Retrieve the number of subnetworks and the number of cpg neurons.
%             num_subnetworks = length( to_neuron_IDs_cell );
%             num_cpg_neurons = length( to_neuron_IDs_cell{ 1 } ) - 1;
%             
%             % Preallocate a cell array to store the applied current IDs.
%             applied_current_IDs_cell = cell( 1, num_subnetworks - 1 );
%             
%             % Create the applied currents for the driven multistate cpg subnetworks.
%             [ self, applied_current_IDs_cell{ 1 } ] = self.create_dmcpg_applied_currents( to_neuron_IDs_cell{ 1 } );
%             [ self, applied_current_IDs_cell{ 2 } ] = self.create_dmcpg_applied_currents( to_neuron_IDs_cell{ 2 } );
% 
%             % Create the applied currents for each of the modulated split subtraction voltage based integration subnetworks.
%             for k = 1:num_cpg_neurons                       % Iterate through each of the CPG neurons...
%                 
%                 % Create the applied currents for each of the modulated split subtraction voltage based integration subnetworks.
%                 [ self, applied_current_IDs_cell{ k + 2 } ] = self.create_mssvbi_applied_currents( to_neuron_IDs_cell{ k + 2 } );
%                 
%             end
%             
%         end
%         
%         
%         % Implement a function to create the applied currents for a driven multistate CPG double centered lead lag subnetwork.
%         function [ self, applied_current_IDs_cell ] = create_dmcpg_dcll_applied_currents( self, to_neuron_IDs_cell )
%             
%             % Create the applied currents for a driven multistate cpg split lead lag subnetwork.
%             [ self, applied_current_IDs_dmcpgsll ] = self.create_dmcpg_sll_applied_currents( to_neuron_IDs_cell{ 1 } );
%             
%             % Create the applied currents for the double centering subnetwork.
%             [ self, applied_current_IDs_dc ] = self.create_double_centering_applied_currents( to_neuron_IDs_cell{ 2 } );
%             
%             % Concatenate the applied current IDs.
%             applied_current_IDs_cell = { applied_current_IDs_dmcpgsll, applied_current_IDs_dc };
%             
%         end
%         
%         
%         % Implement a function to create the applied currents that connect the driven multistate cpg double centered lead lag subnetwork to the centered double subtraction subnetwork.
%         function [ self, applied_current_ID ] = create_dmcpgdcll2cds_applied_current( self, to_neuron_IDs_cell )
%             
%             % Create the applied current.
%             [ self, applied_current_ID ] = self.create_applied_current(  );
%             
%             % Set the name of the applied current.
%             [ applied_currents, applied_current_manager ] = self.set_applied_current_property( applied_current_ID, { 'Desired Lead / Lag' }, 'name', applied_currents, true );
%             
%             % Attach this applied current to a neuron.
%             [ applied_currents, applied_current_manager ] = applied_current_manager.set_applied_current_property( applied_current_ID, to_neuron_IDs_cell{ 3 }, 'to_neuron_ID', applied_currents, true );
%             
%             % Determine whether to update the applied current manager.
%             if set_flag, self = applied_current_manager; end
%             
%         end
%         
%         
%         % Implement a function to create the applied currents for an open loop driven multistate CPG double centered lead lag error subnetwork.
%         function [ self, applied_current_IDs_cell ] = create_ol_dmcpg_dclle_applied_currents( self, to_neuron_IDs_cell )
%             
%             % Create the applied currents for the driven multistate cpg double centered lead lag subnetwork.
%             [ self, applied_current_IDs_dmcpgdcll ] = self.create_dmcpg_dcll_applied_currents( to_neuron_IDs_cell{ 1 } );
%             
%             % Create the applied currents for the centered double subtraction subnetwork.
%             [ self, applied_current_IDs_cds ] = self.create_cds_applied_currents( to_neuron_IDs_cell{ 2 } );
%             
%             % Create the applied currents associated with connecting the driven multistate cpg double centered lead lag subnetwork to the centered double subtraction subnetwork.
%             [ self, applied_current_ID ] = self.create_dmcpgdcll2cds_applied_current( to_neuron_IDs_cell );
%             
%             % Concatenate the applied current IDs.
%             applied_current_IDs_cell = { applied_current_IDs_dmcpgdcll, applied_current_IDs_cds, applied_current_ID };
%             
%         end
%         
%         
%         % Implement a function to create the applied currents for a closed loop P controlled driven multistate CPG double centered lead lag subnetwork.
%         function [ self, applied_current_IDs_cell ] = create_clpc_dmcpg_dcll_applied_currents( self, to_neuron_IDs_cell )
            
            % Create the applied currents for an open loop driven multistate CPG double centered lead lag error subnetwork.
            [ self, applied_current_IDs_cell ] = self.create_ol_dmcpg_dclle_applied_currents( to_neuron_IDs_cell );
            
        end
        %}
        
        
        %% Subnetwork Applied Current Design Functions
                
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to design the applied currents for a transmission subnetwork.
        function Ias2 = design_transmission_applied_current( self, encoding_scheme, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class Object.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Compute the applied current magnitudes of this subnetwork.            
            Ias2 = self.compute_transmission_Ias2( encoding_scheme, applied_currents, array_utilities );
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
            
        % Implement a function to design the applied currents for a addition subnetwork.
        function Ias_n = design_addition_applied_current( self, encoding_scheme, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class Object.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Compute the applied current magnitudes of this subnetwork.            
            Ias_n = self.compute_addition_Iasn( encoding_scheme, applied_currents, array_utilities );
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------

        % Implement a function to design the applied currents for a subtraction subnetwork.
        function Ias_n = design_subtraction_applied_current( self, encoding_scheme, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class Object.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Compute the applied current magnitudes of this subnetwork.
            Ias_n = self.compute_subtraction_Iasn( encoding_scheme, applied_currents, array_utilities );
            
        end
           
        
        % ---------- Inversion Subnetwork Functions ----------

        % Implement a function to design the applied currents for an inversion subnetwork.
        function [ Ias2, applied_currents, self ] = design_inversion_applied_current( self, neuron_IDs, inversion_parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Compute the number of neurons.
            n_neurons = self.num_inversion_neurons;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, inversion_parameters = {  }; end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Retrieve the applied current IDs associated with the provided neuron IDs.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( neuron_IDs, applied_currents, undetected_option );
            
            % Process the parameters.
            inversion_parameters = self.process_inversion_Ias2_parameters( inversion_parameters, encoding_scheme );
            
            % Compute the inversion applied current magnitude outputs.
            [ Ias2, applied_currents, self ] = self.compute_inversion_Ias2( applied_current_IDs, inversion_parameters, encoding_scheme, applied_currents, set_flag, undetected_option );
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to design the applied currents for a reduced inversion subnetwork.
        function [ Ias2, applied_currents, self ] = design_reduced_inversion_applied_current( self, neuron_IDs, reduced_inversion_parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Compute the number of neurons.
            n_neurons = self.num_reduced_inversion_neurons;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, reduced_inversion_parameters = {  }; end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Retrieve the applied current IDs associated with the provided neuron IDs.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( neuron_IDs, applied_currents, undetected_option );
            
            % Process the parameters.
            reduced_inversion_parameters = self.process_reduced_inversion_Ias2_parameters( reduced_inversion_parameters, encoding_scheme );
            
            % Compute the inversion applied current magnitude outputs.
            [ Ias2, applied_currents, self ] = self.compute_reduced_inversion_Ias2( applied_current_IDs, reduced_inversion_parameters, encoding_scheme, applied_currents, set_flag, undetected_option );
            
        end
        
        
        % ---------- Division Subnetwork Functions ----------

        % Implement a function to design the applied currents for a division subnetwork.
        function Ias3 = design_division_applied_current( self, encoding_scheme, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class Object.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Compute the applied current magnitudes of this subnetwork.            
            Ias3 = self.compute_division_Ias3( encoding_scheme, applied_currents, array_utilities );
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions ----------

        % Implement a function to design the reduced applied currents for a reduced division subnetwork.
        function Ias3 = design_reduced_division_applied_current( self, encoding_scheme, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class Object.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Compute the applied current magnitudes of this subnetwork.            
            Ias3 = self.compute_reduced_division_Ias3( encoding_scheme, applied_currents, array_utilities );
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions ----------

        % Implement a function to design the applied currents for a division after inversion subnetwork.
        function Ias3 = design_dai_applied_current( self, encoding_scheme, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class Object.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Compute the applied current magnitudes of this subnetwork.            
            Ias3 = self.compute_dai_Ias3( encoding_scheme, applied_currents, array_utilities );
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------

        % Implement a function to design the applied currents for a reduced division after inversion subnetwork.
        function Ias3 = design_reduced_dai_applied_current( self, encoding_scheme, applied_currents, array_utilities )
        
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class Object.
            if nargin < 3, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Compute the applied current magnitudes of this subnetwork.            
            Ias3 = self.compute_dai_Ias3( encoding_scheme, applied_currents, array_utilities );
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to design the applied currents for a multiplication subnetwork.
        function [ Ias3, applied_currents, self ] = design_multiplication_applied_current( self, neuron_IDs, multiplication_parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Compute the number of multiplication neurons.
            n_neurons = self.num_multiplication_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, multiplication_parameters = {  }; end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Retrieve the applied current IDs associated with the provided neuron IDs.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( neuron_IDs, applied_currents, undetected_option );
            
            % Process the parameters.
            multiplication_parameters = self.process_multiplication_Ias3_parameters( multiplication_parameters, encoding_scheme, applied_currents );
            
            % Compute the multiplication applied current magnitude outputs.
            [ Ias3, applied_currents, self ] = self.compute_multiplication_Ias3( applied_current_IDs, multiplication_parameters, encoding_scheme, applied_currents, set_flag, undetected_option );
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to design the applied currents for a reduced multiplication subnetwork.
        function [ Ias3, applied_currents, self ] = design_reduced_multiplication_applied_current( self, neuron_IDs, reduced_multiplication_parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Compute the number of reduced multiplication neurons.
            n_neurons = self.num_reduced_multiplication_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end                  % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, reduced_multiplication_parameters = {  }; end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Retrieve the applied current IDs associated with the provided neuron IDs.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( neuron_IDs, applied_currents, undetected_option );
            
            % Process the parameters.
            reduced_multiplication_parameters = self.process_reduced_multiplication_Ias3_parameters( reduced_multiplication_parameters, encoding_scheme, applied_currents );
            
            % Compute the multiplication applied current magnitude outputs.
            [ Ias3, applied_currents, self ] = self.compute_reduced_multiplication_Ias3( applied_current_IDs, reduced_multiplication_parameters, encoding_scheme, applied_currents, set_flag, undetected_option );
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        % Implement a function to design the applied currents for an integration subnetwork.
        function [ Ias, applied_currents, self ] = design_integration_applied_currents( self, neuron_IDs, integration_parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Compute the number of neurons.
            n_neurons = self.num_integration_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, integration_parameters = {  }; end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the parameters.
            integration_parameters = self.process_integration_Ias_parameters( integration_parameters, encoding_scheme, applied_currents );
            
            % Get the applied current IDs that comprise this integration subnetwork.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( neuron_IDs, applied_currents, undetected_option );
            
            % Compute the applied current magnitudes associated with this subnetwork.
            [ Ias, applied_currents, self ] = self.compute_integration_Ias( applied_current_IDs, integration_parameters, encoding_scheme, applied_currents, set_flag, undetected_option );
            
        end
        
        
        % Implement a function to design the applied currents for a voltage based integration subnetwork.
        function [ Ias, applied_currents, self ] = design_vbi_applied_currents( self, neuron_IDs, vbi_parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
            % Compute the number of neurons.
            n_neurons = self.num_vbi_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, vbi_parameters = {  }; end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the parameters.
            vbi_parameters = self.process_vbi_Ias_parameters( vbi_parameters, encoding_scheme, applied_currents );
            
            % Get the applied current IDs that comprise this voltage based integration subnetwork.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( neuron_IDs, applied_currents, undetected_option );
            
            % Compute the applied current magnitudes associated with this subnetwork.
            [ Ias, applied_currents, self ] = self.compute_vbi_Ias( applied_current_IDs, vbi_parameters, encoding_scheme, applied_currents, set_flag, undetected_option );
           
        end
        
        
        % Implement a function to design the applied currents for a split voltage based integration subnetwork.
        function [ Ias, applied_currents, self ] = design_svbi_applied_currents( self, neuron_IDs, svbi_parameters, encoding_scheme, applied_currents, set_flag, undetected_option )
            
             % Compute the number of neurons.
            n_neurons = self.num_svbi_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, encoding_scheme = self.encoding_scheme_DEFAULT; end
            if nargin < 3, svbi_parameters = {  }; end
            if nargin < 2, neuron_IDs = 1:n_neurons; end
            
            % Process the parameters.
            parameters_Ias1 = self.process_svbi_Ias1_parameters( svbi_parameters{ 1 }, encoding_scheme, applied_currents );
            parameters_Ias2 = self.process_svbi_Ias2_parameters( svbi_parameters{ 2 }, encoding_scheme, applied_currents );

            % Get the applied current IDs that comprise this split voltage based integration subnetwork.
            applied_current_IDs = self.to_neuron_IDs2applied_current_IDs( neuron_IDs, applied_currents, undetected_option );
            
            % Retrieve the applied current IDs associated with the two operations.
            applied_current_IDs_Ias1 = applied_current_IDs( 1:2 );
            applied_current_IDs_Ias2 = applied_current_IDs( 3 );

            % Compute the applied current magnitudes associated with this subnetwork.
            [ Ias1, applied_currents, applied_current_manager ] = self.compute_svbi_Ias1( applied_current_IDs_Ias1, parameters_Ias1, encoding_scheme, applied_currents, true, undetected_option );
            [ Ias2, applied_currents, applied_current_manager ] = applied_current_manager.compute_svbi_Ias2( applied_current_IDs_Ias2, parameters_Ias2, encoding_scheme, applied_currents, true, undetected_option );
            
            % Concatenate the applied current magnitudes.
            Ias = [ Ias1, Ias2 ];
            
            % Determine whether to update the applied current manager.
            if set_flag, self = applied_current_manager; end      

        end
        
        
        %{
        
        % ---------- Centering Subnetwork Functions ----------

        % Implement a function to design the applied currents for a centering subnetwork.
        function [ Ias, applied_currents, self ] = design_centering_applied_current( self, neuron_IDs, Gm2, R2, applied_currents, set_flag, undetected_option )
            
            % Compute the number of subtraction neurons.
            n_neurons = self.num_centering_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 2, neuron_IDs = 1:n_neurons; end                                                         % [#] Neuron IDs
            
            % Get the applied currents IDs associated with the provided neuron IDs.
            applied_current_ID = self.to_neuron_ID2applied_current_ID( neuron_IDs( 2 ), applied_currents, 'ignore' );
            
            % Compute the applied cuyrrent magnitudes for this subnetwork.
            [ Ias, applied_currents, self ] = self.compute_centering_Ias( applied_current_ID, Gm2, R2, applied_currents, set_flag, undetected_option );
            
        end
        
        %}
        
        
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to design the applied currents for a multistate cpg subnetwork.
        function [ ts, Ias, applied_currents, self ] = design_mcpg_applied_current( self, neuron_IDs, dt, tf, applied_currents, filter_disabled_flag, set_flag, process_option, undetected_option )
            
            % Define the number of mcpg neurons.
            n_neurons = self.num_mcpg_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 9, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 8, process_option = self.process_option_DEFAULT; end                                    % [str] Process Option. (Must be either 'max', 'min', 'mean', or 'none'.)
            if nargin < 7, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 6, filter_disabled_flag = self.filter_disabled_flag_DEFAULT; end            % [T/F] Filter Disabled Flag. (Determines whether to considered disabled applied currents.)  
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, tf = self.tf_DEFAULT; end                                                                % [s] Simulation Duration
            if nargin < 3, dt = self.dt_DEFAULT; end                                                                % [s] Simulation Step Size
            if nargin < 2, neuron_IDs = 1:n_neurons; end                                                         % [#] Neuron IDs
            
            % Retrieve the applied current ID associated with the given final neuron ID.
            applied_current_ID = self.to_neuron_ID2applied_current_ID( neuron_IDs( end ), applied_currents, undetected_option );
            
            % Create an instance of the applied current manager.
            applied_current_manager = self;
            
            % Compute and set the applied current time vector.
            [ ts, applied_currents, applied_current_manager ] = applied_current_manager.compute_mcpg_ts( applied_current_ID, dt, tf, applied_currents, filter_disabled_flag, true, process_option, undetected_option );
            
            % Set the applied current magnitude vector.
            [ Ias, applied_currents, applied_current_manager ] = applied_current_manager.compute_mcpg_Ias( applied_current_ID, dt, tf, applied_currents, filter_disabled_flag, true, process_option, undetected_option );
            
            % Determine whether to update the applied current manager object.
            if set_flag, self = applied_current_manager; end
            
        end
        
        
        % Implement a function to design the applied currents for a driven multistate cpg subnetwork.
        function [ ts, Ias, applied_currents, self ] = design_dmcpg_applied_current( self, neuron_IDs, Gm, R, dt, tf, applied_currents, filter_disabled_flag, set_flag, process_option, undetected_option )
            
            % Define the number of dmcpg neurons.
            n_dmcpg_neurons = self.num_dmcpg_neurons_DEFAULT;
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, R = self.R_DEFAULT; end                                                                  % [V] Activation Domain
            if nargin < 3, Gm = self.Gm_DEFAULT; end                                                                % [S] Membrane Conductance
            if nargin < 2, neuron_IDs = 1:n_dmcpg_neurons; end                                                         % [#] Neuron IDs
            
            % Create an instance of the applied current manager.
            applied_current_manager = self;
            
            % Design the mcpg applied currents.
            [ ts, Ias_mcpg, applied_currents, applied_current_manager ] = applied_current_manager.design_mcpg_applied_current( neuron_IDs( 1:( end - 1 ) ), dt, tf, applied_currents, filter_disabled_flag, true, process_option, undetected_option );
            
            % Retrieve the applied current ID associated with the neuron ID.
            applied_current_ID = self.to_neuron_ID2applied_current_ID( neuron_IDs( end ), applied_currents, undetected_option );
            
            % Set the applied current magnitude vector.
            [ Ias_drive, applied_currents, applied_current_manager ] = applied_current_manager.compute_dmcpg_Ias( applied_current_ID, Gm, R, applied_currents, true, undetected_option );
            
            % Compute the number of mcpg time steps.
            n_timesteps = length( ts );
            
            % Determine whether to augment the drive current magnitudes.
            if size( Ias_drive, 1 ) == 1                              % If the drive current magnitude is a scalar with respect to the temporal dimension...
               
                % Repeat the drive current for each timestep.
                Ias_drive = repmat( Ias_drive, [ n_timesteps, 1 ] );
                
            end
            
            % Concatenate the applied current magnitudes.
            Ias = [ Ias_mcpg, Ias_drive ];
            
            % Determine whether to update the applied current manager object.
            if set_flag, self = applied_current_manager; end
            
        end
        
        
        % Implement a function to design the applied currents that connect a driven multistate cpg double centering lead lag subnetwork to a centered double subtraction subnetwork.
        function [ Ias, applied_currents, self ] = design_dmcpgdcll2cds_applied_current( self, neuron_ID, Gm, R, applied_currents, set_flag, undetected_option )
            
            % Set the default input arguments.
            if nargin < 7, undetected_option = self.undetected_option_DEFAULT; end          % [str] Undetected Option (Determines what to do if neuron ID is not detected.)
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 5, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 4, R = self.R_DEFAULT; end                                                                  % [V] Activation Domain
            if nargin < 3, Gm = self.Gm_DEFAULT; end                                                                % [S] Membrane Conductance
            if nargin < 2, neuron_ID = 1; end                                                         % [#] Neuron ID.
            
            % Retrieve the applied current ID associated with the neuron ID.
            applied_current_ID = self.to_neuron_ID2applied_current_ID( neuron_ID, applied_currents, undetected_option );
            
            % Set the applied current magnitude vector.
            [ Ias, applied_currents, self ] = self.compute_dmcpgdcll2cds_Ias( applied_current_ID, Gm, R, applied_currents, set_flag, undetected_option );
            
        end
        
        
        %% Print Functions.
        
        % Implement a function to print the properties of the constituent applied currents.
        
                
        %% Save & Load Functions.
        
        % Implement a function to save applied current manager data as a matlab object.
        function save( self, directory, file_name )
            
            % Set the default input arguments.
            if nargin < 3, file_name = 'Applied_Current_Manager.mat'; end
            if nargin < 2, directory = '.'; end
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Save the neuron data.
            save( full_path, self )
            
        end
        
        
        % Implement a function to load applied current manager data as a matlab object.
        function [ data, self ] = load( self, directory, file_name, set_flag )
            
            % Set the default input arguments.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag. (Determines whether to updated the applied current manager.)
            if nargin < 3, file_name = 'Applied_Current_Manager.mat'; end
            if nargin < 2, directory = '.'; end
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Load the data.
            data = load( full_path );
            
            % Determine whether to update the applied current manager object.
            if set_flag, self = data; end
            
        end
        
        
        % Implement a function to load applied current data from an xlsx file.
        function [ applied_currents, self ] = load_xlsx( self, file_name, directory, append_flag, verbose_flag, applied_currents, set_flag, data_loader_utilities )
            
            % Set the default input arguments.
            if nargin < 8, data_loader_utilities = self.data_loader_utilities; end          % [class] Data Load Utilities Class.
            if nargin < 7, set_flag  = self.set_flag; end                                   % [T/F] Set Flag (Determines whether output self object is updated.)
            if nargin < 6, applied_currents = self.applied_currents; end                            % [class] Array of Applied Current Class Objects.
            if nargin < 5, verbose_flag = true; end
            if nargin < 4, append_flag = false; end
            if nargin < 3, directory = '.'; end
            if nargin < 2, file_name = 'Applied_Current_Data.xlsx'; end
            
            % Determine whether to print status messages.
            if verbose_flag, fprintf( 'LOADING APPLIED CURRENT DATA. Please Wait...\n' ), end
            
            % Start a timer.
            tic
            
            % Load the applied current data.
            [ applied_current_IDs, applied_current_names, applied_current_to_neuron_IDs, applied_current_ts, applied_current_Ias ] = data_loader_utilities.load_applied_current_data( file_name, directory );
            
            % Define the number of applied currents.
            num_applied_currents_to_load = length( applied_current_IDs );
            
            % Preallocate an array of applied currents.
            applied_currents_to_load = repmat( applied_current_class(  ), 1, num_applied_currents_to_load );
            
            % Create each applied current object.
            for k = 1:num_applied_currents_to_load               % Iterate through each of the applied currents...
                
                % Create this applied current.
                applied_currents_to_load( k ) = applied_current_class( applied_current_IDs( k ), applied_current_names{ k }, applied_current_to_neuron_IDs( k ), applied_current_ts( :, k ), applied_current_Ias( :, k ) );
                
            end
            
            % Determine whether to append the applied currents we just loaded.
            if append_flag                         % If we want to append the applied currents we just loaded...
                
                % Append the applied currents we just loaded to the array of existing applied currents.
                applied_currents = [ applied_currents, applied_currents_to_load ];
                
                % Update the number of applied currents.
                n_applied_currents = length( applied_currents );
                
            else                                % Otherwise...
                
                % Replace the existing applied currents with the applied currents we just loaded.
                applied_currents = applied_currents_to_load;
                
                % Update the number of applied currents.
                n_applied_currents = length( applied_currents );
                
            end
            
            % Determine whether to update the applied current manager properties.
            if set_flag                                             	% If we want to update the applied current manager properties...
                
                % Update the applied currents property.
                self.applied_currents = applied_currents;
                
                % Update the number of applied currents.
                self.num_applied_currents = n_applied_currents;
                
            end
            
            % Retrieve the elapsed time.
            elapsed_time = toc;
            
            % Determine whether to print status messages.
            if verbose_flag, fprintf( 'LOADING APPLIED CURRENT DATA. Please Wait... Done. %0.3f [s] \n\n', elapsed_time ), end
            
        end
        
        
    end
end

