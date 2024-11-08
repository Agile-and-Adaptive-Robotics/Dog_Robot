classdef array_utilities_class
    
    % This class constains properties and methods related to converting data types.
    
    
    %% ARRAY MANAGER PROPERTIES
    
    % Define the class properties.
    properties ( Constant = true )
        
        
        
    end
    
    
    %% ARRAY MANAGER SETUP FUNCTIONS
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = array_utilities_class(  )
            
            
        end
        
        
        %% Array Functions
        
        % Implement a function to remove all specified values from an array.
        function array = remove_entries( ~, array, entries_to_remove )
            
            % Retrieve the number of entries to remove.
            num_entries_to_remove = length( entries_to_remove );
            
            % Remove each of the specified entries.
            for k = 1:num_entries_to_remove                          % Iterate through each of the entries ot remove...
                
                % Remove this entry.
                array( array == entries_to_remove( k ) ) = [  ];
                
            end
            
        end
        
        
        % Implement a function to determine whether a value exists in an array.
        function [ match_found_flag, match_logicals, match_indexes ] = is_value_in_array( ~, value, array )
            
            % Determine how to process the inputs.
            if ( ~isempty( value ) ) && ( ~isempty( array ) )                           % If neither the value nor the array are empty...
            
                % Retrieve the number of array values.
                num_array_values = length( array );

                % Determine which array entries match the value.
                match_logicals = value == array;

                % Define an array of indexes.
                indexes = 1:num_array_values;

                % Retrieve the indexes of the array entries that match the value.
                match_indexes = indexes( match_logicals );

                % Determine whether the value is contained in the array.
                match_found_flag = any( match_logicals );
            
            else                                                                            % Otherwise...
                
                % Set the match found flag to false.
                match_found_flag = false;
                
                % Set the match logicals and match indexes to be empty.
                [ match_logicals, match_indexes ] = deal( [  ] );
                
            end
                
        end
        
        
        % Implement a function to determine whether any of an array of values is contained in a given array.
        function [ match_found_flags, match_logicals, match_indexes ] = are_values_in_array( self, values, array )
            
            % Retrieve the number of values to check.
            num_values = length( values );
            
            % Retrieve the number of array values.
            num_array_values = length( array );
            
            % Preallocate arrays to store the uniqueness flags, match logicals, and match indexes.
            match_found_flags = false( num_values, 1 );
            [ match_logicals, match_indexes ] = deal( zeros( num_values, num_array_values ) );
            
            % Determine whether each value is in the array.
            for k = 1:num_values                % Iterate through each of the values...
                
                % Determine whether this value is in the array.
                [ match_found_flags( k ), match_logicals( k, : ), match_indexes( k, : ) ] = self.is_value_in_array( values( k ), array );
                
            end
            
        end
        
        
        % Implement a function to get the lowest whole number not in an array.
        function n = get_lowest_natural_number( self, exclude_array )
            
            % Determine the lowest natural number (excluding the array of 
            if nargin < 2                           % If no values were given to exclude...
                
                % Set the lowest natural number to one.
                n = 1;
                
            else                                    % Otherwise...
                           
                % Initialize the unique flag to false.
                b_unique = false;
                
                % Initialize the lowest natural number to be zero.
                n = 0;
                
                % Compute the lowest natural number.
                while ~b_unique                                 % While the natural number that we have selected is not unique...
                    
                    % Advance the candidate natural number.
                    n = n + 1;
                    
                    % Determine whether this natural number candidate is in the exclude array.
                    if ~self.is_value_in_array( n, exclude_array )                  % If this natural number candidate is not in the exclude array...
                        
                        % Set the unique flag to true.
                        b_unique = true;
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % Implement a function to compute the step size of an array.
        function step_size = compute_step_size( ~, array )
            
            % Determine how to compute the step size of the array.
           if length( array ) > 1                                       % If there are at least two array elements...
               
               % Compute the step size.
               step_size = mean( array( 2:end ) - array( 1:end-1 ) );
               
           else                                                         % Otherwise...
              
               % Set the step size to be empty.
               step_size = [  ];
               
           end
            
        end
        
        
        % Implement a function to scale a domain in an absolute sense.
        function domain = scale_domain_absolute( ~, domain, scale )
            
            % Determine whether to directly add the scale.
            if domain( 1 ) == 0                                           % If this bound is zero...
                
                % Directly add the scale to the bound.
                domain( 1 ) = domain( 1 ) - scale;
                
            else                                                        % Otherwise...
                
                % Scale the domain lower bound.
                domain( 1 ) = domain( 1 ) - abs( scale*domain( 1 ) );
                
            end
            
            % Determine whether to directly add the scale.
            if domain( 2 ) == 0                                           % If this bound is zero...
                
                % Directly add the scale to the bound.
                domain( 2 ) = domain( 2 ) + scale;
                
            else                                                        % Otherwise...
                
                % Scale the domain upper bound.
                domain( 2 ) = domain( 2 ) + abs( scale*domain( 2 ) );
                
            end

        end
        
        
        % Implement a function to scale a domain in a relative sense.
        function domain = scale_domain_relative( ~, domain, scale )
           
            % Compute the middle of the domain.
            middle = ( domain( 2 ) - domain( 1 ) )/2 + domain( 1 );

            % Center the domain at the origin.
            domain_centered = domain - middle;
            
            % Scale the domain.
            domain_scaled = scale*domain_centered;
            
            % Uncenter the scaled domain.
            domain = domain_scaled + middle; 
            
        end
        
        
        % Implement a function to expand a given domain.
        function domain = scale_domain( self, domain, scale, type )
        
            % Define the default input arguments.
            if nargin < 4, type = 'relative'; end
            
            % Determine how to perform the scaling.
            if strcmpi( type, 'relative' )                                      % If the scaling type is relative...
                
                % Scale the domain in a relative sense.
                domain = self.scale_domain_relative( domain, scale );
                
            elseif strcmpi( type, 'absolute' )                                  % If the scaling type is absolute...
                
                % Scale the domain in an absolute sense.
                domain = self.scale_domain_absolute( domain, scale );
                
            else                                                                % Otherwise...
                
                % Throw an error.
                error( 'Unrecognized scaling type: %s.', type )
                
            end
            
        end
            
        
        % Implement a function to return the nearest lower value in an array.
        function [ x_low_nearest, index_low_nearest ] = get_lower_nearest_value( ~, xs, xs_sample )
            
            % Compute the difference between the provided array and the sample.
            dxs = xs - xs_sample;
            
            % Create a copy of the difference array.
            dxs_neg = dxs;
            
            % Set the positive difference value to be the maximum value in the array (this ensures that they are ignored by the upcoming min call).
            dxs_neg( dxs > 0 ) = max( dxs );
            
            % Compute the absolute value of the difference.
            dxs_neg = abs( dxs_neg );
            
            % Retrieve the index of the nearest lower value in the array.
            [ ~, index_low_nearest ] = min( dxs_neg );
            
            % Retrieve the nearest lower value in the array.
            x_low_nearest = xs( index_low_nearest );
            
        end
        
        
        % Implement a function to return the nearest lower value in an array.
        function [ x_high_nearest, index_high_nearest ] = get_higher_nearest_value( ~, xs, xs_sample )
            
            % Compute the difference between the provided array and the sample.
            dxs = xs - xs_sample;
            
            % Create a copy of the difference array.
            dxs_pos = dxs;
            
            % Set the negative difference values to be the maximum value in the array (this ensures that they are ignored by the upcoming min call).
            dxs_pos( dxs < 0 ) = max( dxs );
            
            % Retrieve the index of the nearest lower value in the array.
            [ ~, index_high_nearest ] = min( dxs_pos );
            
            % Retrieve the nearest higher value in the array.
            x_high_nearest = xs( index_high_nearest );
            
        end
        
        
        % Implement a function to return the nearest value in an array.
        function [ xs_nearest, indexes_nearest ] = get_nearest_value( self, xs, xs_sample, direction )
            
            % Set the default input arguments.
            if nargin < 4, direction = 'Both'; end
            
            % Determine how to return the nearest values.
            if strcmpi( direction, 'left' ) || strcmpi( direction, 'down' ) || strcmpi( direction, 'floor' )
                
                % Compute the nearest lower value.
                [ xs_nearest, indexes_nearest ] = self.get_lower_nearest_value( xs, xs_sample );
                
            elseif strcmpi( direction, 'right' ) || strcmpi( direction, 'up' ) || strcmpi( direction, 'ceil' )
            
                % Compute the nearest higher value.
                [ xs_nearest, indexes_nearest ] = self.get_higher_nearest_value( xs, xs_sample );
                
            else
                
                % Compute the nearest lower value.
                [ x_low_nearest, index_low_nearest ] = self.get_lower_nearest_value( xs, xs_sample );
                
                % Compute the nearest higher value.
                [ x_high_nearest, index_high_nearest ] = self.get_higher_nearest_value( xs, xs_sample );
               
                % Concatenate the nearest values.
                indexes_nearest = [ index_low_nearest, index_high_nearest ];
                xs_nearest = [ x_low_nearest, x_high_nearest ];
                
            end
            
        end
        
        
        %% Cell Array Functions.
        
        % Implement a function to determine whether any of the entries in a cell array are themselves cells.
        function cell_present = contains_cell( ~, array )
            
            %{
            Input(s):
                array           =   [variable] Cell array that may itself contain cells.
            
            Output(s):
                cell_present    =   [T/F] Flag that describes whether the provided array contains cells.
            %}
            
            % Initialize the cell present flag to false.
            cell_present = false;
            
            % Determine whether the input array is in fact a cell.
            if iscell( array )                                      % If the provided array is a cell...
            
                % Compute the number of array entries.
                n_entries = length( array );
               
                % Initialize a loop counter variable.
                k = 0;

                % Determine whether
                while ( k < n_entries ) && ~cell_present             % While we have not yet checked every entry and have not yet found a cell...

                    % Advance the loop counter variable.
                    k = k + 1;

                    % Determine whetheer this entry is a cell.
                    cell_present = iscell( array{ k } );

                end
                
            end
            
        end
        
        
        % Implement a function to convert a cell array to a regular array.
        function array = cell2array( self, array )
            
            % Extract the contents of the array as long as it still contains
            while self.contains_cell( array )               % While the provided array contains at least one cell entry...
            
                % Extract and concatenate the contents of the cell entries.
                array = [ array{ : } ];
            
            end
            
        end
        
        
        % Implement a function to interpolate a cell array that may contain empty entries.
        function ys_cell_sample = interp1_cell( self, xs_cell, ys_cell, xs_cell_sample )
        
            % Determine whether it is necessary to convert the reference and sample inputs to matrices.
            if isa( xs_cell, 'cell' ), xs = cell2mat( xs_cell ); else, xs = xs_cell; end
            if isa( xs_cell_sample, 'cell' ), xs_sample = cell2mat( xs_cell_sample ); else, xs_sample = xs_cell_sample; end
            
            % Compute the number of samples to evaluate.
            num_samples = length( xs_sample );
            
            % Preallocate an array to store the sample outputs.
            ys_cell_sample = cell( size( xs_sample ) );
            
            % Ensure that the reference inputs are unique.
            assert( all( unique( xs ) == xs ), 'xs entries must be unique.')
            
            % Compute the output value at each of the sample inputs.
            for k = 1:num_samples                % Iterate through each of the sample inputs...
                
                % Determine how to compute the sample output.
                if ( xs_sample( k ) < min( xs ) ) || ( xs_sample( k ) > max( xs ) )             % If this sample input is output of bounds...
                    
                    % Set the sample output to be empty.
                    ys_cell_sample( k ) = { [  ] };
                    
                elseif any( xs_sample( k ) == xs )                                              % If this sample input is equal to one of the reference inputs...
                    
                    % Set the sample output to be the reference output.
                    ys_cell_sample( k ) = ys_cell( xs_sample( k ) == xs );
                    
                else                                                                            % Otherwise...
                    
                    % Compute the nearest reference inputs.
                    [ ts_nearest, indexes_nearest ] = self.get_nearest_value( xs, xs_sample( k ) );
                    
                    % Compute the nearest reference outputs.
                    ys_cell_nearest = ys_cell( indexes_nearest );
                    
                    % Determine how to set the sample output.
                    if isempty( ys_cell_nearest{ 1 } ) || isempty( ys_cell_nearest{ 2 } )             % If either of the nearest reference outputs are empty...
                        
                        % Set the sample output to be empty.
                        ys_cell_sample( k ) = { [  ] };
                        
                    else                                                                            % Otherwise...
                        
                        % Compute the sample output via interpolation.
                        ys_cell_sample( k ) = num2cell( interp1( ts_nearest, cell2mat( ys_cell( indexes_nearest ) ), xs_sample( k ) ) );
                        
                    end
                    
                end
                
            end
            
        end
        
    end
end

