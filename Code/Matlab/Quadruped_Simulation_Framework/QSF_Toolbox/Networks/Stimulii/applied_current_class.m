classdef applied_current_class
    
    % This class contains properties and methods related to applied currents.
    
    
    %% APPLIED CURRENT PROPERTIES
    
    % Define the class properties.
    properties
        
        ID                                              % [#] Applied Current ID.
        name                                            % [str] Applied Current Name.
        to_neuron_ID                                    % [#] To Neuron ID.
        
        ts                                              % [s] Time Vector.
        Ias                                             % [A] Applied Current Vector.
        
        num_timesteps                                   % [#] Number of Timesteps.
        dt                                              % [s] Time Step Duration.
        tf                                              % [s] Simulation Duration.
        
        enabled_flag                                    % [T/F] Enabled Flag.
        
        array_utilities                                 % [class] Array Utilities.
        applied_current_utilities                       % [class] Applied Current Utilities.
        
    end
    
    
    % Define private, constant class properties.
    properties ( Access = private, Constant = true )
        
        % Define the neuron parameters.
        R_DEFAULT = 20e-3;                              % [V] Activation Domain.
        Gm_DEFAULT = 1e-6;                              % [S] Membrane Conductance.
        
        % Define the applied current properties.
        Ia_DEFAULT = 0e-9;                              % [A] Applied Current Magnitudes.
        ts_DEFAULT = 0;                                 % [s] Applied Current Times.
        
        % Set the default encoding scheme.
        encoding_scheme_DEFAULT = 'absolute';           % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
        
        % Set the default flag values.
        enabled_flag_DEFAULT = true;                    % [T/F] Enabled Flag (Determines whether this applied current is active during simulations). 
        
    end
    
    
    %% APPLIED CURRENT METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = applied_current_class( ID, name, to_neuron_ID, ts, Ias, enabled_flag, array_utilities, applied_current_utilities )
            
            % Set the default properties.
            if nargin < 8, applied_current_utilities = applied_current_utilities_class(  ); end             % [class] Applied Current Utilities.
            if nargin < 7, array_utilities = array_utitlies_class(  ); end                                  % [class] Array Utilities.
            if nargin < 6, enabled_flag = self.enabled_flag_DEFAULT; end                                    % [T/F] Applied Current Enabled Flag.
            if nargin < 5, Ias = self.Ia_DEFAULT; end                                                       % [A] Applied Current Magnitudes.
            if nargin < 4, ts = self.ts_DEFAULT; end                                                        % [s] Applied Current Times.
            if nargin < 3, to_neuron_ID = 0; end                                                         	% [#] Application Neuron ID.
            if nargin < 2, name = ''; end                                                                   % [-] Applied Current Name.
            if nargin < 1, ID = 0; end                                                                      % [#] Applied Current ID.
            
            % Store an instance of the utility classes.
            self.applied_current_utilities = applied_current_utilities;
            self.array_utilities = array_utilities;
            
            % Store the applied current flags.
            self.enabled_flag = enabled_flag;
            
            % Store the applied current properties.
            self.Ias = Ias;
            self.ts = ts;
            
            % Store the applied current information.
            self.to_neuron_ID = to_neuron_ID;
            self.name = name;
            self.ID = ID;
            
            % Validate the applied current.
            assert( self.is_applied_current_valid( ts, Ias ), 'The lengths of the time vector and applied current vectors must be equal.' )

            % Compute the number of timesteps.
            [ ~, self ] = self.compute_num_timesteps( ts, true );
            
            % Compute the step size.
            [ ~, self ] = self.compute_dt( ts, true, array_utilities );
                        
            % Compute the final time.
            [ ~, self ] = self.compute_tf( ts, true );
            
        end
        

        %% Name Functions.
        
        % Implement a function to generate a name for this applied current.
        function [ name, self ] = generate_name( self, ID, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end
            if nargin < 2, ID = self.ID; end
            
            % Generate a name for the applied current.
            name = applied_current_utilities.ID2name( ID );
            
            % Determine whether to update the name.
            if set_flag, self.name = name; end
            
        end
        
        
        %% Validation Functions.
        
        % Implement a function to validate the applied current.
        function valid_flag = is_applied_current_valid( self, ts, Ias )
           
            % Set the default input arguments.
            if nargin < 3, Ias = self.Ias; end          % [A] Applied Current Magnitudes.
            if nargin < 2, ts = self.ts; end            % [s] Applied Current Times.
            
            % Validate the applied current.
            valid_flag = length( ts ) == length( Ias );
            
        end
        
        
        %% Compute Time Functions.
        
        % Implement a function to compute the number of time steps.
        function [ n_timesteps, self ] = compute_num_timesteps( self, ts, set_flag )
           
            % Set the default input arguments.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                  	% [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, ts = self.ts; end                                                        % [s] Applied Current Times.
            
            % Compute the number of time steps.
            n_timesteps = length( ts );                                                             % [#] Number of Simulation Timesteps.
            
            % Determine whether to update the applied current object.
            if set_flag, self.num_timesteps = n_timesteps; end
            
        end
        
        
        % Implement a function to compute the time step.
        function [ dt, self ] = compute_dt( self, ts, set_flag, array_utilities )
            
            % Set the default input arguments.
            if nargin < 4, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, ts = self.ts; end                                                        % [s] Applied Current Times.
            
            % Compute the step size.
            dt = array_utilities.compute_step_size( ts );                                           % [s] Simulation Time Step Size.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dt = dt; end
            
        end
        
        
        % Implement a function to compute the final time.
        function [ tf, self ]= compute_tf( self, ts, set_flag )
            
            % Set the default input arguments.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, ts = self.ts; end                                                        % [s] Applied Current Times.
            
            % Compute the final time.
            tf = max( ts );                                                                         % [s] Final Simulation Time.
            
           % Determine whether to update the synapse object.
            if set_flag, self.tf = tf; end
           
        end
        
        
        % Implement a function to compute the properties associated with an applied current vector.
        function [ n_timesteps, dt, tf, self ] = compute_applied_current_properties( self, ts, Ias, set_flag, array_utilities )
            
            % Set the default input arugments.
            if nargin < 5, array_utilities = self.array_utilities; end                              % [class] Array Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, Ias = self.Ia_DEFAULT; end                                               % [A] Applied Current Magnitudes.
            if nargin < 2, ts = self.ts_DEFAULT; end                                                % [s] Applied Current Times.
            
            % Ensure that there are the same number of time steps as applied currents.
            assert( self.is_applied_current_valid( ts, Ias ), 'The lengths of the time vector and applied current vectors must be equal.' )
            
            % Create an instance of the applied current object that can be updated.
            applied_current = self;
            
            % Update the applied current time and magnitude vectors.
            applied_current.ts = ts;
            applied_current.Ias = Ias;
            
            % Compute the number of timesteps.
            [ n_timesteps, applied_current ] = applied_current.compute_num_timesteps( ts, true );
            
            % Compute the step size.
            [ dt, applied_current ] = applied_current.compute_dt( ts, true, array_utilities );
                        
            % Compute the final time.
            [ tf, applied_current ] = applied_current.compute_tf( ts, true );
            
            % Determine whether to update the applied current object.
            if set_flag, self = applied_current; end
            
        end
        
        
        % Implement a function to compute the time vector of multistate cpg subnetwork applied currents.
        function [ ts, self ] = compute_mcpg_ts( self, dt, tf, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities.
            if anrgin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, tf = self.tf; end                                                        % [s] Final Simulation Time.
            if nargin < 2, dt = self.dt; end                                                        % [s] Simulation Timestep.
            
           % Compute the time vector of multistate cpg subnetwork applied currents.
           ts = applied_current_utilities.compute_mcpg_ts( dt, tf );                              	% [s] Applied Current Times.
            
           % Determine whether to update the applied current object.
           if set_flag, self.ts = ts; end
           
        end
        
        
        %% Parameter Unpacking Functions.
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 2 of an absolute inversion subnetwork.
        function [ Gm2, R2 ] = unpack_absolute_inversion_Ias2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                Gm2 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                R2 = self.R_DEFAULT;                                 	% [V] Activation Domain.
                
            elseif length( parameters ) == 2                           	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                Gm2 = parameters{ 1 };                                  % [S] Membrane Conductance.
                R2 = parameters{ 2 };                               	% [V] Activation Domain.
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
            
            end      
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 2 of a relative inversion subnetwork.
        function [ Gm2, R2 ] = unpack_relative_inversion_Ias2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                R2 = self.R_DEFAULT;                                    % [V] Activation Domain.
                Gm2 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                
            elseif length( parameters ) == 2                          	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                R2 = parameters{ 1 };                                 	% [V] Activation Domain.
                Gm2 = parameters{ 2 };                                  % [S] Membrane Conductance.
                
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end      
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 2 of a reduced absolute inversion subnetwork.
        function [ Gm2, R2 ] = unpack_reduced_absolute_inversion_Ias2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                Gm2 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                R2 = self.R_DEFAULT;                                 	% [V] Activation Domain.
                
            elseif length( parameters ) == 2                           	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                Gm2 = parameters{ 1 };                                  % [S] Membrane Conductance.
                R2 = parameters{ 2 };                               	% [V] Activation Domain.
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end      
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 2 of a reduced relative inversion subnetwork.
        function [ Gm2, R2 ] = unpack_reduced_relative_inversion_Ias2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                Gm2 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                R2 = self.R_DEFAULT;                                    % [V] Activation Domain.
                
            elseif length( parameters ) == 1                          	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                Gm2 = parameters{ 1 };                                  % [S] Membrane Conductance.
                R2 = parameters{ 2 };                                 	% [V] Activation Domain.
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end      
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------

        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 3 of an absolute multiplication subnetwork.
        function [ Gm3, R3 ] = unpack_absolute_multiplication_Ias3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                Gm3 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                R3 = self.R_DEFAULT;                                 	% [V] Activation Domain.
                
            elseif length( parameters ) == 2                           	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                Gm3 = parameters{ 1 };                                  % [S] Membrane Conductance.
                R3 = parameters{ 2 };                               	% [V] Activation Domain.
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end      
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 3 of a relative multiplication subnetwork.
        function [ Gm3, R3 ] = unpack_relative_multiplication_Ias3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                Gm3 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                R3 = self.R_DEFAULT;                                 	% [V] Activation Domain.
                
            elseif length( parameters ) == 2                           	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                Gm3 = parameters{ 1 };                                  % [S] Membrane Conductance.
                R3 = parameters{ 2 };                               	% [V] Activation Domain.
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end      
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 3 of a reduced absolute multiplication subnetwork.
        function [ Gm3, R3 ] = unpack_reduced_absolute_multiplication_Ias3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                Gm3 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                R3 = self.R_DEFAULT;                                 	% [V] Activation Domain.
                
            elseif length( parameters ) == 2                           	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                Gm3 = parameters{ 1 };                                  % [S] Membrane Conductance.
                R3 = parameters{ 2 };                               	% [V] Activation Domain.
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end      
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the applied current magnitudes for neuron 3 of a reduced relative multiplication subnetwork.
        function [ Gm3, R3 ] = unpack_reduced_relative_multiplication_Ias3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [cell] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                Gm3 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                R3 = self.R_DEFAULT;                                 	% [V] Activation Domain.
                
            elseif length( parameters ) == 2                           	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                Gm3 = parameters{ 1 };                                  % [S] Membrane Conductance.
                R3 = parameters{ 2 };                               	% [V] Activation Domain.
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end      
            
        end
        
        
        %% Compute Applied Current Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron 2 of a transmission subnetwork.
        function [ Ias2, self ] = compute_transmission_Ias2( self, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
                            
                % Compute the applied current magnitudes.
                Ias2 = applied_current_utilities.compute_absolute_transmission_Ias2(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Compute the applied current magnitudes.
                Ias2 = applied_current_utilities.compute_relative_transmission_Ias2(  );
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias2; end

        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron n of an addition subnetwork.
        function [ Ias_n, self ] = compute_addition_Iasn( self, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
                            
                % Compute the applied current magnitudes.
                Ias_n = applied_current_utilities.compute_absolute_addition_Iasn(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Compute the applied current magnitudes.
                Ias_n = applied_current_utilities.compute_relative_addition_Iasn(  );
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias_n; end

        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron n of a subtraction subnetwork.
        function [ Ias_n, self ] = compute_subtraction_Iasn( self, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
                            
                % Compute the applied current magnitudes.
                Ias_n = applied_current_utilities.compute_absolute_subtraction_Iasn(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Compute the applied current magnitudes.
                Ias_n = applied_current_utilities.compute_relative_subtraction_Iasn(  );
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias_n; end

        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron 2 of an inversion subnetwork.
        function [ Ias2, self ] = compute_inversion_Ias2( self, parameters, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            if nargin < 2, parameters = {  }; end                                                   % [cell] Parmaeters Cell.
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
            
                % Unpack the parameters required to compute the absolute inversion applied current magnitudes.
                [ Gm2, R2 ] = self.unpack_absolute_inversion_Ias2_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias2 = applied_current_utilities.compute_absolute_inversion_Ias2( Gm2, R2 );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Unpack the parameters required to compute the relative inversion applied current magnitudes.
                [ Gm2, R2 ] = self.unpack_relative_inversion_Ias2_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias2 = applied_current_utilities.compute_relative_inversion_Ias2( Gm2, R2 );  
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias2; end

        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to compute the applied current magnitude for neuron 2 of a reduced inversion subnetwork.
        function [ Ias2, self ] = compute_reduced_inversion_Ias2( self, parameters, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            if nargin < 2, parameters = {  }; end                                                   % [cell] Parmaeters Cell.
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
            
                % Unpack the parameters required to compute the absolute inversion applied current magnitudes.
                [ Gm2, R2 ] = self.unpack_reduced_absolute_inversion_Ias2_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias2 = applied_current_utilities.compute_reduced_absolute_inversion_Ias2( Gm2, R2 );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Unpack the parameters required to compute the relative inversion applied current magnitudes.
                [ Gm2, R2 ] = self.unpack_reduced_relative_inversion_Ias2_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias2 = applied_current_utilities.compute_reduced_relative_inversion_Ias2( Gm2, R2 );  
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias2; end

        end
        
        
        % ---------- Division Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron 3 of a division subnetwork.
        function [ Ias3, self ] = compute_division_Ias3( self, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
                            
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_absolute_division_Ias3(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_relative_division_Ias3(  );
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias3; end

        end
        
        
        % ---------- Reduced Division Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron 3 of a reduced division subnetwork.
        function [ Ias3, self ] = compute_reduced_division_Ias3( self, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
                            
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_reduced_absolute_division_Ias3(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_reduced_relative_division_Ias3(  );
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias3; end

        end
        
        
        % ---------- Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron 3 of a division after inversion subnetwork.
        function [ Ias3, self ] = compute_dai_Ias3( self, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
                            
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_absolute_dai_Ias3(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_relative_dai_Ias3(  );
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias3; end

        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron 3 of a reduced division after inversion subnetwork.
        function [ Ias3, self ] = compute_reduced_dai_Ias3( self, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 4, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
                            
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_reduced_absolute_dai_Ias3(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_reduced_relative_dai_Ias3(  );
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias3; end

        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the applied current magnitude for neuron 3 of a multiplication subnetwork.
        function [ Ias3, self ] = compute_multiplication_Ias3( self, parameters, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            if nargin < 2, parameters = {  }; end                                                   % [cell] Parmaeters Cell.
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
            
                % Unpack the parameters required to compute the absolute multiplication applied current magnitudes.
                [ Gm3, R3 ] = self.unpack_absolute_multiplication_Ias3_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_absolute_multiplication_Ias3( Gm3, R3 );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Unpack the parameters required to compute the relative multiplication applied current magnitudes.
                [ Gm3, R3 ] = self.unpack_relative_multiplication_Ias3_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_relative_multiplication_Ias3( Gm3, R3 );  
                
            else                                                                                  	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias3; end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to compute the applied current magnitude for neuron 3 of a reduced multiplication subnetwork.
        function [ Ias3, self ] = compute_reduced_multiplication_Ias3( self, parameters, encoding_scheme, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'absolute' or 'relative'.)
            if nargin < 2, parameters = {  }; end                                                   % [cell] Parmaeters Cell.
            
            % Determine how to compute the applied current magnitude.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is absolute...
            
                % Unpack the parameters required to compute the absolute multiplication applied current magnitudes.
                [ Gm3, R3 ] = self.unpack_reduced_absolute_multiplication_Ias3_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_reduced_absolute_multiplication_Ias3( Gm3, R3 );
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is relative...
               
                % Unpack the parameters required to compute the relative multiplication applied current magnitudes.
                [ Gm3, R3 ] = self.unpack_reduced_relative_multiplication_Ias3_parameters( parameters );
                
                % Compute the applied current magnitudes.
                Ias3 = applied_current_utilities.compute_reduced_relative_multiplication_Ias3( Gm3, R3 );  
                
            else                                                                                  	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias3; end
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of integration subnetwork applied currents.
        function [ Ias, self ] = compute_integration_Ias( self, Gm, R, set_flag, applied_current_utilities )
        
            % Define the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, R = self.R_DEFAULT; end                                                  % [V] Activation Domain.
            if nargin < 2, Gm = self.Gm_DEFAULT; end                                                % [S] Membrane Conductance.
            
            % Compute the magnitude of integration subnetwork applied currents.
            Ias = applied_current_utilities.compute_integration_Iapps( Gm, R );
           
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias; end
            
        end
            
        
        % Implement a function to compute the magnitude of voltage based integration subnetwork applied currents.
        function [ Ias, self ] = compute_vbi_Ias( self, Gm, R, set_flag, applied_current_utilities )
        
            % Define the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, R = self.R_DEFAULT; end                                                  % [V] Activation Domain.
            if nargin < 2, Gm = self.Gm_DEFAULT; end                                                % [S] Membrane Conductance.
            
            % Compute the magnitude of voltage based integration subnetwork applied currents.
            Ias = applied_current_utilities.compute_vb_integration_Iapps( Gm, R );
            
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias; end
            
        end
        
        
        % Implement a function to compute the first magnitude of split voltage based integration subnetwork applied currents.
        function [ Ias, self ] = compute_svbi_Ias1( self, Gm, R, set_flag, applied_current_utilities )
        
            % Define the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                  	% [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, R = self.R_DEFAULT; end                                                	% [V] Activation Domain.
            if nargin < 2, Gm = self.Gm_DEFAULT; end                                               	% [S] Membrane Conductance.
            
            % Compute the first magnitude of split voltage based integration subnetwork applied currents.
            Ias = applied_current_utilities.compute_svbi_Ias1( Gm, R );                             % [A] Applied Current.
            
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias; end
            
        end
        
        
        % Implement a function to compute the second magnitude of split voltage based integration subnetwork applied currents.
        function [ Ias, self ] = compute_svbi_Ias2( self, Gm, R, set_flag, applied_current_utilities )
        
            % Define the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, R = self.R_DEFAULT; end                                               	% [V] Activation Domain.
            if nargin < 2, Gm = self.Gm_DEFAULT; end                                              	% [S] Membrane Conductance.
            
            % Compute the second magnitude of split voltage based integration subnetwork applied currents.
            Ias = applied_current_utilities.compute_svbi_Ias2( Gm, R );                             % [A] Applied Current.
            
            % Determine whether to update the applied current object.
            if set_flag, self.Ias = Ias; end
            
        end
        
        
        % ---------- Centering Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of centering subnetwork applied currents.
        function [ Ias, self ] = compute_centering_Ias( self, Gm, R, set_flag, applied_current_utilities )
            
            % Define the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, R = self.R_DEFAULT; end                                                  % [V] Activation Domain.
            if nargin < 2, Gm = self.Gm_DEFAULT; end                                                % [S] Membrane Conductance.
            
            % Compute the magnitude of centering subnetwork applied currents.
            Ias = applied_current_utilities.compute_centering_Ias( Gm, R );
            
            % Determine whether to update the applied current object.
           if set_flag, self.Ias = Ias; end
            
        end
        
        
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to compute the magnitude of multistate cpg subnetwork applied currents.
        function [ Ias, self ] = compute_mcpg_Ias( self, dt, tf, set_flag, applied_current_utilities )
            
            % Set the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, tf = self.tf; end                                                        % [s] Final Simulation Time.
            if nargin < 2, dt = self.dt; end                                                        % [s] Simulation Timestep.
            
           % Compute the magnitude of multistate cpg subnetwork applied currents.
           Ias = applied_current_utilities.compute_mcpg_Ias( dt, tf );
            
           % Determine whether to update the applied current object.
           if set_flag, self.Ias = Ias; end
           
        end
        
        
        % Implement a function to compute the magnitude of driven multistate cpg subnetwork applied currents.
        function [ Ias, self ] = compute_dmcpg_Ias( self, Gm, R, set_flag, applied_current_utilities )
            
            % Define the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, R = self.R_DEFAULT; end                                                  % [V] Activation Domain.
            if nargin < 2, Gm = self.Gm_DEFAULT; end                                                % [S] Membrane Conductance.
            
           % Compute the magnitude of driven multistate cpg subnetwork applied currents.
           Ias = applied_current_utilities.compute_dmcpg_Ias( Gm, R );
            
           % Determine whether to update the applied current object.
           if set_flag, self.Ias = Ias; end
           
        end
        
        
        % Implement a function to compute the applied current magnitude that connects the dmcpgdcll and cds subnetworks.
        function [ Ias, self ] = compute_dmcpgdcll2cds_Ias( self, Gm, R, set_flag, applied_current_utilities )
            
            % Define the default input arguments.
            if nargin < 5, applied_current_utilities = self.applied_current_utilities; end          % [class] Applied Current Utilities Class.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the applied current manager.)
            if nargin < 3, R = self.R_DEFAULT; end                                                  % [V] Activation Domain.
            if nargin < 2, Gm = self.Gm_DEFAULT; end                                                % [S] Membrane Conductance.
            
           % Compute the magnitude of these applied currents.
           Ias = applied_current_utilities.compute_dmcpgdcll2cds_Ias( Gm, R );
            
           % Determine whether to update the applied current object.
           if set_flag, self.Ias = Ias; end
           
        end
                
        
        %% Sampling Functions.
        
        % Implement a function to sample an applied current.
        function [ ts_sample, Ias_sample ] = sample_Ias( self, dt_sample, tf_sample, ts_reference, Ias_reference )
            
            % Set the default input arguments.
            if nargin < 5, Ias_reference = self.Ias; end                        % [A] Applied Current Magnitudes.
            if nargin < 4, ts_reference = self.ts; end                          % [s] Applied Current Times.
            if nargin < 3, tf_sample = self.tf; end                             % [s] Final Simulation Time.
            if nargin < 2, dt_sample = self.dt; end                             % [s] Simulation Timestep.
            
            % Compute the number of reference timesteps.
            n_timesteps_reference = length( Ias_reference );
            
            % Determine how to sample the applied current.
            if ~isempty( dt_sample ) && ( n_timesteps_reference > 0 )           % If the sample spacing and existing applied current data is not empty...
                
                % Create the sampled time vector.
                ts_sample = ( 0:dt_sample:tf_sample )';
                
                % Determine how to sample the applied current.
                if n_timesteps_reference == 1                               	% If the number of timesteps is one.. ( The applied current is constant. )
                    
                    % Create the applied current sample as a constant vector.
                    Ias_sample = Ias_reference*ones( n_timesteps_reference, 1 );
                    
                else                                                        	% Otherwise...
                    
                    % Create the applied current sample via interpolation.
                    Ias_sample = interp1( ts_reference, Ias_reference, ts_sample );
                    
                end
                
            else                                                             	% Otherwise...
                
                % Set the sampled time vector to be the existing time vector (perhaps empty, perhaps a complete time vector).
                ts_sample = ts_reference;
                
                % Set the sampled applied current vector to be the existing time vector (perhaps empty, perhaps a complete time vector).
                Ias_sample = Ias_reference;
                
            end
            
        end
        
        
        %% Enable & Disable Functions.
        
        % Implement a function to toogle whether this applied current is enabled.
        function [ enabled_flag, self ] = toggle_enabled( self, enabled_flag, set_flag )
            
            % Set the default input arguments.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end            % [T/F] Set Flag (Determines whether to update the neuron object.)
            if narign < 2, enabled_flag = self.enabled_flag; end          	% [T/F] Enabled Flag.
            
            % Toggle whether the applied current is enabled.
            enabled_flag = ~enabled_flag;                                 	% [T/F] Applied Current Enabled Flag.
            
            % Determine whether to update the applied current object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        % Implement a function to enable this applied current.
        function [ enabled_flag, self ] = enable( self, set_flag )
            
            % Set the default input arguments.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end            % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Enable this applied current.
            enabled_flag = true;                                         	% [T/F] Applied Current Enabled Flag.
            
            % Determine whether to update the applied current object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        % Implement a function to disable this applied current.
        function [ enabled_flag, self ] = disable( self, set_flag )
            
            % Set the default input arguments.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end            % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Disable this applied current.
            enabled_flag = false;                                        	% [T/F] Applied Current Enabled Flag.
            
            % Determine wehther to update the applied current object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        %% Print Functions.
        
        % Implement a function to print the properties associated with this applied current.
        function print( self, verbose_flag )
           
            % Set the default input arguments.
            if nargin < 2, verbose_flag = self.verbose_flag_DEFAULT; end
            
            % Print out the information for this neuron.
            self.applied_current_utilities.print( self.ID, self.name, self.to_neuron_ID, self.ts, self.Ias, self.num_timesteps, self.dt, self.tf, self.enabled_flag, verbose_flag );
            
        end
        
        
        %% Save & Load Functions.
        
        % Implement a function to save applied current data as a matlab object.
        function save( self, directory, file_name, applied_current )
            
            % Set the default input arguments.
            if nargin < 4, applied_current = self; end                      % [class] Array of Applied Current Class Objects.
            if nargin < 3, file_name = 'Applied_Current.mat'; end           % [str]  File Name.
            if nargin < 2, directory = '.'; end                             % [str] Save Directory.
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Save the neuron data.
            save( full_path, applied_current )
            
        end
        
        
        % Implement a function to load applied current data as a matlab object.
        function applied_current = load( ~, directory, file_name )
            
            % Set the default input arguments.
            if nargin < 3, file_name = 'Applied_Current.mat'; end           % [str] File Name.
            if nargin < 2, directory = '.'; end                             % [str] Load Directory.
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Load the data.
            data = load( full_path );
            
            % Retrieve the desired variable from the loaded data structure.
            applied_current = data.applied_current;
            
        end
        
        
    end
end

