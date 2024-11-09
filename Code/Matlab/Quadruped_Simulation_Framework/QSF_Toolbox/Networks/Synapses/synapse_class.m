classdef synapse_class
    
    % This class contains properties and methods related to synapses.
    
    %% SYNAPSE PROPERTIES
    
    % Define the class properties.
    properties
        
        ID                                                    	% [#] Synapse ID.
        name                                                  	% [-] Synapse Name.
        
        dEs                                                 	% [V] Synaptic Reversal Potential.
        gs                                                      % [S] Maximum Synaptic Conductance.
        Gs                                                      % [S] Synaptic Conductance.
        
        from_neuron_ID                                       	% [#] From Neuron ID.
        to_neuron_ID                                           	% [#] To Neuron ID.
        
        delta                                                 	% [V] CPG Equilibrium Offset.
        
        enabled_flag                                            % [T/F] Synapse Enabled Flag.
        
        synapse_utilities                                       % [-] Synapse Utilities Class.
        
    end
    
    
    % Define private, constant class properties.
    properties ( Access = private, Constant = true )
        
        % ---------- Neuron Properties ----------
        
        % Define the neuron parameters.
        R_DEFAULT = 20e-3;                                    	% [V] Activation Domain.
        Gm_DEFAULT = 1e-6;                                     	% [S] Membrane Conductance.
       
        
        % ---------- Synapse Properties ----------

        % Define the maximum synaptic conductance.
        gs_DEFAULT = 1e-6;                                     	% [S] Maximum Synaptic Conductance.
        Gs_DEFAULT = 0.0;                                     	% [S] Synaptic Conductance.
        
        % Define the synaptic reversal potential parameters.
        dEs_maximum_DEFAULT = 194e-3;                        	% [V] Maximum Synaptic Reversal Potential.
        dEs_minimum_DEFAULT = -40e-3;                         	% [V] Minimum Synaptic Reversal Potential.
        dEs_small_negative_DEFAULT = -1e-3;                     % [V] Small Negative Synaptic Reversal Potential.
        
        % Define the synapse identification information.
        to_neuron_ID_DEFAULT = 0;                            	% [#] To Neuron ID.
        from_neuron_ID_DEFAULT = 0;                           	% [#] From Neuron ID.
        ID_DEFAULT = 0;                                       	% [#] Synapse ID.
        name_DEFAULT = '';                                      % [str] Synapse Name.
        
        
        % ---------- Transmission Subnetwork Properties ----------
        
        % Define the applied current magnitude properties.
        Ia2_absolute_transmission_DEFAULT = 0.0;              	% [A] Absolute Transmission Applied Current 2.
        Ia2_relative_transmission_DEFAULT = 0.0;              	% [A] Relative Transmission Applied Current 2.
    
        % Define the transmission subnetwork gains.
        c_absolute_transmission_DEFAULT = 1.0;                  % [-] Absolute Transmission Subnetwork Gain.
        c_relative_tranmission_DEFAULT = 1.0;                   % [-] Relative Transmission Subnetwork Gain.
        
        
        % ---------- Addition Subnetwork Properties ----------

        % Define the applied current magnitude properties.
        Ia_absolute_addition_DEFAULT = 0.0;                   	% [A] Absolute Addition Applied Current.
        Ia_relative_addition_DEFAULT = 0.0;                  	% [A] Relative Addition Applied Current.
        
        % Define the subnetwork gain properties.
        c_absolute_addition_DEFAULT = 1.0;                  	% [-] Absolute Addition Subnetwork Gain.
        c_relative_addition_DEFAULT = 1.0;                   	% [-] Relative Addition Subnetwork Gain.
       
        
        % ---------- Subtraction Subnetwork Properties ----------

        % Define the applied current magnitude properties.
        Ia_absolute_subtraction_DEFAULT = 0.0;              	% [A] Absolute Subtraction Applied Current.
        Ia_relative_subtraction_DEFAULT = 0.0;                	% [A] Relative Subtraction Applied Current.
        
        % Define the subnetwork gain properties.
        c_absolute_subtraction_DEFAULT = 1.0;                	% [-] Absolute Subtraction Subnetwork Gain.
        c_relative_subtraction_DEFAULT = 1.0;                  	% [-] Relative Subtraction Subnetwork Gain.
        
        
        % ---------- Inversion Subnetwork Properties ----------

        % Define the applied current magnitude properties.
        Ia2_absolute_inversion_DEFAULT = 20e-9;                	% [A] Absolute Inversion Applied Current 2.
        Ia2_relative_inversion_DEFAULT = 20e-9;                	% [A] Relative Inversion Applied Current 2.
        
        % Define the subnetwork gain properties.
        c_absolute_inversion_DEFAULT = 1.0;                   	% [-] Absolute Inversion Subnetwork Gain.
        c_relative_inversion_DEFAULT = 1.0;                   	% [-] Relative Inversion Subnetwork Gain.
        
        
        % ---------- Reduced Inversion Subnetwork Properties ----------

        % Define the reduced inversion applied current magnitudes.
        Ia2_reduced_absolute_inversion_DEFAULT = 20e-9;         % [A] Reduced Absolute Inversion Applied Current 2.
        Ia2_reduced_relative_inversion_DEFAULT = 20e-9;         % [A] Reduced Relative Inversion Applied Current 2.
        
        % Define the reduced inversion subnetwork gains.
        c_reduced_absolute_inversion_DEFAULT = 1.0;             % [-] Reduced Absolute Inversion Subnetwork Gain.
        c_reduced_relative_inversion_DEFAULT = 1.0;             % [-] Reduced Relative Inversion Subnetwork Gain.
        
        
        % ---------- Division Subnetwork Properties ----------

        % Define the applied current magnitude properties.
        Ia_absolute_division_DEFAULT = 0.0;                    	% [A] Absolute Division Applied Current.
        Ia_relative_division_DEFAULT = 0.0;                   	% [A] Relative Division Applied Current.
        
        % Define the subnetwork gain properties.
        c_absolute_division_DEFAULT = 1.0;                    	% [-] Absolute Division Subnetwork Gain.
        c_relative_division_DEFAULT = 1.0;                    	% [-] Relative Division Subnetwork Gain.
        
        
        % ---------- Reduced Division Subnetwork Properties ----------

        % Define the reduced division applied current magnitudes.
        Ia3_reduced_absolute_division_DEFAULT = 0.0;          	% [A] Reduced Absolute Division Applied Current.
        Ia3_reduced_relative_division_DEFAULT = 0.0;          	% [A] Reduced Relative Division Applied Current.
        
        % Define the reduced division subnetwork gains.
        c_reduced_absolute_division_DEFAULT = 1.0;              % [-] Reduced Absolute Division Subnetwork Gain.
        c_reduced_relative_division_DEFAULT = 1.0;              % [-] Reduced Relative Division Subnetwork Gain.
        
        
        % ---------- Division After Inversion Subnetwork Properties ----------

        % Define the division after inversion applied current magnitudes.
        Ia3_absolute_dai_DEFAULT = 0.0;                         % [A] Absolute Division After Inversion Applied Current.
        Ia3_relative_dai_DEFAULT = 0.0;                         % [A] Relative Division After Inversion Applied Current.
        
        % Define the division after inversion subnetwork gains.
        c_absolute_dai_DEFAULT = 1.0;                       	% [-] Absolute Division After Inversion Subnetwork Gain.
        c_relative_dai_DEFAULT = 1.0;                           % [-] Relative Division After Inversion Subnetwork Gain.
        
        
        % ---------- Reduced Division After Inversion Subnetwork Properties ----------

        % Define the reduced division after inversion applied current magnitudes.
        Ia3_reduced_absolute_dai_DEFAULT = 0.0;                	% [A] Reduced Absolute Division After Inversion Applied Current.
        Ia3_reduced_relative_dai_DEFAULT = 0.0;                	% [A] Reduced Relative Division After Inversion Applied Current.
        
        % Define the reduced division after inversion subnetwork gains.
        c_reduced_absolute_dai_DEFAULT = 1.0;                  	% [-] Reduced Absolute Division After Inversion Subnetwork Gain.
        c_reduced_relative_dai_DEFAULT = 1.0;                	% [-] Reduced Relative Division After Inversion Subnetwork Gain.

        
        % ---------- Central Pattern Generator Subnetwork Properties ----------

        % Define applied current properties.
        Id_max_DEFAULT = 1.25e-9;                             	% [A] Maximum Drive Current.
        
        % Define the CPG offset parameters.
        delta_oscillatory_DEFAULT = 0.01e-3;                  	% [V] CPG Oscillatory Delta.
        delta_bistable_DEFAULT = -10e-3;                      	% [V] CPG Bistable Delta.
        delta_noncpg_DEFAULT = 0.0;                            	% [V] CPG Delta.
                
        
        % ---------- Synapse Design Properties ----------

        % Define the default encoding scheme.
        encoding_scheme_DEFAULT = 'Absolute';               	% [-] Encoding Scheme.
        
        % Set the default flag values.
        enabled_flag_DEFAULT = true;                            % [T/F] Enabled Flag (Determines whether this synapse is used when simulating.)
        set_flag_DEFAULT = true;                                % [T/F] Set Flag (Determines whther to update the synapse object.)
        
    end
    
    
    %% SYNAPSE METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = synapse_class( ID, name, dEs, gs, from_neuron_ID, to_neuron_ID, delta, enabled_flag, synapse_utilities )
            
            % Set the default synapse properties.
            if nargin < 9, synapse_utilities = synapse_utilities_class(  ); end             % [class] Synapse Utiliities.
            if nargin < 8, enabled_flag = self.enabled_flag_DEFAULT; end                   	% [T/F] Synapse Enabled Flag.
            if nargin < 7, delta = self.delta_noncpg_DEFAULT; end                           % [V] CPG Equilibrium Offset.
            if nargin < 6, to_neuron_ID = self.to_neuron_ID_DEFAULT; end                  	% [#] To Neuron ID.
            if nargin < 5, from_neuron_ID = self.from_neuron_ID_DEFAULT; end             	% [S] Synaptic Conductance.
            if nargin < 4, gs = self.gs_DEFAULT; end                                        % [S] Maximum Synaptic Conductance.
            if nargin < 3, dEs = self.dEs_minimum_DEFAULT; end                              % [V] Synaptic Reversal Potential.
            if nargin < 2, name = ''; end                                                 	% [-] Synapse Name.
            if nargin < 1, ID = self.ID_DEFAULT; end                                      	% [#] Synapse ID.
            
            % Store an instance of the synapse utilities class.
            self.synapse_utilities = synapse_utilities;
            
            % Store whether this synapse is active.
            self.enabled_flag = enabled_flag;
            
            % Store the synapse connectivity information.
            self.from_neuron_ID = from_neuron_ID;
            self.to_neuron_ID = to_neuron_ID;
            
            % Store the synapse properties.
            self.delta = delta;
            self.gs = gs;
            self.Gs = self.Gs_DEFAULT;
            self.dEs = dEs;
            
            % Store the synapse identification information.
            self.name = name;
            self.ID = ID;
                        
        end
        
        
        %% Name Functions.
        
        % Implement a function to generate a name for this synapse.
        function [ name, self ] = generate_name( self, ID, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end
            if nargin < 2, ID = self.ID; end
            
            % Generate a name for the synapse.
            name = synapse_utilities.ID2name( ID );
            
            % Determine whether to update the name.
            if set_flag, self.name = name; end
            
        end
        
                
        %% Synaptic Conductance Unpacking Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute transmission synaptic conductance.
        function [ R2, Gm2, dEs21, Ia2 ] = unpack_absolute_transmission_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                R2 = self.R_DEFAULT;                                    % [V] Activation Domain.
                Gm2 = self.Gm_DEFAULT;                                  % [S] Membrane Conductance.
                dEs21 = self.dEs;                                       % [V] Synaptic Reversal Potential.
                Ia2 = self.Ia_absolute_transmission_DEFAULT;            % [A] Applied Current.
            
            elseif length( parameters ) == 4                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R2 = parameters{ 1 };                                   % [V] Activation Domain.
                Gm2 = parameters{ 2 };                                  % [S] Membrane Conductance.
                dEs21 = parameters{ 3 };                                % [V] Synaptic Reversal Potential.
                Ia2 = parameters{ 4 };                                  % [A] Applied Current.
                
            else                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end  
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative transmission synaptic conductance.
        function [ R2, Gm2, dEs21, Ia2 ] = unpack_relative_transmission_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                R2 = self.R_DEFAULT;                                  	% [V] Activation Domain.
                Gm2 = self.Gm_DEFAULT;                                	% [S] Membrane Conductance.
                dEs21 = self.dEs;                                       % [V] Synaptic Reversal Potential.
                Ia2 = self.Ia_absolute_transmission_DEFAULT;            % [A] Applied Current.
            
            elseif length( parameters ) == 4                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R2 = parameters{ 1 };                                   % [V] Activation Domain.
                Gm2 = parameters{ 2 };                                  % [S] Membrane Conductance.
                dEs21 = parameters{ 3 };                                % [V] Synaptic Reversal Potential.
                Ia2 = parameters{ 4 };                                  % [A] Applied Current.
                
            else                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end  
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
            
        % Implement a function to unpack the parameters required to compute the absolute addition synaptic conductance.
        function [ c_k, R_k, Gm_n, dEs_nk, Ia_n ] = unpack_absolute_addition_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                            % If the parameters are empty...
            
                % Set the parameters to default values.
                c_k = self.c_absolute_addition_DEFAULT;        	% [-] Absolute Addition Subnetwork Gain.
                R_k = self.R_DEFAULT;                          	% [V] Activation Domain.
                Gm_n = self.Gm_DEFAULT;                        	% [S] Membrane Conductance.
                dEs_nk = self.dEs;                              % [V] Synaptic Reversal Potential.
                Ia_n = self.Ia_absolute_addition_DEFAULT;       % [A] Applied Current.
            
            elseif length( parameters ) == 5                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c_k = parameters{ 1 };                       	% [-] Absolute Addition Subnetwork Gain.
                R_k = parameters{ 2 };                          % [V] Activation Domain.
                Gm_n = parameters{ 3 };                         % [S] Membrane Conductance.
                dEs_nk = parameters{ 4 };                       % [V] Synaptic Reversal Potential.
                Ia_n = parameters{ 5 };                         % [A] Applied Current.
                
            else                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative addition synaptic conductance.
        function [ c_k, R_n, Gm_n, dEs_nk, Ia_n ] = unpack_relative_addition_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                            % If the parameters are empty...
            
                % Set the parameters to default values.
                c_k = self.c_relative_addition_DEFAULT;       	% [-] Absolute Addition Subnetwork Gain.
                R_n = self.R_DEFAULT;                          	% [V] Activation Domain.
                Gm_n = self.Gm_DEFAULT;                        	% [S] Membrane Conductance.
                dEs_nk = self.dEs;                              % [V] Synaptic Reversal Potential.
                Ia_n = self.Ia_relative_addition_DEFAULT;       % [A] Applied Current.
            
            elseif length( parameters ) == 5                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c_k = parameters{ 1 };                         	% [-] Absolute Addition Subnetwork Gain.
                R_n = parameters{ 2 };                          % [V] Activation Domain.
                Gm_n = parameters{ 3 };                         % [S] Membrane Conductance.
                dEs_nk = parameters{ 4 };                       % [V] Synaptic Reversal Potential.
                Ia_n = parameters{ 5 };                         % [A] Applied Current.
                
            else                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute subtraction synaptic conductance.
        function [ c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n ] = unpack_absolute_subtraction_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...
            
                % Set the parameters to default values.
                c_k = self.c_absolute_subtraction_DEFAULT;         	% [-] Absolute Subtraction Subnetwork Gain.
                s_k = self.signature_DEFAULT;                     	% [-] Excitation / Inhibition Sign.
                R_k = self.R_DEFAULT;                             	% [V] Activation Domain.
                Gm_n = self.Gm_DEFAULT;                          	% [S] Membrane Conductance.
                dEs_nk = self.dEs;                                  % [V] Synaptic Reversal Potential.
                Ia_n = self.Ia_absolute_subtraction_DEFAULT;        % [A] Applied Current.
                
            elseif length( parameters ) == 6                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c_k = parameters{ 1 };                           	% [-] Absolute Subtraction Subnetwork Gain.
                s_k = parameters{ 2 };                            	% [-] Excitation / Inhibition Sign.
                R_k = parameters{ 3 };                              % [V] Activation Domain.
                Gm_n = parameters{ 4 };                             % [S] Membrane Conductance.
                dEs_nk = parameters{ 5 };                           % [V] Synaptic Reversal Potential.
                Ia_n = parameters{ 6 };                             % [A] Applied Current.
            
            else                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative subtraction synaptic conductance.
        function [ c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n ] = unpack_relative_subtraction_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...
            
                % Set the parameters to default values.
                c_k = self.c_relative_subtraction_DEFAULT;        	% [-] Absolute Subtraction Subnetwork Gain.
                s_k = self.signature_DEFAULT;                    	% [-] Excitation / Inhibition Sign.
                R_k = self.R_DEFAULT;                             	% [V] Activation Domain.
                Gm_n = self.Gm_DEFAULT;                          	% [S] Membrane Conductance.
                dEs_nk = self.dEs;                                  % [V] Synaptic Reversal Potential.
                Ia_n = self.Ia_relative_subtraction_DEFAULT;        % [A] Applied Current.

            elseif length( parameters ) == 6                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c_k = parameters{ 1 };                             	% [-] Absolute Subtraction Subnetwork Gain.
                s_k = parameters{ 2 };                            	% [-] Excitation / Inhibition Sign.
                R_k = parameters{ 3 };                              % [V] Activation Domain.
                Gm_n = parameters{ 4 };                             % [S] Membrane Conductance.
                dEs_nk = parameters{ 5 };                           % [V] Synaptic Reversal Potential.
                Ia_n = parameters{ 6 };                             % [A] Applied Current.
            
            else                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute inversion synaptic conductance.
        function [ delta, Gm2, dEs21, Ia2 ] = unpack_absolute_inversion_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_absolute_inversion_DEFAULT;          % [V] Absolute Inversion Offset.
                Gm2 = self.Gm_DEFAULT;                                  % [S] Membrane Conductance.
                dEs21 = self.dEs;                                       % [V] Synaptic Reversal Potential.
                Ia2 = self.Ia2_absolute_inversion_DEFAULT;              % [A] Applied Current.

            elseif length( parameters ) == 4                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                % [V] Absolute Inversion Offset.
                Gm2 = parameters{ 2 };                                  % [S] Membrane Conductance.
                dEs21 = parameters{ 3 };                                % [V] Synaptic Reversal Potential.
                Ia2 = parameters{ 4 };                                 	% [A] Applied Current.
            
            else                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
            
        % Implement a function to unpack the parameters required to compute the relative inversion synaptic conductance.
        function [ delta, Gm2, dEs21, Ia2 ] = unpack_relative_inversion_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_relative_inversion_DEFAULT;          % [V] Relative Inversion Offset.
                Gm2 = self.Gm_DEFAULT;                                  % [S] Membrane Conductance.
                dEs21 = self.dEs;                                       % [V] Synaptic Reversal Potential.
                Ia2 = self.Ia2_relative_inversion_DEFAULT;              % [A] Applied Current.

            elseif length( parameters ) == 4                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                % [V] Relative Inversion Offset.
                Gm2 = parameters{ 2 };                                  % [S] Membrane Conductance.
                dEs21 = parameters{ 3 };                                % [V] Synaptic Reversal Potential.
                Ia2 = parameters{ 4 };                                  % [A] Applied Current.
            
            else                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the reduced absolute inversion synaptic conductance.
        function [ delta, Gm2, dEs21, Ia2 ] = unpack_reduced_absolute_inversion_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_reduced_absolute_inversion_DEFAULT;          % [V] Reduced Absolute Inversion Offset.
                Gm2 = self.Gm_DEFAULT;                                          % [S] Membrane Conductance.
                dEs21 = self.dEs;                                               % [V] Synaptic Reversal Potential.
                Ia2 = self.Ia2_reduced_absolute_inversion_DEFAULT;              % [A] Applied Current.

            elseif length( parameters ) == 4                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                        % [V] Reduced Absolute Inversion Offset.
                Gm2 = parameters{ 2 };                                          % [S] Membrane Conductance.
                dEs21 = parameters{ 3 };                                        % [V] Synaptic Reversal Potential.
                Ia2 = parameters{ 4 };                                          % [A] Applied Current.
            
            else                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
            
        % Implement a function to unpack the parameters required to compute the reduced relative inversion synaptic conductance.
        function [ delta, Gm2, dEs21, Ia2 ] = unpack_reduced_relative_inversion_gs_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_reduced_relative_inversion_DEFAULT;          % [V] Reduced Relative Inversion Offset.
                Gm2 = self.Gm_DEFAULT;                                          % [S] Membrane Conductance.
                dEs21 = self.dEs;                                               % [V] Synaptic Reversal Potential.
                Ia2 = self.Ia2_reduced_relative_inversion_DEFAULT;              % [A] Applied Current.

            elseif length( parameters ) == 4                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                        % [V] Reduced Relative Inversion Offset.
                Gm2 = parameters{ 2 };                                          % [S] Membrane Conductance.
                dEs21 = parameters{ 3 };                                        % [V] Synaptic Reversal Potential.
                Ia2 = parameters{ 4 };                                          % [A] Applied Current.
            
            else                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Division Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to unpack the parameters required to compute the absolute division synaptic conductance.
        function [ R3, Gm3, dEs31, Ia3 ] = unpack_absolute_division_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                              % [S] Membrane Conductance.
                dEs31 = self.dEs;                                   % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_absolute_division_DEFAULT;            % [A] Applied Current.
                
            elseif length( parameters ) == 4                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = parameters{ 1 };                               % [V] Activation Domain.
                Gm3 = parameters{ 2 };                              % [S] Membrane Conductance.
                dEs31 = parameters{ 3 };                          	% [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                              % [A] Applied Current.
            
            else                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative division synaptic conductance.
        function [ R3, Gm3, dEs31, Ia3 ] = unpack_relative_division_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                              % [S] Membrane Conductance.
                dEs31 = self.dEs;                                   % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_relative_division_DEFAULT;            % [A] Applied Current.

            elseif length( parameters ) == 4                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = parameters{ 1 };                               % [V] Activation Domain.
                Gm3 = parameters{ 2 };                             	% [S] Membrane Conductance.
                dEs31 = parameters{ 3 };                            % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                              % [A] Applied Current.
            
            else                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Division Subnetwork Functions (Synapse 32) ----------

        % Implement a function to unpack the parameters required to compute the absolute division synaptic conductance.
        function [ delta, Gm3, gs31, dEs31, dEs32, Ia3 ] = unpack_absolute_division_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_absolute_division_DEFAULT;           % [V] Absolute Division Offset.
                Gm3 = self.Gm_DEFAULT;                                  % [S] Membrane Conductance.
                gs31 = self.gs_DEFAULT;                                 % [S] Synaptic Conductance.
                dEs31 = self.dEs_DEFAULT;                               % [V] Synaptic Reversal Potential.
                dEs32 = self.dEs;                                       % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_absolute_division_DEFAULT;                % [A] Applied Current.
                
            elseif length( parameters ) == 6                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                % [V] Absolute Division Offset.
                Gm3 = parameters{ 2 };                                  % [S] Membrane Conductance.
                gs31 = parameters{ 3 };                                 % [S] Synaptic Conductance.
                dEs31 = parameters{ 4 };                                % [V] Synaptic Reversal Potential.
                dEs32 = parameters{ 5 };                                % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 6 };                                  % [A] Applied Current.
                
            else                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative division synaptic conductance.
        function [ delta, Gm3, gs31, dEs31, dEs32, Ia3 ] = unpack_relative_division_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_absolute_division_DEFAULT;           % [V] Absolute Division Offset.
                Gm3 = self.Gm_DEFAULT;                                  % [S] Membrane Conductance.
                gs31 = self.gs_DEFAULT;                                 % [S] Synaptic Conductance.
                dEs31 = self.dEs_DEFAULT;                               % [V] Synaptic Reversal Potential.
                dEs32 = self.dEs;                                       % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_absolute_division_DEFAULT;                % [A] Applied Current.
                
            elseif length( parameters ) == 6                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                % [V] Absolute Division Offset.
                Gm3 = parameters{ 2 };                                  % [S] Membrane Conductance.
                gs31 = parameters{ 3 };                                 % [S] Synaptic Conductance.
                dEs31 = parameters{ 4 };                               	% [V] Synaptic Reversal Potential.
                dEs32 = parameters{ 5 };                                % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 6 };                                  % [A] Applied Current.
                
            else                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to unpack the parameters required to compute the reduced absolute division synaptic conductance.
        function [ R3, Gm3, dEs31, Ia3 ] = unpack_reduced_absolute_division_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                        % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                      % [S] Membrane Conductance.
                dEs31 = self.dEs;                                           % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_reduced_absolute_division_DEFAULT;            % [A] Applied Current.
                
            elseif length( parameters ) == 4                                % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = parameters{ 1 };                                       % [V] Activation Domain.
                Gm3 = parameters{ 2 };                                      % [S] Membrane Conductance.
                dEs31 = parameters{ 3 };                                    % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                                      % [A] Applied Current.
            
            else                                                            % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the reduced relative division synaptic conductance.
        function [ R3, Gm3, dEs31, Ia3 ] = unpack_reduced_relative_division_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...
            
                % Set the parameters to default values.
                R3 = self.R_DEFAULT;                                        % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                      % [S] Membrane Conductance.
                dEs31 = self.dEs;                                           % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_reduced_relative_division_DEFAULT;            % [A] Applied Current.

            elseif length( parameters ) == 4                                % If there are a specific number of parameters...
                
                % Unpack the parameters.
                R3 = parameters{ 1 };                                       % [V] Activation Domain.
                Gm3 = parameters{ 2 };                                      % [S] Membrane Conductance.
                dEs31 = parameters{ 3 };                                    % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                                      % [A] Applied Current.
            
            else                                                            % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Synapse 32) ----------

        % Implement a function to unpack the parameters required to compute the reduced absolute division synaptic conductance.
        function [ delta, Gm3, gs31, dEs31, dEs32, Ia3 ] = unpack_reduced_absolute_division_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_reduced_absolute_division_DEFAULT;           % [V] Reduced Absolute Division Offset.
                Gm3 = self.Gm_DEFAULT;                                          % [S] Membrane Conductance.
                gs31 = self.gs_DEFAULT;                                         % [S] Synaptic Conductance.
                dEs31 = self.dEs_DEFAULT;                                       % [V] Synaptic Reversal Potential.
                dEs32 = self.dEs;                                               % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_reduced_absolute_division_DEFAULT;                % [A] Applied Current.
                
            elseif length( parameters ) == 6                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                        % [V] Reduced Absolute Division Offset.
                Gm3 = parameters{ 2 };                                          % [S] Membrane Conductance.
                gs31 = parameters{ 3 };                                         % [S] Synaptic Conductance.
                dEs31 = parameters{ 4 };                                        % [V] Synaptic Reversal Potential.
                dEs32 = parameters{ 5 };                                        % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 6 };                                          % [A] Applied Current.
                
            else                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the reduced relative division synaptic conductance.
        function [ delta, Gm3, gs31, dEs31, dEs32, Ia3 ] = unpack_reduced_relative_division_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta = self.delta_reduced_relative_division_DEFAULT;           % [V] Reduced Relative Division Offset.
                Gm3 = self.Gm_DEFAULT;                                          % [S] Membrane Conductance.
                gs31 = self.gs_DEFAULT;                                         % [S] Synaptic Conductance.
                dEs31 = self.dEs_DEFAULT;                                       % [V] Synaptic Reversal Potential.
                dEs32 = self.dEs;                                               % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_DEFAULT;                                          % [A] Applied Current.
                
            elseif length( parameters ) == 6                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta = parameters{ 1 };                                        % [V] Reduced Relative Division Offset.
                Gm3 = parameters{ 2 };                                          % [S] Membrane Conductance.
                gs31 = parameters{ 3 };                                         % [S] Synaptic Conductance.
                dEs31 = parameters{ 4 };                                        % [V] Synaptic Reversal Potential.
                dEs32 = parameters{ 5 };                                        % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 6 };                                          % [A] Applied Current.
                
            else                                                              	% Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to unpack the parameters required to compute the absolute division after inversion synaptic conductance.
        function [ c1, c3, delta1, delta2, R1, R2 ] = unpack_absolute_dai_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                c1 = self.c1_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 1.
                c3 = self.c3_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 3.
                delta1 = self.delta_absolute_inversion_DEFAULT;                         % [V] Absolute Inversion Offset.
                delta2 = self.delta_absolute_dai_DEFAULT;                               % [V] Absolute Division After Inversion Offset.
                R1 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                R2 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                
            elseif length( parameters ) == 6                                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c1 = parameters{ 1 };                                                   % [-] Absolute Division After Inversion Gain 1.
                c3 = parameters{ 2 };                                                   % [-] Absolute Division After Inversion Gain 3.
                delta1 = parameters{ 3 };                                               % [V] Absolute Inversion Offset.
                delta2 = parameters{ 4 };                                               % [V] Absolute Division After Inversion Offset.
                R1 = parameters{ 5 };                                                   % [V] Activation Domain.
                R2 = parameters{ 6 };                                                   % [V] Activation Domain.
            
            else                                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative division after inversion synaptic conductance.
        function [ c1, c3, delta1, delta2, R2, dEs31 ] = unpack_relative_dai_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                % If the parameters are empty...
            
                % Set the parameters to default values.
                c1 = c1_relative_dai_DEFAULT;                                       % [-] Relative Division After Inversion Gain 1.
                c3 = c3_relative_dai_DEFAULT;                                     	% [-] Relative Division After Inversion Gain 3.
                delta1 = delta_relative_inversion_DEFAULT;                          % [V] Relative Inversion Offset.
                delta2 = delta_relative_dai_DEFAULT;                                % [V] Relative Division After Inversion Offset.
                R2 = self.R_DEFAULT;                                                % [V] Activation Domain.
                dEs31 = self.dEs;                                                   % [V] Synaptic Reversal Potential.

            elseif length( parameters ) == 6                                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c1 = parameters{ 1 };                                               % [-] Relative Division After Inversion Gain 1.
                c3 = parameters{ 2 };                                               % [-] Relative Division After Inversion Gain 3.
                delta1 = parameters{ 3 };                                           % [V] Relative Inversion Offset.
                delta2 = parameters{ 4 };                                           % [V] Relative Division After Inversion Offset.
                R2 = parameters{ 5 };                                               % [V] Activation Domain.
                dEs31 = parameters{ 6 };                                            % [V] Synaptic Reversal Potential.
            
            else                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions (Synapse 32) ----------

        % Implement a function to unpack the parameters required to compute the absolute division after inversion synaptic conductance.
        function [ c1, c3, delta2, R1, R2, dEs31 ] = unpack_absolute_dai_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                c1 = self.c1_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 1.
                c3 = self.c3_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 3.
                delta2 = self.delta_absolute_dai_DEFAULT;                            	% [V] Absolute Division After Inversion Offset.
                R1 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                R2 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                dEs31 = self.dEs_DEFAULT;                                               % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c1 = parameters{ 1 };                                                   % [-] Absolute Division After Inversion Gain 1.
                c3 = parameters{ 2 };                                                   % [-] Absolute Division After Inversion Gain 3.
                delta2 = parameters{ 3 };                                             	% [V] Absolute Division After Inversion Offset.
                R1 = parameters{ 4 };                                                   % [V] Activation Domain.
                R2 = parameters{ 5 };                                               	% [V] Activation Domain.
                dEs31 = parameters{ 6 };                                               	% [V] Synaptic Reversal Potential.
                
            else                                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative division after inversion synaptic conductance.
        function [ c1, c3, delta1, delta2, R2, dEs31 ] = unpack_relative_dai_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                c1 = self.c1_relative_dai_DEFAULT;                                      % [-] Relative Division After Inversion Gain 1.
                c3 = self.c3_relative_dai_DEFAULT;                                     	% [-] Relative Division After Inversion Gain 3.
                delta1 = self.delta_relative_inversion_DEFAULT;                         % [V] Relative Inversion Offset.
                delta2 = self.delta_relative_dai_DEFAULT;                               % [V] Relative Division After Inversion Offset.
                R2 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                dEs31 = self.dEs_DEFAULT;                                               % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c1 = parameters{ 1 };                                                   % [-] Relative Division After Inversion Gain 1.
                c3 = parameters{ 2 };                                                   % [-] Relative Division After Inversion Gain 3.
                delta1 = parameters{ 3 };                                               % [V] Relative Inversion Offset.
                delta2 = parameters{ 4 };                                               % [V] Relative Division After Inversion Offset.
                R2 = parameters{ 5 };                                                  	% [V] Activation Domain.
                dEs31 = parameters{ 6 };                                                % [V] Synaptic Reversal Potential.
                
            else                                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to unpack the parameters required to compute the reduced absolute division after inversion synaptic conductance.
        function [ delta1, delta2, R2, R3, Gm3, dEs31 ] = unpack_reduced_absolute_dai_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_absolute_inversion_DEFAULT;                         % [V] Reduced Absolute Inversion Offset.
                delta2 = self.delta_reduced_absolute_dai_DEFAULT;                               % [V] Reduced Absolute Division After Inversion Offset.
                R2 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                                          % [S] Membrane Conductance.
                dEs31 = self.dEs;                                                               % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                                	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Reduced Absolute Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Reduced Absolute Division After Inversion Offset.
                R2 = parameters{ 3 };                                                           % [V] Activation Domain.
                R3 = parameters{ 4 };                                                         	% [V] Activation Domain.
                Gm3 = parameters{ 5 };                                                          % [S] Membrane Conductance.
                dEs31 = parameters{ 6 };                                                        % [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative division after inversion synaptic conductance.
        function [ delta1, delta2, R2, R3, dEs31 ] = unpack_reduced_relative_dai_gs31_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_relative_inversion_DEFAULT;                         % [V] Reduced Relative Inversion Offset.
                delta2 = self.delta_reduced_relative_dai_DEFAULT;                               % [V] Reduced Rleative Division After Inversion Offset.
                R2 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                dEs31 = self.dEs;                                                               % [V] Synaptic Reversal Potential.

            elseif length( parameters ) == 5                                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Reduced Relative Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Reduced Rleative Division After Inversion Offset.
                R2 = parameters{ 3 };                                                           % [V] Activation Domain.
                R3 = parameters{ 4 };                                                           % [V] Activation Domain.
                dEs31 = parameters{ 5 };                                                        % [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Synapse 32) ----------
        
        % Implement a function to unpack the parameters required to compute the reduced absolute division after inversion synapse 32.
        function [ delta1, delta2, R2, R3, Gm3, dEs31 ] = unpack_reduced_absolute_dai_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_absolute_inversion_DEFAULT;                         % [V] Reduced Absolute Inversion Offset.
                delta2 = self.delta_reduced_absolute_dai_DEFAULT;                               % [V] Reduced Absolute Division After Inversion Offset.
                R2 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                                          % [S] Membrane Conductance.
                dEs31 = self.dEs_DEFAULT;                                                       % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                                  	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Reduced Absolute Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Reduced Absolute Division After Inversion Offset.
                R2 = parameters{ 3 };                                                           % [V] Activation Domain.
                R3 = parameters{ 4 };                                                           % [V] Activation Domain.
                Gm3 = parameters{ 5 };                                                          % [S] Membrane Conductance.
                dEs31 = parameters{ 6 };                                                        % [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative division after inversion synapse 32.
        function [ delta1, delta2, R2, R3, Gm3, dEs31 ] = unpack_reduced_relative_dai_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_relative_inversion_DEFAULT;                         % [V] Reduced Relative Inversion Offset.
                delta2 = self.delta_reduced_relative_dai_DEFAULT;                               % [V] Reduced Relative Division After Inversion Offset.
                R2 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                Gm3 = self.Gm_DEFAULT;                                                          % [S] Membrane Conductance.
                dEs31 = self.dEs_DEFAULT;                                                       % [V] Synaptic Reversal Potential.

            elseif length( parameters ) == 6                                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Reduced Relative Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Reduced Relative Division After Inversion Offset.
                R2 = parameters{ 3 };                                                           % [V] Activation Domain.
                R3 = parameters{ 4 };                                                           % [V] Activation Domain.
                Gm3 = parameters{ 5 };                                                          % [S] Membrane Conductance.
                dEs31 = parameters{ 6 };                                                     	% [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 41) ----------
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 41 of an absolute multiplication subnetwork.
        function [ c4, c6, delta1, delta2, R1, R3 ] = unpack_absolute_multiplication_gs41_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                c4 = self.c1_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 1.
                c6 = self.c3_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 3.
                delta1 = self.delta_absolute_inversion_DEFAULT;                         % [V] Absolute Inversion Offset.
                delta2 = self.delta_absolute_dai_DEFAULT;                               % [V] Absolute Division Offset.
                R1 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                R3 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                
            elseif length( parameters ) == 6                                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c4 = parameters{ 1 };                                                   % [-] Absolute Division After Inversion Gain 1.
                c6 = parameters{ 2 };                                                   % [-] Absolute Division After Inversion Gain 3.
                delta1 = parameters{ 3 };                                               % [V] Absolute Inversion Offset.
                delta2 = parameters{ 4 };                                               % [V] Absolute Division Offset.
                R1 = parameters{ 5 };                                                   % [V] Activation Domain.
                R3 = parameters{ 6 };                                                   % [V] Activation Domain.
            
            else                                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 41 of a relative multiplication subnetwork.
        function [ c4, c6, delta1, delta2, R3, dEs41 ] = unpack_relative_multiplication_gs41_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                % If the parameters are empty...
            
                % Set the parameters to default values.
                c4 = self.c1_relative_dai_DEFAULT;                                  % [-] Absolute Division After Inversion Gain 1.
                c6 = self.c3_relative_dai_DEFAULT;                                  % [-] Absolute Division After Inversion Gain 3.
                delta1 = self.delta_relative_inversion_DEFAULT;                     % [V] Absolute Inversion Offset.
                delta2 = self.delta_relative_division_DEFAULT;                      % [V] Absolute Division Offset.
                R3 = self.R_DEFAULT;                                                % [V] Activation Domain.
                dEs41 = self.dEs;                                                   % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c4 = parameters{ 1 };                                               % [-] Absolute Division After Inversion Gain 1.
                c6 = parameters{ 2 };                                               % [-] Absolute Division After Inversion Gain 3.
                delta1 = parameters{ 3 };                                           % [V] Absolute Inversion Offset.
                delta2 = parameters{ 4 };                                           % [V] Absolute Division Offset.
                R3 = parameters{ 5 };                                               % [V] Activation Domain.
                dEs41 = parameters{ 6 };                                           	% [V] Synaptic Reversal Potential.
            
            else                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 32) ----------

        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 32 of an absolute multiplication subnetwork.
        function [ delta1, Gm3, dEs32, Ia3 ] = unpack_absolute_multiplication_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_absolute_inversion_DEFAULT;             % [V] Absolute Inversion Offset.
                Gm3 = self.Gm_DEFAULT;                                      % [S] Membrane Conductance.
                dEs32 = self.dEs;                                           % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_absolute_inversion_DEFAULT;                   % [A] Applied Current.
                
            elseif length( parameters ) == 4                               	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                   % [V] Absolute Inversion Offset.
                Gm3 = parameters{ 2 };                                      % [S] Membrane Conductance.
                dEs32 = parameters{ 3 };                                    % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                                      % [A] Applied Current.
            
            else                                                            % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 32 of a relative multiplication subnetwork.
        function [ delta1, Gm3, dEs32, Ia3 ] = unpack_relative_multiplication_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_relative_inversion_DEFAULT;             % [V] Relative Inversion Offset.
                Gm3 = self.Gm_DEFAULT;                                      % [S] Membrane Conductance.
                dEs32 = self.dEs;                                           % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_relative_inversion_DEFAULT;                   % [A] Applied Current.
                
            elseif length( parameters ) == 4                                % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                   % [V] Relative Inversion Offset.
                Gm3 = parameters{ 2 };                                   	% [S] Membrane Conductance.
                dEs32 = parameters{ 3 };                                    % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                                      % [A] Applied Current.
            
            else                                                            % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 43) ----------

        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 43 of an absolute multiplication subnetwork.
        function [ c4, c6, delta2, R1, R3, dEs41 ] = unpack_absolute_multiplication_gs43_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                c4 = self.c1_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 1.
                c6 = self.c3_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 3.
                delta2 = self.delta_absolute_dai_DEFAULT;                               % [V] Absolute Division Offset.
                R1 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                R3 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                dEs41 = self.dEs_DEFAULT;                                               % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c4 = parameters{ 1 };                                                   % [-] Absolute Division After Inversion Gain 1.
                c6 = parameters{ 2 };                                                   % [-] Absolute Division After Inversion Gain 3.
                delta2 = parameters{ 3 };                                               % [V] Absolute Division Offset.
                R1 = parameters{ 4 };                                                   % [V] Activation Domain.
                R3 = parameters{ 5 };                                                   % [V] Activation Domain.
                dEs41 = parameters{ 6 };                                                % [V] Synaptic Reversal Potential.
            
            else                                                                        % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 43 of a relative multiplication subnetwork.
        function [ c4, c6, delta1, delta2, R3, dEs41 ] = unpack_relative_multiplication_gs43_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                       % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                    % If the parameters are empty...
            
                % Set the parameters to default values.
                c4 = self.c1_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 1.
                c6 = self.c3_absolute_dai_DEFAULT;                                      % [-] Absolute Division After Inversion Gain 3.
                delta1 = self.delta_inversion_division_DEFAULT;                         % [V] Absolute Inversion Offset.
                delta2 = self.delta_absolute_dai_DEFAULT;                               % [V] Absolute Division After Inversion Offset.
                R3 = self.R_DEFAULT;                                                    % [V] Activation Domain.
                dEs41 = self.dEs_DEFAULT;                                               % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                c4 = parameters{ 1 };                                                   % [-] Absolute Division After Inversion Gain 1.
                c6 = parameters{ 2 };                                                  	% [-] Absolute Division After Inversion Gain 3.
                delta1 = parameters{ 3 };                                               % [V] Absolute Inversion Offset.
                delta2 = parameters{ 4 };                                               % [V] Absolute Division After Inversion Offset.
                R3 = parameters{ 5 };                                                   % [V] Activation Domain.
                dEs41 = parameters{ 6 };                                                % [V] Synaptic Reversal Potential.
            
            else                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 41) ----------
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 41 of a reduced absolute multiplication subnetwork.
        function [ delta1, delta2, R3, R4, Gm4, dEs41 ] = unpack_reduced_absolute_multiplication_gs41_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_absolute_inversion_DEFAULT;                         % [V] Absolute Inversion Offset.
                delta2 = self.delta_reduced_absolute_dai_DEFAULT;                               % [V] Absolute Division After Inversion Offset.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R4 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                Gm4 = self.Gm_DEFAULT;                                                          % [S] Membrane Conductance.
                dEs41 = self.dEs;                                                               % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Absolute Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Absolute Division Offset.
                R3 = parameters{ 3 };                                                           % [V] Activation Domain.
                R4 = parameters{ 4 };                                                           % [V] Activation Domain.
                Gm4 = parameters{ 5 };                                                          % [S] Membrane Conductance.
                dEs41 = parameters{ 6 };                                                        % [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 41 of a relative multiplication subnetwork.
        function [ delta1, delta2, R3, R4, dEs41 ] = unpack_reduced_relative_multiplication_gs41_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_relative_inversion_DEFAULT;                         % [V] Relative Inversion Offset.
                delta2 = self.delta_reduced_relative_dai_DEFAULT;                               % [V] Relative Division After Inversion Offset.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R4 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                dEs41 = self.dEs_DEFAULT;                                                       % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 5                                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Relative Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Relative Division After Inversion Offset.
                R3 = parameters{ 3 };                                                           % [V] Activation Domain.
                R4 = parameters{ 4 };                                                           % [V] Activation Domain.
                dEs41 = parameters{ 5 };                                                        % [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 32) ----------

        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 32 of a reduced absolute multiplication subnetwork.
        function [ delta1, Gm3, dEs32, Ia3 ] = unpack_reduced_absolute_multiplication_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_absolute_inversion_DEFAULT;             % [V] Reduced Absolute Inversion Offset.
                Gm3 = self.Gm_DEFAULT;                                              % [S] Membrane Conductance.
                dEs32 = self.dEs;                                                   % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_reduced_absolute_inversion_DEFAULT;                   % [A] Applied Current.
                
            elseif length( parameters ) == 4                                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                           % [V] Reduced Absolute Inversion Offset.
                Gm3 = parameters{ 2 };                                              % [S] Membrane Conductance.
                dEs32 = parameters{ 3 };                                            % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                                              % [A] Applied Current.
            
            else                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 32 of a reduced relative multiplication subnetwork.
        function [ delta1, Gm3, dEs32, Ia3 ] = unpack_reduced_relative_multiplication_gs32_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                   % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_relative_inversion_DEFAULT;             % [V] Reduced Relative Inversion Offset.
                Gm3 = self.Gm_DEFAULT;                                              % [S] Membrane Conductance.
                dEs32 = self.dEs;                                                   % [V] Synaptic Reversal Potential.
                Ia3 = self.Ia_reduced_relative_inversion_DEFAULT;                   % [A] Applied Current.
                
            elseif length( parameters ) == 4                                        % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                           % [V] Reduced Relative Inversion Offset.
                Gm3 = parameters{ 2 };                                              % [S] Membrane Conductance.
                dEs32 = parameters{ 3 };                                            % [V] Synaptic Reversal Potential.
                Ia3 = parameters{ 4 };                                              % [A] Applied Current.
            
            else                                                                    % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
                
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 43) ----------
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 43 of a reduced absolute multiplication subnetwork.
        function [ delta1, delta2, R3, R4, Gm4, dEs41 ] = unpack_reduced_absolute_multiplication_gs43_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_absolute_inversion_DEFAULT;                         % [V] Absolute Inversion Offset.
                delta2 = self.delta_reduced_absolute_dai_DEFAULT;                               % [V] Absolute Division After Inversion Offset.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R4 = self.R_DEFAULT;                                                            % [V} Activation Domain.
                Gm4 = self.Gm_DEFAULT;                                                          % [S] Membrane Conductance.
                dEs41 = self.dEs_DEFAULT;                                                       % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Absolute Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Absolute Division After Inversion Offset.
                R3 = parameters{ 3 };                                                           % [V] Activation Domain.
                R4 = parameters{ 4 };                                                           % [V] Activation Domain.
                Gm4 = parameters{ 5 };                                                          % [S] Membrane Conductance.
                dEs41 = parameters{ 6 };                                                        % [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the synaptic conductance for synapse 43 of a reduced relative multiplication subnetwork.
        function [ delta1, delta2, R3, R4, Gm4, dEs41 ] = unpack_reduced_relative_multiplication_gs43_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                               % [-] Input Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                            % If the parameters are empty...
            
                % Set the parameters to default values.
                delta1 = self.delta_reduced_relative_inversion_DEFAULT;                         % [V] Relative Inversion Offset.
                delta2 = self.delta_reduced_relative_dai_DEFAULT;                             	% [V] Relative Division After Inversion Offset.
                R3 = self.R_DEFAULT;                                                            % [V] Activation Domain.
                R4 = self.R_DEFAULT;                                                            % [V} Activation Domain.
                Gm4 = self.Gm_DEFAULT;                                                          % [S] Membrane Conductance.
                dEs41 = self.dEs_DEFAULT;                                                       % [V] Synaptic Reversal Potential.
                
            elseif length( parameters ) == 6                                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };                                                       % [V] Relative Inversion Offset.
                delta2 = parameters{ 2 };                                                       % [V] Relative Division After Inversion Offset.
                R3 = parameters{ 3 };                                                           % [V] Activation Domain.
                R4 = parameters{ 4 };                                                           % [V] Activation Domain.
                Gm4 = parameters{ 5 };                                                          % [S] Membrane Conductance.
                dEs41 = parameters{ 6 };                                                        % [V] Synaptic Reversal Potential.
            
            else                                                                                % Otherwise...
               
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end 
            
        end
        
        
        %% Synaptic Reversal Potential Compute Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential of a transmission subnetwork.
        function [ dEs21, self ] = compute_transmission_dEs21( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                  % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                             % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a transmission subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue transmission subnetwork.
                dEs21 = synapse_utilities.compute_absolute_transmission_dEs21(  );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative transmission subnetwork.
                dEs21 = synapse_utilities.compute_relative_transmission_dEs21(  );          % [V] Synaptic Reversal Potential.
                
            else                                                                            % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs21; end
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential of addition subnetwork synapses.
        function [ dEs, self ] = compute_addition_dEs( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                     % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end      % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for an addition subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                           	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue addition subnetwork.
                dEs = synapse_utilities.compute_absolute_addition_dEs(  );          % [V] Synaptic Reversal Potential.
                                
            elseif strcmpi( encoding_scheme, 'relative' )                           % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative addition subnetwork.
                dEs = synapse_utilities.compute_relative_addition_dEs(  );          % [V] Synaptic Reversal Potential.
                
            else                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs; end
            
        end

        
        % ---------- Subtraction Subnetwork Functions ----------

        % Implement a function to compute the synaptic reversal potential of subtraction subnetwork excitatory synapses.
        function [ dEs, self ] = compute_subtraction_dEs_excitatory( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                      % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                 % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                  % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a subtraction subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                           % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue subtraction subnetwork.
                dEs = synapse_utilities.compute_absolute_subtraction_dEs_excitatory(  );        % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                       % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative subtraction subnetwork.
                dEs = synapse_utilities.compute_relative_subtraction_dEs_excitatory(  );        % [V] Synaptic Reversal Potential.
                
            else                                                                              	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of subtraction subnetwork inhibitory synapses.
        function [ dEs, self ] = compute_subtraction_dEs_inhibitory( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                          % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                     % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a subtraction subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue subtraction subnetwork.
                dEs = synapse_utilities.compute_absolute_subtraction_dEs_inhibitory(  );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative subtraction subnetwork.
                dEs = synapse_utilities.compute_relative_subtraction_dEs_inhibitory(  );            % [V] Synaptic Reversal Potential.
                
            else                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs; end
            
        end
        
                
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential of inversion subnetwork synapses.
        function [ dEs21, self ] = compute_inversion_dEs21( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                 	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                           	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end            	% [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for an inversion subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                      	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue inversion subnetwork.
                dEs21 = synapse_utilities.compute_absolute_inversion_dEs21(  );             % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                 	% If the encoding scheme is set to relative...

                % Compute the synaptic reversal potential for a relative inversion subnetwork.
                dEs21 = synapse_utilities.compute_relative_inversion_dEs21(  );            	% [V] Synaptic Reversal Potential.
                
            else                                                                           	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs21; end

        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential of reduced inversion subnetwork synapses.
        function [ dEs21, self ] = compute_reduced_inversion_dEs21( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                        	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                    	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                     	% [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for an inversion subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                             	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue inversion subnetwork.
                dEs21 = synapse_utilities.compute_reduced_absolute_inversion_dEs21(  );           	% [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                        	% If the encoding scheme is set to relative...

                % Compute the synaptic reversal potential for a relative inversion subnetwork.
                dEs21 = synapse_utilities.compute_reduced_relative_inversion_dEs21(  );             % [V] Synaptic Reversal Potential.
                
            else                                                                                  	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs21; end

        end
        
        
        % ---------- Division Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential of a division subnetwork.
        function [ dEs31, self ] = compute_division_dEs31( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end            	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                        	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end        	% [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs31 = synapse_utilities.compute_absolute_division_dEs31(  );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                             	% If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs31 = synapse_utilities.compute_relative_division_dEs31(  );       	% [V] Synaptic Reversal Potential.
                
            else                                                                      	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs31; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a division subnetwork.
        function [ dEs32, self ] = compute_division_dEs32( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end            	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                        	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end        	% [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                 	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs32 = synapse_utilities.compute_absolute_division_dEs32(  );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                            	% If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs32 = synapse_utilities.compute_absolute_division_dEs32(  );         	% [V] Synaptic Reversal Potential.
                
            else                                                                      	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs32; end
            
        end
        
    
        % ---------- Reduced Division Subnetwork Functions ----------

        % Implement a function to compute the synaptic reversal potential of a reduced division subnetwork.
        function [ dEs31, self ] = compute_reduced_division_dEs31( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                    	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                               	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                 	% [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                         	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs31 = synapse_utilities.compute_reduced_absolute_division_dEs31(  );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                    	% If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs31 = synapse_utilities.compute_reduced_relative_division_dEs31(  );      	% [V] Synaptic Reversal Potential.
                
            else                                                                               	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs31; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a reduced division subnetwork.
        function [ dEs32, self ] = compute_reduced_division_dEs32( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                      % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                 % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                	% [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                           % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs32 = synapse_utilities.compute_reduced_absolute_division_dEs32(  );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                       % If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs32 = synapse_utilities.compute_reduced_relative_division_dEs32(  );       	% [V] Synaptic Reversal Potential.
                
            else                                                                                % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs32; end
            
        end
                
        
        % ---------- Division After Inversion Subnetwork Functions ----------

        % Implement a function to compute the synaptic reversal potential of a division after inversion subnetwork.
        function [ dEs31, self ] = compute_dai_dEs31( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end            	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                        	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end          % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs31 = synapse_utilities.compute_absolute_dai_dEs31(  );               % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                               % If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs31 = synapse_utilities.compute_relative_dai_dEs31(  );               % [V] Synaptic Reversal Potential.
                
            else                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs31; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a division after inversion subnetwork.
        function [ dEs32, self ] = compute_dai_dEs32( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end             	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                        	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end          % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                  	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs32 = synapse_utilities.compute_absolute_dai_dEs32(  );             	% [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                              	% If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs32 = synapse_utilities.compute_absolute_dai_dEs32(  );              	% [V] Synaptic Reversal Potential.
                
            else                                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs32; end
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------

        % Implement a function to compute the synaptic reversal potential of a reduced division after inversion subnetwork.
        function [ dEs31, self ] = compute_reduced_dai_dEs31( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                 	% [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                           	% [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end             	% [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                      	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs31 = synapse_utilities.compute_reduced_absolute_dai_dEs31(  );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                  	% If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs31 = synapse_utilities.compute_reduced_relative_dai_dEs31(  );        	% [V] Synaptic Reversal Potential.
                
            else                                                                         	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs31; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a reduced division after inversion subnetwork.
        function [ dEs32, self ] = compute_reduced_dai_dEs32( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                  % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                             % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                     	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue division subnetwork.
                dEs32 = synapse_utilities.compute_reduced_absolute_dai_dEs32(  );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
                
                % Compute the synaptic reversal potential for a relative division subnetwork.
                dEs32 = synapse_utilities.compute_reduced_absolute_dai_dEs32(  );        	% [V] Synaptic Reversal Potential.
                
            else                                                                          	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs32; end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------

        % Implement a function to compute the synaptic reversal potential of a multiplication subnetwork.
        function [ dEs41, self ] = compute_multiplication_dEs41( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                      % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                 % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                  % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a multiplication subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                        	% If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue multiplication subnetwork.
                dEs41 = synapse_utilities.compute_absolute_multiplication_dEs41(  );           	% [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                       % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative multiplication subnetwork.
                dEs41 = synapse_utilities.compute_relative_multiplication_dEs41(  );          	% [V] Synaptic Reversal Potential.
                
            else                                                                              	% Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs41; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a multiplication subnetwork.
        function [ dEs32, self ] = compute_multiplication_dEs32( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                      % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                 % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                  % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a multiplication subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                           % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue multiplication subnetwork.
                dEs32 = synapse_utilities.compute_absolute_multiplication_dEs32(  );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                       % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative multiplication subnetwork.
                dEs32 = synapse_utilities.compute_relative_multiplication_dEs32(  );            % [V] Synaptic Reversal Potential.
                
            else                                                                                % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs32; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a multiplication subnetwork.
        function [ dEs43, self ] = compute_multiplication_dEs43( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                      % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                 % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                  % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a multiplication subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                           % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue multiplication subnetwork.
                dEs43 = synapse_utilities.compute_absolute_multiplication_dEs43(  );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                       % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative multiplication subnetwork.
                dEs43 = synapse_utilities.compute_relative_multiplication_dEs43(  );            % [V] Synaptic Reversal Potential.
                
            else                                                                                % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs43; end

        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to compute the synaptic reversal potential of a reduced multiplication subnetwork.
        function [ dEs41, self ] = compute_reduced_multiplication_dEs41( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                              % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                         % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                          % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a multiplication subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                   % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue multiplication subnetwork.
                dEs41 = synapse_utilities.compute_reduced_absolute_multiplication_dEs41(  );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                               % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative multiplication subnetwork.
                dEs41 = synapse_utilities.compute_reduced_relative_multiplication_dEs41(  );          	% [V] Synaptic Reversal Potential.
                
            else                                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs41; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a reduced multiplication subnetwork.
        function [ dEs32, self ] = compute_reduced_multiplication_dEs32( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                              % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                         % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                          % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a multiplication subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                   % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue multiplication subnetwork.
                dEs32 = synapse_utilities.compute_reduced_absolute_multiplication_dEs32(  );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                               % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative multiplication subnetwork.
                dEs32 = synapse_utilities.compute_reduced_relative_multiplication_dEs32(  );          	% [V] Synaptic Reversal Potential.
                
            else                                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs32; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a reduced multiplication subnetwork.
        function [ dEs43, self ] = compute_reduced_multiplication_dEs43( self, encoding_scheme, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 4, synapse_utilities = self.synapse_utilities; end                              % [class] Synapse Utilities.
            if nargin < 3, set_flag = true; end                                                         % [T/F] Set Flag.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end                          % [str] Encoding Scheme.
            
            % Determine how to compute the synaptic reversal potential for a multiplication subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                   % If the encoding scheme is set to absolute...
                
                % Compute the synaptic reversal potential for an absolue multiplication subnetwork.
                dEs43 = synapse_utilities.compute_reduced_absolute_multiplication_dEs43(  );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                               % If the encoding scheme is set to relative...
            
                % Compute the synaptic reversal potential for a relative multiplication subnetwork.
                dEs43 = synapse_utilities.compute_reduced_relative_multiplication_dEs43(  );          	% [V] Synaptic Reversal Potential.
                
            else                                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs43; end

        end
        
        
        % ---------- Derivation Subnetwork Functions ----------

        % Implement a function to compute the synaptic reversal potential of a derivation subnetwork.
        function [ dEs31, self ] = compute_derivation_dEs31( self, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 3, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 2, set_flag = true; end                                     % [T/F] Set Flag.
            
            % Compute the synaptic reversal potential.
            dEs31 = synapse_utilities.compute_derivation_dEs31(  );              	% [V] Synaptic Reversal Potential.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs31; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a derivation subnetwork.
        function [ dEs32, self ] = compute_derivation_dEs32( self, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 3, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 2, set_flag = true; end                                     % [T/F] Set Flag.
            
            % Compute the synaptic reversal potential.
            dEs32 = synapse_utilities.compute_derivation_dEs32(  );                	% [V] Synaptic Reversal Potential.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs32; end
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential of a voltage based integration subnetwork.
        function [ dEs1, self ] = compute_integration_dEs1( self, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 3, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 2, set_flag = true; end                                     % [T/F] Set Flag.
            
            % Compute the synaptic reversal potential.
            dEs1 = synapse_utilities.compute_integration_dEs1(  );                  % [V] Synaptic Reversal Potential.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs1; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a voltage based integration subnetwork.
        function [ dEs2, self ] = compute_integration_dEs2( self, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 3, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 2, set_flag = true; end                                     % [T/F] Set Flag.
            
            % Compute the synaptic reversal potential.
            dEs2 = synapse_utilities.compute_integration_dEs2(  );                  % [V] Synaptic Reversal Potential.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs2; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a voltage based integration subnetwork.
        function [ dEs1, self ] = compute_vbi_dEs1( self, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 3, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 2, set_flag = true; end                                     % [T/F] Set Flag.
            
            % Compute the synaptic reversal potential.
            dEs1 = synapse_utilities.compute_vbi_dEs1(  );                          % [V] Synaptic Reversal Potential.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs1; end
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential of a voltage based integration subnetwork.
        function [ dEs2, self ] = compute_vbi_dEs2( self, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 3, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 2, set_flag = true; end                                     % [T/F] Set Flag.
            
            % Compute the synaptic reversal potential.
            dEs2 = synapse_utilities.compute_vbi_dEs2(  );                          % [V] Synaptic Reversal Potential.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs2; end
            
        end
        
        
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential of a driven multistate cpg subnetwork.
        function [ dEs, self ] = compute_dmcpg_dEs( self, set_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 3, synapse_utilities = self.synapse_utilities; end          % [class] Synapse Utilities.
            if nargin < 2, set_flag = true; end                                     % [T/F] Set Flag.
            
            % Compute the synaptic reversal potential.
            dEs = synapse_utilities.compute_dmcpg_dEs(  );                          % [V] Synaptic Reversal Potential.
            
            % Determine whether to update the synapse object.
            if set_flag, self.dEs = dEs; end
            
        end
        
        
        %% Maximum Synaptic Conductance Compute Functions.
                
        % ---------- Transmission Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 21 of a transmission subnetwork.
        function [ gs21, self ] = compute_transmission_gs21( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                      % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                  % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                 % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                	% [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                               % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an addition subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                           % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute addition subnetwork.                
                [ R2, Gm2, dEs21, Ia2 ] = self.unpack_absolute_transmission_gs_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue addition subnetwork.
                gs21 = synapse_utilities.compute_absolute_transmission_gs21( R2, Gm2, dEs21, Ia2, validation_flag );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                       % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative addition subnetwork.                
                [ R2, Gm2, dEs21, Ia2 ] = self.unpack_relative_transmission_gs_parameters( parameters );
                
                % Compute the synaptic conductance for a relative addition subnetwork.
                gs21 = synapse_utilities.compute_relative_transmission_gs21( R2, Gm2, dEs21, Ia2, validation_flag );            % [V] Synaptic Reversal Potential.
                
            else                                                                                                                % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs21; end

        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance for the synapses of an addition subnetwork.
        function [ gs_nk, self ] = compute_addition_gs( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                          % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                      % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                     % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                      % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                   % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an addition subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute addition subnetwork.
                [ c_k, R_k, Gm_n, dEs_nk, Ia_n ] = self.unpack_absolute_addition_gs_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue addition subnetwork.
                gs_nk = synapse_utilities.compute_absolute_addition_gs( c_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                           % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative addition subnetwork.
                [  c_k, R_n, Gm_n, dEs_nk, Ia_n ] = self.unpack_relative_addition_gs_parameters( parameters );
                
                % Compute the synaptic conductance for a relative addition subnetwork.
                gs_nk = synapse_utilities.compute_relative_addition_gs( c_k, R_n, Gm_n, dEs_nk, Ia_n, validation_flag );            % [V] Synaptic Reversal Potential.
                
            else                                                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs_nk; end

        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for the synapses of a subtraction subnetwork.
        function [ gs_nk, self ] = compute_subtraction_gs( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                  % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                              % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                             % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                              % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                           % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for a subtraction subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                       % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute subtraction subnetwork.
                [ c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n ] = self.unpack_absolute_subtraction_gs_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue subtraction subnetwork.
                gs_nk = synapse_utilities.compute_absolute_subtraction_gs( c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                   % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative subtraction subnetwork.
                [  c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n ] = self.unpack_relative_subtraction_gs_parameters( parameters );
                
                % Compute the synaptic conductance for a relative subtraction subnetwork.
                gs_nk = synapse_utilities.compute_relative_subtraction_gs( c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag );          	% [V] Synaptic Reversal Potential.
                
            else                                                                                                                            % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs_nk; end
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance for synapse 21 of an inversion subnetwork.
        function [ gs21, self ] = compute_inversion_gs21( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                      % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                 	% [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                 % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                  % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                               % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an inversion subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                           % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute inversion subnetwork.
                [ delta1, Gm2, dEs21, Ia2 ] = self.unpack_absolute_inversion_gs_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue inversion subnetwork.
                gs21 = synapse_utilities.compute_absolute_inversion_gs21( delta1, Gm2, dEs21, Ia2, validation_flag );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                       % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative inversion subnetwork.
                [  delta1, Gm2, dEs21, Ia2 ] = self.unpack_relative_inversion_gs_parameters( parameters );
                
                % Compute the synaptic conductance for a relative inversion subnetwork.                
                gs21 = synapse_utilities.compute_relative_inversion_gs21( delta1, Gm2, dEs21, Ia2, validation_flag );          	% [V] Synaptic Reversal Potential.
                 
            else                                                                                                                % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs21; end

        end


        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 21 of a reduced inversion subnetwork.
        function [ gs21, self ] = compute_reduced_inversion_gs21( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                              % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                          % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                         % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                          % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                       % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an inversion subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                   % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute inversion subnetwork.
                [ delta1, Gm2, dEs21, Ia2 ] = self.unpack_reduced_absolute_inversion_gs_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue inversion subnetwork.                
                gs21 = synapse_utilities.compute_reduced_absolute_inversion_gs21( delta1, Gm2, dEs21, Ia2, validation_flag );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                               % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative inversion subnetwork.
                [  delta1, Gm2, dEs21, Ia2 ] = self.unpack_reduced_relative_inversion_gs_parameters( parameters );
                
                % Compute the synaptic conductance for a relative inversion subnetwork.                
                gs21 = synapse_utilities.compute_reduced_relative_inversion_gs21( delta1, Gm2, dEs21, Ia2, validation_flag );         	% [V] Synaptic Reversal Potential.
                 
            else                                                                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs21; end

        end
        
        
        % ---------- Division Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 31 of a division subnetwork.
        function [ gs31, self ] = compute_division_gs31( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                  % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                            	% [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                             % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                              % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                           % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                       % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ R3, Gm3, dEs31, Ia3 ] = self.unpack_absolute_division_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs31 = synapse_utilities.compute_absolute_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                   % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  R3, Gm3, dEs31, Ia3 ] = self.unpack_relative_division_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs31 = synapse_utilities.compute_relative_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );          	% [V] Synaptic Reversal Potential.
                
            else                                                                                                            % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs31; end
            
        end


        % Implement a function to compute the maximum synaptic conductance for synapse 32 of a division subnetwork.
        function [ gs32, self ] = compute_division_gs32( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                  % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                              % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                             % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                              % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                           % [-] Input Parameter Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                       % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta2, Gm3, gs31, dEs31, dEs32, Ia3 ] = self.unpack_absolute_division_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs32 = synapse_utilities.compute_absolute_division_gs32( delta2, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                   % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  delta2, Gm3, gs31, dEs31, dEs32, Ia3 ] = self.unpack_relative_division_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs32 = synapse_utilities.compute_relative_division_gs32( delta2, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );           % [V] Synaptic Reversal Potential.
                
            else                                                                                                                            % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs32; end
           
        end

        
        % ---------- Reduced Division Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 31 of a reduced division subnetwork.
        function [ gs31, self ] = compute_reduced_division_gs31( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                          % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                      % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                     % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                   	% [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                   % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ R3, Gm3, dEs31, Ia3 ] = self.unpack_reduced_absolute_division_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs31 = synapse_utilities.compute_reduced_absolute_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );            % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                           % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  R3, Gm3, dEs31, Ia3 ] = self.unpack_reduced_relative_division_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs31 = synapse_utilities.compute_reduced_relative_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );          	% [V] Synaptic Reversal Potential.
                
            else                                                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs31; end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance for synapse 32 of a reduced division subnetwork.
        function [ gs32, self ] = compute_reduced_division_gs32( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                          % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                                      % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                                     % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                                     	% [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                                   % [-] Input Parameter Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta2, Gm3, gs31, dEs31, dEs32, Ia3 ] = self.unpack_reduced_absolute_division_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs32 = synapse_utilities.compute_reduced_absolute_division_gs32( delta2, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                           % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  delta2, Gm3, gs31, dEs31, dEs32, Ia3 ] = self.unpack_reduced_relative_division_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs32 = synapse_utilities.compute_reduced_relative_division_gs32( delta2, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );         	% [V] Synaptic Reversal Potential.
                
            else                                                                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs32; end
           
        end

        
        % ---------- Division After Inversion Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 31 of a division after inversion subnetwork.
        function [ gs31, self ] = compute_dai_gs31( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                          % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                      % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                     % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                      % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                   % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ c1, c3, delta1, delta2, R1, R2 ] = self.unpack_absolute_dai_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs31 = synapse_utilities.compute_absolute_dai_gs31( c1, c3, delta1, delta2, R1, R2, validation_flag );              % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                           % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  c1, c3, delta1, delta2, R2, dEs31 ] = self.unpack_relative_dai_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs31 = synapse_utilities.compute_relative_dai_gs31( c1, c3, delta1, delta2, R2, dEs31, validation_flag );           % [V] Synaptic Reversal Potential.
                
            else                                                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs31; end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance for synapse 32 of a division after inversion subnetwork.
        function [ gs32, self ] = compute_dai_gs32( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                          % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                      % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                     % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                      % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                   % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ c1, c3, delta2, R1, R2, dEs31 ] = self.unpack_absolute_dai_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs32 = synapse_utilities.compute_absolute_dai_gs32( c1, c3, delta2, R1, R2, dEs31, validation_flag );              	% [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                           % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  c1, c3, delta1, delta2, R2, dEs31 ] = self.unpack_relative_dai_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs32 = synapse_utilities.compute_relative_dai_gs32( c1, c3, delta1, delta2, R2, dEs31, validation_flag );           % [V] Synaptic Reversal Potential.
                
            else                                                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs32; end
            
        end

        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 31 of a reduced division after inversion subnetwork.
        function [ gs31, self ] = compute_reduced_dai_gs31( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                  % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                              % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                             % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                              % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                           % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                       % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta1, delta2, R2, R3, Gm3, dEs31 ] = self.unpack_reduced_absolute_dai_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs31 = synapse_utilities.compute_reduced_absolute_dai_gs31( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                   % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  delta1, delta2, R2, R3, dEs31 ] = self.unpack_reduced_relative_dai_gs31_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs31 = synapse_utilities.compute_reduced_relative_dai_gs31( delta1, delta2, R2, R3, dEs31, validation_flag );              	% [V] Synaptic Reversal Potential.
                
            else                                                                                                                            % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs31; end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance for synapse 32 of a reduced division after inversion subnetwork.
        function [ gs32, self ] = compute_reduced_dai_gs32( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                  % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                              % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                             % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                              % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                           % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                       % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta1, delta2, R2, R3, Gm3, dEs31 ] = self.unpack_reduced_absolute_dai_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs32 = synapse_utilities.compute_reduced_absolute_dai_gs32( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                   % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [ delta1, delta2, R2, R3, Gm3, dEs31 ] = self.unpack_reduced_relative_dai_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs32 = synapse_utilities.compute_reduced_relative_dai_gs32( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );          % [V] Synaptic Reversal Potential.
                
            else                                                                                                                            % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs32; end
            
        end

        
        % ---------- Multiplication Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 41 of a multiplication subnetwork.
        function [ gs41, self ] = compute_multiplication_gs41( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                      % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                                  % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                                 % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                                  % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                               % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                           % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ c4, c6, delta1, delta2, R1, R3 ] = self.unpack_absolute_multiplication_gs41_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs41 = synapse_utilities.compute_absolute_multiplication_gs41( c4, c6, delta1, delta2, R1, R3, validation_flag );               % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                       % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  c4, c6, delta1, delta2, R3, dEs41 ] = self.unpack_relative_multiplication_gs41_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs41 = synapse_utilities.compute_relative_multiplication_gs41( c4, c6, delta1, delta2, R3, dEs41, validation_flag );            % [V] Synaptic Reversal Potential.
                
            else                                                                                                                                % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs41; end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance for synapse 32 of a multiplication subnetwork.
        function [ gs32, self ] = compute_multiplication_gs32( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                          % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                      % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                     % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                      % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                   % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta1, Gm3, dEs32, Ia3 ] = self.unpack_absolute_multiplication_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs32 = synapse_utilities.compute_absolute_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                           % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  delta1, Gm3, dEs32, Ia3 ] = self.unpack_relative_multiplication_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs32 = synapse_utilities.compute_relative_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );          % [V] Synaptic Reversal Potential.
                
            else                                                                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs32; end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance for synapse 43 of a multiplication subnetwork.
        function [ gs43, self ] = compute_multiplication_gs43( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                      % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                                  % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                                 % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                                  % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                               % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                           % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ c4, c6, delta2, R1, R3, dEs41 ] = self.unpack_absolute_multiplication_gs43_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs43 = synapse_utilities.compute_absolute_multiplication_gs43( c4, c6, delta2, R1, R3, dEs41, validation_flag );                % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                       % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  c4, c6, delta1, delta2, R3, dEs41 ] = self.unpack_relative_multiplication_gs43_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs43 = synapse_utilities.compute_relative_multiplication_gs43( c4, c6, delta1, delta2, R3, dEs41, validation_flag );            % [V] Synaptic Reversal Potential.
                
            else                                                                                                                                % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs43; end
            
        end

        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance for synapse 41 of a reduced multiplication subnetwork.
        function [ gs41, self ] = compute_reduced_multiplication_gs41( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                              % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                                          % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                                         % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                                          % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                                       % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                                   % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta1, delta2, R3, R4, Gm4, dEs41 ] = self.unpack_reduced_absolute_multiplication_gs41_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs41 = synapse_utilities.compute_reduced_absolute_multiplication_gs41( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                               % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  delta1, delta2, R3, R4, dEs41 ] = self.unpack_reduced_relative_multiplication_gs41_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs41 = synapse_utilities.compute_reduced_relative_multiplication_gs41( delta1, delta2, R3, R4, dEs41, validation_flag );              	% [V] Synaptic Reversal Potential.
                
            else                                                                                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs41; end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance for synapse 32 of a reduced multiplication subnetwork.
        function [ gs32, self ] = compute_reduced_multiplication_gs32( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                  % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                              % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                             % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                            	% [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                           % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                       % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta1, Gm3, dEs32, Ia3 ] = self.unpack_reduced_absolute_multiplication_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs32 = synapse_utilities.compute_reduced_absolute_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );          % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                  	% If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  delta1, Gm3, dEs32, Ia3 ] = self.unpack_reduced_relative_multiplication_gs32_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs32 = synapse_utilities.compute_reduced_relative_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );          % [V] Synaptic Reversal Potential.
                
            else                                                                                                                            % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs32; end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance for synapse 43 of a reduced multiplication subnetwork.
        function [ gs43, self ] = compute_reduced_multiplication_gs43( self, parameters, encoding_scheme, set_flag, validation_flag, synapse_utilities )
            
            % Set the default input arguments.
            if nargin < 6, synapse_utilities = self.synapse_utilities; end                                                                              % [class] Synapse Utilities.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end                                                                          % [T/F] Validation Flag.
            if nargin < 4, set_flag = true; end                                                                                                         % [T/F] Set Flag.
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                                          % [str] Encoding Scheme.
            if nargin < 2, parameters = {  }; end                                                                                                       % [-] Input Parameters Cell.
            
            % Determine how to compute the synaptic conductance for an division subnetwork.
            if strcmpi( encoding_scheme, 'absolute' )                                                                                                   % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the synaptic conductance for an absolute division subnetwork.
                [ delta1, delta2, R3, R4, Gm4, dEs41 ] = self.unpack_reduced_absolute_multiplication_gs43_parameters( parameters );
                
                % Compute the synaptic conductance for an absolue division subnetwork.
                gs43 = synapse_utilities.compute_reduced_absolute_multiplication_gs43( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );           % [V] Synaptic Reversal Potential.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                                                               % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the synaptic conductance for a relative division subnetwork.
                [  delta1, delta2, R3, R4, Gm4, dEs41 ] = self.unpack_reduced_relative_multiplication_gs43_parameters( parameters );
                
                % Compute the synaptic conductance for a relative division subnetwork.
                gs43 = synapse_utilities.compute_reduced_relative_multiplication_gs43( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );       	% [V] Synaptic Reversal Potential.
                
            else                                                                                                                                        % Otherwise...
                
                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
                
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs43; end
            
        end

        
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance of a driven multistate cpg subnetwork.
        function [ gs, self ] = compute_dmcpg_gs( self, dEs, delta_oscillatory, Id_max, set_flag, synapse_utilities )
            
            % Define the default input arguments.
            if nargin < 5, synapse_utilities = self.synapse_utilities; end                              % [class] Synapse Utilities.
            if nargin < 4, Id_max = self.Id_max_DEFAULT; end                                            % [A] Maximum Drive Current.
            if nargin < 3, delta_oscillatory = self.delta_oscillatory_DEFAULT; end                      % [-] Oscillatory Delta.
            if nargin < 2, dEs = self.dEs; end                                                          % [V] Synaptic Reversal Potential.
            
            % Compute the maximum synaptic conductance.
            gs = synapse_utilities.compute_dmcpg_gsynmax( dEs, delta_oscillatory, Id_max );             % [S] Maximum Synaptic Conductance
            
            % Determine whether to update the synapse object.
            if set_flag, self.gs = gs; end
            
        end
        
        
        %% Enable & Disable Functions
        
        % Implement a function to toogle whether this syanpse is enabled.
        function [ enabled_flag, self ] = toggle_enabled( self, enabled_flag, set_flag )
            
            % Set the default input arguments.
            if nargin < 3, set_flag = true; end                             % [T/F] Set Flag.
            if narign < 2, enabled_flag = self.enabled_flag; end            % [T/F] Enabled Flag.
            
            % Toggle whether the neuron is enabled.
            enabled_flag = ~enabled_flag;                                   % [T/F] Synapse Enabled Flag.
            
            % Determine whether to update the neuron object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        % Implement a function to enable this neuron.
        function [ enabled_flag, self ] = enable( self, set_flag )
            
            % Set the default input arguments.
            if nargin < 2, set_flag = true; end                     % [T/F] Set Flag.
            
            % Enable this neuron.
            enabled_flag = true;                                   	% [T/F] Synapse Enabled Flag.
            
            % Determine whether to update the neuron object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        % Implement a function to disable this neuron.
        function [ enabled_flag, self ] = disable( self, set_flag )
            
            % Set the default input arguments.
            if anrgin < 2, set_flag = true; end                     % [T/F] Set Flag.
            
            % Disable this neuron.
            enabled_flag = false;                                  	% [T/F] Synapse Enabled Flag.
            
            % Determine wehther to update the neuron object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        %% Print Functions.
        
        % Implement a function to print the synapse information.
        function print( self, ID, name, dEs, gs, Gs, from_neuron_ID, to_neuron_ID, delta, enabled_flag, verbose_flag )
           
            % Set the default input arguments.
            if nargin < 11, verbose_flag = false; end
            if nargin < 10, enabled_flag = self.enabled_flag; end
            if nargin < 9, delta = self.delta; end
            if nargin < 8, to_neuron_ID = self.to_neuron_ID; end
            if nargin < 7, from_neuron_ID = self.from_neuron_ID; end
            if nargin < 6, Gs = self.Gs; end
            if nargin < 5, gs = self.gs; end
            if nargin < 4, dEs = self.dEs; end
            if nargin < 3, name = self.name; end
            if nargin < 2, ID = self.ID; end
            
            % Print out the information for this neuron.
            synapse_utilies.print( ID, name, dEs, gs, Gs, from_neuron_ID, to_neuron_ID, delta, enabled_flag, verbose_flag );
            
        end
        
        
        %% Save & Load Functions.
        
        % Implement a function to save synapse data as a matlab object.
        function save( self, directory, file_name, synapse )
            
            % Set the default input arguments.
            if nargin < 4, synapse = self; end                      % [class] Synapse Object to Save.
            if nargin < 3, file_name = 'Synapse.mat'; end           % [str] File Name to Save.
            if nargin < 2, directory = '.'; end                     % [str] Save Directory.
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Save the neuron data.
            save( full_path, synapse )
            
        end
        
        
        % Implement a function to load synapse data as a matlab object.
        function synapse = load( ~, directory, file_name )
            
            % Set the default input arguments.
            if nargin < 3, file_name = 'Synapse.mat'; end           % [str] File Name to Load.
            if nargin < 2, directory = '.'; end                     % [str] Load Directory.
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];
            
            % Load the data.
            data = load( full_path );
            
            % Retrieve the synapse object from the loaded data object.
            synapse = data.synapse;
            
        end
        
        
    end
end

