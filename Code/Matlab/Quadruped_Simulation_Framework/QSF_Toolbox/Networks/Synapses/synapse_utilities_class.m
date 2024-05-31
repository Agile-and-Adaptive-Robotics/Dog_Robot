classdef synapse_utilities_class
    
    % This class contains properties and methods related to synapse utilities.
    
    
    %% SYNAPSE UTILITIES PROPERTIES
    
    % Define the class properties.
    properties
        
        
        
    end
    
    
    % Define private, constant class properties.
    properties ( Access = private, Constant = true )
        
        % Define the neuron parameters.
        R_DEFAULT = 20e-3;                                     	% [V] Activation Domain.
        Gm_DEFAULT = 1e-6;                                  	% [S] Membrane Conductance.
        
        % Define the synaptic reversal potential parameters.
        dEs_max_DEFAULT = 194e-3;                             	% [V] Maximum Synaptic Reversal Potential.
        dEs_min_DEFAULT = -194e-3;                             	% [V] Minimum Synaptic Reversal Potential. (-40e-3 [mV] is a common choice.)
        dEs_small_negative_DEFAULT = -1e-3;                   	% [V] Small Negative Synaptic Reversal Potential.
        
        % Define the applied current parameters.
        Id_max_DEFAULT = 1.25e-9;                              	% [A] Maximum Drive Current.
        Ia_absolute_addition_DEFAULT = 0;                      	% [A] Absolute Addition Applied Current.
        Ia_relative_addition_DEFAULT = 0;                      	% [A] Relative Addition Applied Current.
        Ia_absolute_subtraction_DEFAULT = 0;                  	% [A] Absolute Subtraction Applied Current.
        Ia_relative_subtraction_DEFAULT = 0;                  	% [A] Relative Subtraction Applied Current.
        Ia1_absolute_inversion_DEFAULT = 0;                    	% [A] Absolute Inversion Applied Current 1.
        Ia2_absolute_inversion_DEFAULT = 2e-8;              	% [A] Absolute Inversion Applied Current 2.
        Ia1_relative_inversion_DEFAULT = 0;                  	% [A] Relative Inversion Applied Current 1.
        Ia2_relative_inversion_DEFAULT = 2e-8;                	% [A] Relative Inversion Applied Current 2.
        Ia_absolute_division_DEFAULT = 0;                      	% [A] Absolute Division Applied Current.
        Ia_relative_division_DEFAULT = 0;                    	% [A] Relative Division Applied Current.
        
        % Define the CPG parameters.
        delta_oscillatory_DEFAULT = 0.01e-3;                  	% [-] CPG Oscillatory Delta.
        delta_bistable_DEFAULT = -10e-3;                     	% [-] CPG Bistable Delta.
        delta_noncpg_DEFAULT = 0;                            	% [-] CPG Delta.

        % Define the subnetwork design parameters.
        c_absolute_addition_DEFAULT = 1;                      	% [-] Absolute Addition Subnetwork Gain.
        c_relative_addition_DEFAULT = 1;                       	% [-] Relative Addition Subnetwork Gain.
        c_absolute_subtraction_DEFAULT = 1;                    	% [-] Absolute Subtraction Subnetwork Gain.
        c_relative_subtraction_DEFAULT = 1;                    	% [-] Relative Subtraction Subnetwork Gain.
        c_absolute_inversion_DEFAULT = 1;                      	% [-] Absolute Inversion Subnetwork Gain.
        c_relative_inversion_DEFAULT = 1;                      	% [-] Relative Inversion Subnetwork Gain.
        c_absolute_division_DEFAULT = 1;                       	% [-] Absolute Division Subnetwork Gain.
        c_relative_division_DEFAULT = 1;                    	% [-] Relative Division Subnetwork Gain.
        c_absolute_multiplication_DEFAULT = 1;                	% [-] Absolute Multiplication Subnetwork Gain.
        c_relative_multiplication_DEFAULT = 1;                  % [-] Relative Multiplication Subnetwork Gain.
        
        % Define the subnetwork neuron numbers.
        num_addition_neurons_DEFAULT = 3;                      	% [#] Number of Addition Neurons.
        
        % Define inversion and division subnetwork parameters.
        epsilon_DEFAULT = 1e-6;                               	% [V] Subnetwork Input Offset.
        delta_DEFAULT = 1e-6;                                 	% [V] Subnetwork Output Offset.
        alpha_DEFAULT = 1e-6;                                 	% [V] Division Subnetwork Denominator Offset.
        
        % Set the default flags.
        validation_flag_DEFAULT = true;                         % [T/F] Validation Flag. (Determines whether to validate computed values.)
        
    end
    
    
    %% SYNAPSE UTILITIES METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = synapse_utilities_class(  )
            
            
            
        end
        
        
        %% Name Functions.
        
        % Implement a function to generate a name from an ID.
        function name = ID2name( ~, ID )
            
            % Generate a name for the synapse.
            name = sprintf( 'Synapse %s', ID );
            
        end
        
        
        % Implement a function to generate names from IDs.
        function names = IDs2names( self, IDs )
        
            % Compute the number of IDs.
            num_IDs = length( IDs );
            
            % Preallocate a cell array to store the names.
            names = cell( 1, num_IDs );
            
            % Generate a name for each ID.
            for k = 1:num_IDs                 % Iterate through each of the IDs...
                
                % Generate the name associated with this ID.
                names{ k } = self.ID2name( IDs( k ) );
                
            end
            
        end
        
        
        %% Synaptic Validation Functions.
        
        % Implement a function to validate synaptic conductance results.
        function valid_flag = validate_gs( self, gs )

            % Set the default input arguments.
            if nargin < 2, gs = self.gs_DEFAULT; end
            
            % Determine whether the provided synaptic conductances are valid.
            if all( gs >= 0 )           % If all of the synaptic conductances are greater than or equal to zero...
               
                % Set the valid flag to true;
                valid_flag = true;
                
            else                        % Otherwise....
                
                % Set the valid flag to false.
                valid_flag = false;
                
            end

        end
        
        
        %% Synaptic Reversal Potential Compute Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential for an absolute transmission subnetwork.
        function dEs21 = compute_absolute_transmission_dEs21( self )
            
            % Compute the synaptic reversal potential.
            dEs21 = self.dEs_max_DEFAULT;                                     % [V] Synaptic Reversal Potential.
            
        end
        
                
        % Implement a function to compute the synaptic reversal potential for a relative transmission subnetwork.
        function dEs21 = compute_relative_transmission_dEs21( self )
            
            % Compute the synaptic reversal potential.
            dEs21 = self.dEs_max_DEFAULT;                                     % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a transmission subnetwork.
        function dEs21 = compute_transmission_dEs21( self, encoding_scheme )
           
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs21 = self.compute_absolute_transmission_dEs21(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs21 = self.compute_relative_transmission_dEs21(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential for absolute addition subnetwork synapses.
        function dEs = compute_absolute_addition_dEs( self )
            
            % Compute the synaptic reversal potential.
            dEs = self.dEs_max_DEFAULT;                                     % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative addition subnetwork synapses.
        function dEs = compute_relative_addition_dEs( self )
            
            % Compute the synaptic reversal potential.
            dEs = self.dEs_max_DEFAULT;                                     % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for addition subnetwork synapses.
        function dEs = compute_addition_dEs( self, encoding_scheme )
           
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs = self.compute_absolute_addition_dEs(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs = self.compute_relative_addition_dEs(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions (Excitatory) ----------
        
        % Implement a function to compute the synaptic reversal potential for absolute subtraction subnetwork excitatory synapses.
        function dEs = compute_absolute_subtraction_dEs_excitatory( self )
            
            % Compute the synaptic reversal potential.
            dEs = self.dEs_max_DEFAULT;                                     % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative subtraction subnetwork excitatory synapses.
        function dEs = compute_relative_subtraction_dEs_excitatory( self )
            
            % Compute the synaptic reversal potential.
            dEs = self.dEs_max_DEFAULT;                                     % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for subtraction subnetwork excitatory synapses.
        function dEs = compute_subtraction_dEs_excitatory( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs = self.compute_absolute_subtraction_dEs_excitatory(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs = self.compute_relative_subtraction_dEs_excitatory(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions (Inhibitory) ----------
        
        % Implement a function to compute the synaptic reversal potential for absolute subtraction subnetwork inhibitory synapses.
        function dEs = compute_absolute_subtraction_dEs_inhibitory( self )
            
            % Compute the synaptic reversal potential.
%             dEs = self.dEs_min_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            dEs = self.dEs_min_DEFAULT;                                    % [V] Synaptic Reversal Potential.

        end
        
                
        % Implement a function to compute the synaptic reversal potential for relative subtraction subnetwork inhibitory synapses.
        function dEs = compute_relative_subtraction_dEs_inhibitory( self )
            
            % Compute the synaptic reversal potential.
%             dEs = self.dEs_min_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            dEs = self.dEs_min_DEFAULT;                                    % [V] Synaptic Reversal Potential.

        end
                
        
        % Implement a function to compute the synaptic reversal potential for subtraction subnetwork inhibitory synapses.
        function dEs = compute_subtraction_dEs_inhibitory( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs = self.compute_absolute_subtraction_dEs_inhibitory(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs = self.compute_relative_subtraction_dEs_inhibitory(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
                
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential for absolute inversion subnetwork synapses.
        function dEs21 = compute_absolute_inversion_dEs21( ~ )
            
            % Compute the synaptic reversal potential.
            dEs21 = 0;            % [V] Synaptic Reversal Potential.

        end
         
        
        % Implement a function to compute the synaptic reversal potential for relative inversion subnetwork synapses.
        function dEs21 = compute_relative_inversion_dEs21( ~ )
            
            % Compute the synaptic reversal potential.
            dEs21 = 0;         	% [V] Synaptic Reversal Potential.

        end
        
        
        % Implement a function to compute the synaptic reversal potential for inversion subnetwork synapses.
        function dEs21 = compute_inversion_dEs21( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs21 = self.compute_absolute_inversion_dEs21(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs21 = self.compute_relative_inversion_dEs21(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential for reduced absolute inversion subnetwork synapses.
        function dEs21 = compute_reduced_absolute_inversion_dEs21( ~ )
           
            % Compute the synaptic reversal potential.
            dEs21 = 0;          % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced relative inversion subnetwork synapses.
        function dEs21 = compute_reduced_relative_inversion_dEs21( ~ )
            
            % Compute the synaptic reversal potential.
            dEs21 = 0;          % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced inversion subnetwork synapses.
        function dEs21 = compute_reduced_inversion_dEs21( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs21 = self.compute_reduced_absolute_inversion_dEs21(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs21 = self.compute_reduced_relative_inversion_dEs21(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to compute the synaptic reversal potential for absolute division numerator synapses.
        function dEs31 = compute_absolute_division_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative division numerator synapses.
        function dEs31 = compute_relative_division_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for division subnetwork numerator synapses.
        function dEs31 = compute_division_dEs31( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs31 = self.compute_absolute_division_dEs31(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs31 = self.compute_relative_division_dEs31(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division Subnetwork Functions (Synapse 32) ----------
        
        % Implement a function to compute the synaptic reversal potential for absolute division denominator synapses.
        function dEs32 = compute_absolute_division_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
            % dEs = self.dEs_small_negative_DEFAULT;                      	% [V] Synaptic Reversal Potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative division denominator synapses.
        function dEs32 = compute_relative_division_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
%             dEs = self.dEs_small_negative_DEFAULT;                     	% [V] Synaptic Reversal Potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end

        
        % Implement a function to compute the synaptic reversal potential for division subnetwork denominator synapses.
        function dEs32 = compute_division_dEs32( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs32 = self.compute_absolute_division_dEs32(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs32 = self.compute_relative_division_dEs32(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division Subnetwork Functions (Combined) ----------

        % Implement a function to compute the synaptic reversal potential for absolute division synapses.
        function [ dEs31, dEs32 ] = compute_absolute_division_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_absolute_division_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_absolute_division_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative division synapses.
        function [ dEs31, dEs32 ] = compute_relative_division_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_relative_division_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_relative_division_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for division synapes.
        function [ dEs31, dEs32 ] = compute_division_dEs( self, encoding_scheme )
        
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_division_dEs31( encoding_scheme );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_division_dEs32( encoding_scheme );
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Synapse 31) ----------

        % Implement a function to compute the synaptic reversal potential for reduced absolute division numerator synapses.
        function dEs31 = compute_reduced_absolute_division_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced relative division numerator synapses.
        function dEs31 = compute_reduced_relative_division_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced division subnetwork numerator synapses.
        function dEs31 = compute_reduced_division_dEs31( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs31 = self.compute_reduced_absolute_division_dEs31(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs31 = self.compute_reduced_relative_division_dEs31(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Synapse 32) ----------
        
        % Implement a function to compute the synaptic reversal potential for reduced absolute division denominator synapses.
        function dEs32 = compute_reduced_absolute_division_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end
        
                
        % Implement a function to compute the synaptic reversal potential for reduced relative division denominator synapses.
        function dEs32 = compute_reduced_relative_division_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced division subnetwork denominator synapses.
        function dEs32 = compute_reduced_division_dEs32( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs32 = self.compute_reduced_absolute_division_dEs32(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs32 = self.compute_reduced_relative_division_dEs32(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Combined) ----------

        % Implement a function to compute the synaptic reversal potential for absolute division synapses.
        function [ dEs31, dEs32 ] = compute_reduced_absolute_division_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_reduced_absolute_division_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_reduced_absolute_division_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative division synapses.
        function [ dEs31, dEs32 ] = compute_reduced_relative_division_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_reduced_relative_division_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_reduced_relative_division_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for division synapes.
        function [ dEs31, dEs32 ] = compute_reduced_division_dEs( self, encoding_scheme )
        
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_reduced_division_dEs31( encoding_scheme );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_reduced_division_dEs32( encoding_scheme );
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions (Synapse 31) ----------

        % Implement a function to compute the synaptic reversal potential for absolute division after inversion synapse 31.
        function dEs31 = compute_absolute_dai_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative division after inversion synapse 31.
        function dEs31 = compute_relative_dai_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for division after inversion subnetwork synapse 31.
        function dEs31 = compute_dai_dEs31( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs31 = self.compute_absolute_dai_dEs31(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs31 = self.compute_relative_dai_dEs31(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions (Synapse 32) ----------

        % Implement a function to compute the synaptic reversal potential for absolute division after inversion synapse 32.
        function dEs32 = compute_absolute_dai_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative division after inversion synapse 32.
        function dEs32 = compute_relative_dai_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end

        
        % Implement a function to compute the synaptic reversal potential for division after inversion subnetwork synapse 32.
        function dEs32 = compute_dai_dEs32( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs32 = self.compute_absolute_dai_dEs32(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs32 = self.compute_relative_dai_dEs32(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions (Combined) ----------

        % Implement a function to compute the synaptic reversal potential for absolute division after inversion synapses.
        function [ dEs31, dEs32 ] = compute_absolute_dai_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_absolute_dai_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_absolute_dai_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for relative division after inversion synapses.
        function [ dEs31, dEs32 ] = compute_relative_dai_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_relative_dai_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_relative_dai_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for division after inversion synapes.
        function [ dEs31, dEs32 ] = compute_dai_dEs( self, encoding_scheme )
        
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_dai_dEs31( encoding_scheme );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_dai_dEs32( encoding_scheme );
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Synapse 31) ----------

        % Implement a function to compute the synaptic reversal potential for reduced absolute division after inversion synapse 31.
        function dEs31 = compute_reduced_absolute_dai_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced relative division after inversion synapse 31.
        function dEs31 = compute_reduced_relative_dai_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                   % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced division after inversion subnetwork synapse 31.
        function dEs31 = compute_reduced_dai_dEs31( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs31 = self.compute_reduced_absolute_dai_dEs31(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs31 = self.compute_reduced_relative_dai_dEs31(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Synapse 32) ----------

        % Implement a function to compute the synaptic reversal potential for reduced absolute division after inversion synapse 32.
        function dEs32 = compute_reduced_absolute_dai_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced relative division after inversion synapse 32.
        function dEs32 = compute_reduced_relative_dai_dEs32( ~ )
            
            % Compute the synaptic reversal potential.
            dEs32 = 0;                                                        % [V] Synaptic Reversal Potential.

        end

        
        % Implement a function to compute the synaptic reversal potential for reduced division after inversion subnetwork synapse 32.
        function dEs32 = compute_reduced_dai_dEs32( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs32 = self.compute_reduced_absolute_dai_dEs32(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs32 = self.compute_reduced_relative_dai_dEs32(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Combined) ----------
        
        % Implement a function to compute the synaptic reversal potential for reduced absolute division after inversion synapses.
        function [ dEs31, dEs32 ] = compute_reduced_absolute_dai_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_reduced_absolute_dai_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_reduced_absolute_dai_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced relative division after inversion synapses.
        function [ dEs31, dEs32 ] = compute_reduced_relative_dai_dEs( self )
        
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_reduced_relative_dai_dEs31(  );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_reduced_relative_dai_dEs32(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for reduced division after inversion synapes.
        function [ dEs31, dEs32 ] = compute_reduced_dai_dEs( self, encoding_scheme )
        
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Compute the synaptic reversal potential of the numerator synapse.
            dEs31 = self.compute_reduced_dai_dEs31( encoding_scheme );
            
            % Compute the synaptic reversal potential of the denominator synapse.
            dEs32 = self.compute_reduced_dai_dEs32( encoding_scheme );
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 41) ----------
        
        % Implement a function to compute the synaptic reversal potential for an absolute multiplication subnetwork (synapse 41).
        function dEs41 = compute_absolute_multiplication_dEs41( self )
            
            % Compute the synaptic reversal potential.
            dEs41 = self.compute_absolute_dai_dEs31( self );                                    % [V] Synaptic Reversal Potential.
            
        end
        

        % Implement a function to compute the synaptic reversal potential for a relative multiplication subnetwork (synapse 41).
        function dEs41 = compute_relative_multiplication_dEs41( self )
            
            % Compute the synaptic reversal potential.
            dEs41 = self.compute_relative_dai_dEs31( self );                                    % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a multiplication subnetwork (synapse 41).
        function dEs41 = compute_multiplication_dEs41( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs41 = self.compute_absolute_multiplication_dEs41(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs41 = self.compute_relative_multiplication_dEs41(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 32) ----------
        
        % Implement a function to compute the synaptic reversal potential for an absolute multiplication subnetwork (synapse 32).
        function dEs32 = compute_absolute_multiplication_dEs32( self )
            
            % Compute the synaptic reversal potential.
            dEs32 = self.compute_absolute_inversion_dEs21(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a relative multiplication subnetwork (synapse 32).
        function dEs32 = compute_relative_multiplication_dEs32( self )
            
            % Compute the synaptic reversal potential.
            dEs32 = self.compute_relative_inversion_dEs21(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a multiplication subnetwork (synapse 32).
        function dEs32 = compute_multiplication_dEs32( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs32 = self.compute_absolute_multiplication_dEs32(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs32 = self.compute_relative_multiplication_dEs32(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 43) ----------
        
        % Implement a function to compute the synaptic reversal potential for an absolute multiplication subnetwork (synapse 43).
        function dEs43 = compute_absolute_multiplication_dEs43( self )
            
            % Compute the synaptic reversal potential.
            dEs43 = self.compute_absolute_dai_dEs32(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a relative multiplication subnetwork (synapse 43).
        function dEs43 = compute_relative_multiplication_dEs43( self )
            
            % Compute the synaptic reversal potential.
            dEs43 = self.compute_relative_dai_dEs32(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a multiplication subnetwork (synapse 43).
        function dEs43 = compute_multiplication_dEs43( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs43 = self.compute_absolute_multiplication_dEs43(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs43 = self.compute_relative_multiplication_dEs43(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Combined) ----------

        % Implement a function to compute the synaptic reversal potential for an absolute multiplication subnetwork (combined).
        function [ dEs41, dEs32, dEs43 ] = compute_absolute_multiplication_dEs( self )
            
            % Compute the synaptic reversal potential of the 41 synapse.
            dEs41 = self.compute_absolute_multiplication_dEs41(  );
            
            % Compute the synaptic reversal potential of the 32 synapse.
            dEs32 = self.compute_absolute_multiplication_dEs32(  );
            
            % Compute the synaptic reversal potential of the 43 synapse.
            dEs43 = self.compute_absolute_multiplication_dEs43(  );
            
        end
                
        
        % Implement a function to compute the synaptic reversal potential for a relative multiplication subnetwork (combined).
        function [ dEs41, dEs32, dEs43 ] = compute_relative_multiplication_dEs( self )
            
            % Compute the synaptic reversal potential of the 41 synapse.
            dEs41 = self.compute_relative_multiplication_dEs41(  );
            
            % Compute the synaptic reversal potential of the 32 synapse.
            dEs32 = self.compute_relative_multiplication_dEs32(  );
            
            % Compute the synaptic reversal potential of the 43 synapse.
            dEs43 = self.compute_relative_multiplication_dEs43(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a multiplication subnetwork (combined).
        function [ dEs41, dEs32, dEs43 ] = compute_multiplication_dEs( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                [ dEs41, dEs32, dEs43 ] = self.compute_absolute_multiplication_dEs(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                [ dEs41, dEs32, dEs43 ] = self.compute_relative_multiplication_dEs(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 41) ----------
        
        % Implement a function to compute the synaptic reversal potential for a reduced absolute multiplication subnetwork (synapse 41).
        function dEs41 = compute_reduced_absolute_multiplication_dEs41( self )
            
            % Compute the synaptic reversal potential.
            dEs41 = self.compute_reduced_absolute_dai_dEs31( self );                                    % [V] Synaptic Reversal Potential.
            
        end
        

        % Implement a function to compute the synaptic reversal potential for a reduced relative multiplication subnetwork (synapse 41).
        function dEs41 = compute_reduced_relative_multiplication_dEs41( self )
            
            % Compute the synaptic reversal potential.
            dEs41 = self.compute_reduced_relative_dai_dEs31( self );                                    % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a reduced multiplication subnetwork (synapse 41).
        function dEs41 = compute_reduced_multiplication_dEs41( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs41 = self.compute_reduced_absolute_multiplication_dEs41(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs41 = self.compute_reduced_relative_multiplication_dEs41(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 32) ----------
        
        % Implement a function to compute the synaptic reversal potential for a reduced absolute multiplication subnetwork (synapse 32).
        function dEs32 = compute_reduced_absolute_multiplication_dEs32( self )
            
            % Compute the synaptic reversal potential.
            dEs32 = self.compute_reduced_absolute_inversion_dEs21(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a reduced relative multiplication subnetwork (synapse 32).
        function dEs32 = compute_reduced_relative_multiplication_dEs32( self )
            
            % Compute the synaptic reversal potential.
            dEs32 = self.compute_reduced_relative_inversion_dEs21(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a reduced multiplication subnetwork (synapse 32).
        function dEs32 = compute_reduced_multiplication_dEs32( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs32 = self.compute_reduced_absolute_multiplication_dEs32(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs32 = self.compute_reduced_relative_multiplication_dEs32(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 43) ----------
        
        % Implement a function to compute the synaptic reversal potential for a reduced absolute multiplication subnetwork (synapse 43).
        function dEs43 = compute_reduced_absolute_multiplication_dEs43( self )
            
            % Compute the synaptic reversal potential.
            dEs43 = self.compute_reduced_absolute_dai_dEs32(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a reduced relative multiplication subnetwork (synapse 43).
        function dEs43 = compute_reduced_relative_multiplication_dEs43( self )
            
            % Compute the synaptic reversal potential.
            dEs43 = self.compute_reduced_relative_dai_dEs32(  );                         % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a reduced multiplication subnetwork (synapse 43).
        function dEs43 = compute_reduced_multiplication_dEs43( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                dEs43 = self.compute_reduced_absolute_multiplication_dEs43(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                dEs43 = self.compute_reduced_relative_multiplication_dEs43(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Combined) ----------

        % Implement a function to compute the synaptic reversal potential for a reduced absolute multiplication subnetwork (combined).
        function [ dEs41, dEs32, dEs43 ] = compute_reduced_absolute_multiplication_dEs( self )
            
            % Compute the synaptic reversal potential of the 41 synapse.
            dEs41 = self.compute_reduced_absolute_multiplication_dEs41(  );
            
            % Compute the synaptic reversal potential of the 32 synapse.
            dEs32 = self.compute_reduced_absolute_multiplication_dEs32(  );
            
            % Compute the synaptic reversal potential of the 43 synapse.
            dEs43 = self.compute_reduced_absolute_multiplication_dEs43(  );
            
        end
                
        
        % Implement a function to compute the synaptic reversal potential for a reduced relative multiplication subnetwork (combined).
        function [ dEs41, dEs32, dEs43 ] = compute_reduced_relative_multiplication_dEs( self )
            
            % Compute the synaptic reversal potential of the 41 synapse.
            dEs41 = self.compute_reduced_relative_multiplication_dEs41(  );
            
            % Compute the synaptic reversal potential of the 32 synapse.
            dEs32 = self.compute_reduced_relative_multiplication_dEs32(  );
            
            % Compute the synaptic reversal potential of the 43 synapse.
            dEs43 = self.compute_reduced_relative_multiplication_dEs43(  );
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a reduced multiplication subnetwork (combined).
        function [ dEs41, dEs32, dEs43 ] = compute_reduced_multiplication_dEs( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                [ dEs41, dEs32, dEs43 ] = self.compute_reduced_absolute_multiplication_dEs(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                [ dEs41, dEs32, dEs43 ] = self.compute_reduced_relative_multiplication_dEs(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end

        
        % ---------- Derivation Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential for a derivation subnetwork.
        function dEs31 = compute_derivation_dEs31( self )
            
            % Compute the synaptic reversal potential.
            dEs31 = self.dEs_max_DEFAULT;                                    % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a derivation subnetwork.
        function dEs32 = compute_derivation_dEs32( self )
            
            % Compute the synaptic reversal potential.
            dEs32 = self.dEs_min_DEFAULT;                                	% [V] Synaptic Reversal Potential.
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------

        %{
        
        % Implement a function to compute the synaptic reversal potentials for an integration subnetwork.
        function dEs = compute_integration_dEs( ~, Gm, R, gs )
        
            %{
            Input(s):
                Gm      =   [S] Membrane Conductance.
                R       =   [V] Maximum Membrane Voltage.
                gs      =   [S] Maximum Synaptic Conductance.
            
            Output(s):
                dEs   =   [S] Synaptic Reversal Potential.
            %}
            
           % Compute the synaptic reversal potentials for an integration subnetwork.
           dEs = - ( Gm*R )./gs;
            
        end
        
        %}
        
        % Implement a function to compute the synaptic reversal potential for a integration subnetwork.
        function dEs1 = compute_integration_dEs1( self )
            
            % Compute the synaptic reversal potential.
            dEs1 = self.dEs_min_DEFAULT;                                    % [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a integration subnetwork.
        function dEs2 = compute_integration_dEs2( self )
            
            % Compute the synaptic reversal potential.
            dEs2 = self.dEs_min_DEFAULT;                                	% [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a voltage based integration subnetwork.
        function dEs1 = compute_vbi_dEs1( self )
            
            % Compute the synaptic reversal potential.
            dEs1 = self.dEs_max_DEFAULT;                                 	% [V] Synaptic Reversal Potential.
            
        end
        
        
        % Implement a function to compute the synaptic reversal potential for a voltage based integration subnetwork.
        function dEs2 = compute_vbi_dEs2( self )
            
            % Compute the synaptic reversal potential.
            dEs2 = -self.dEs_max_DEFAULT;                                  	% [V] Synaptic Reversal Potential.
            
        end
        
        
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to compute the synaptic reversal potential for a driven multistate cpg subnetwork.
        function dEs = compute_dmcpg_dEs( self )
            
            % Compute the synaptic reversal potential.
            dEs = self.dEs_max_DEFAULT;                                     % [V] Synaptic Reversal Potential.
            
        end
        
        
        %% Maximum Synaptic Conductance Compute Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance of absolute transmission subnetwork synapses.
        function gs21 = compute_absolute_transmission_gs21( self, R2, Gm2, dEs21, Ia2, validation_flag )
        
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia2 = self.Ia_DEFAULT; end
            if nargin < 4, dEs21 = self.dEs_DEFAULT; end
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, R2 = self.R_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs21 = ( R2*Gm2 - Ia2 )/( dEs21 - R2 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs21 ), 'Invalid synaptic conductance detected.' )
            
            end
                
        end
        
        
        % Implement a function to compute the maximm synaptic conductance of relative transmission subnetwork synapses.
        function gs21 = compute_relative_transmission_gs21( self, R2, Gm2, dEs21, Ia2, validation_flag )
            
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia2 = self.Ia_DEFAULT; end
            if nargin < 4, dEs21 = self.dEs_DEFAULT; end
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, R2 = self.R_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs21 = ( R2*Gm2 - Ia2 )/( dEs21 - R2 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs21 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of transmission subnetwork synapses.
        function gs21 = compute_transmission_gs21( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                R2 = parameters{ 1 };
                Gm2 = parameters{ 2 };
                dEs21 = parameters{ 3 };
                Ia2 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs21 = self.compute_absolute_transmission_gs21( R2, Gm2, dEs21, Ia2, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                R2 = parameters{ 1 };
                Gm2 = parameters{ 2 };
                dEs21 = parameters{ 3 };
                Ia2 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs21 = self.compute_relative_transmission_gs21( R2, Gm2, dEs21, Ia2, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
            
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance of absolute addition subnetwork synapses.
        function gs_nk = compute_absolute_addition_gs( self, c_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag )
            
            % Define the default input arguments.
            if nargin < 7, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 6, Ia_n = self.Ia_absolute_addition_DEFAULT; end                                % [A] Applied Current.
            if nargin < 5, dEs_nk = self.dEs_max_DEFAULT; end                                           % [V] Synaptic Reversal Potential.
            if nargin < 4, Gm_n = self.Gm_DEFAULT; end                                                  % [S] Membrane Conductance.
            if nargin < 3, R_k = self.R_DEFAULT; end                                                    % [V] Activation Domain.
            if nargin < 2, c_k = self.c_absolute_addition_DEFAULT; end                               	% [-] Absolute Addition Subnetwork Gain.
            
            % Compute the maximum synaptic conductance.
            gs_nk = ( Ia_n - c_k.*R_k.*Gm_n )./( c_k.*R_k - dEs_nk );                                 	% [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs_nk ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of relative addition subnetwork synapses.
        function gs_nk = compute_relative_addition_gs( self, c_k, R_n, Gm_n, dEs_nk, Ia_n, validation_flag )
            
            % Define the default input arguments.
            if nargin < 7, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 6, Ia_n = self.Ia_relative_addition_DEFAULT; end                                % [A] Applied Current.
            if nargin < 5, dEs_nk = self.dEs_max_DEFAULT; end                                           % [V] Synaptic Reversal Potential.
            if nargin < 4, Gm_n = self.Gm_DEFAULT; end                                                  % [S] Membrane Conductance.
            if nargin < 3, R_n = self.R_DEFAULT; end                                                    % [V] Activation Domain.
            if nargin < 2, c_k = self.c_relative_addition_DEFAULT; end                                	% [-] Absolute Addition Subnetwork Gain.
            
            % Compute the maximum synaptic conductance.
            gs_nk = ( Ia_n - c_k*R_n*Gm_n )/( c_k*R_n - dEs_nk );                                       % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs_nk ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of addition subnetwork synapses.
        function gs_nk = compute_addition_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c_k = parameters{ 1 };
                R_k = parameters{ 2 };
                Gm_n = parameters{ 3 };
                dEs_nk = parameters{ 4 };
                Ia_n = parameters{ 5 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs_nk = self.compute_absolute_addition_gs( c_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c_k = parameters{ 1 };
                R_n = parameters{ 2 };
                Gm_n = parameters{ 3 };
                dEs_nk = parameters{ 4 };
                Ia_n = parameters{ 5 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs_nk = compute_relative_addition_gs( c_k, R_n, Gm_n, dEs_nk, Ia_n, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance of absolute subtraction subnetwork synapses.
        function gs_nk = compute_absolute_subtraction_gs( self, c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag )
            
            % Define the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia_n = self.Ia_absolute_subtraction_DEFAULT; end                             % [A] Applied Current.
            if nargin < 6, dEs_nk = self.dEs_max_DEFAULT; end                                           % [V] Synaptic Reversal Potential.
            if nargin < 5, Gm_n = self.Gm_DEFAULT; end                                                  % [S] Membrane Conductance.
            if nargin < 4, R_k = self.R_DEFAULT; end                                                    % [V] Activation Domain.
            if nargin < 3, s_k = 1; end                                                                 % [-] Excitation / Inhibition Sign.
            if nargin < 2, c_k = self.c_absolute_subtraction_DEFAULT; end                                 % [-] Absolute Subtraction Subnetwork Gain.
            
            % Compute the maximum synaptic conductance.
            gs_nk = ( Ia_n - c_k.*s_k.*R_k.*Gm_n )./( c_k.*s_k.*R_k - dEs_nk );                             % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs_nk ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of relative subtraction subnetwork synapses.
        function gs_nk = compute_relative_subtraction_gs( self, c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag )
            
            % Define the default input arguments.
            if nargin < 9, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 8, Ia_n = self.Ia_relative_subtraction_DEFAULT; end                             % [A] Applied Current.
            if nargin < 7, dEs_nk = self.dEs_max_DEFAULT; end                                           % [V] Synaptic Reversal Potential.
            if nargin < 6, Gm_n = self.Gm_DEFAULT; end                                                  % [S] Membrane Conductance.
            if nargin < 5, R_k = self.R_DEFAULT; end                                                    % [V] Activation Domain.
            if nargin < 4, s_k = 1; end                                                                 % [-] Excitation / Inhibition Sign.
            if nargin < 2, c_k = self.c_relative_subtraction_DEFAULT; end                            	% [-] Relative Subtraction Subnetwork Gain.
            
            % Compute the maximum synaptic conductance.
            % gs_nk = ( npm_k.*Ia_n - c.*s_k.*Gm_n.*R_n )./( c.*s_k.*R_n - npm_k.*dEs_nk );	% [S] Maximum Synaptic Conductance.
            gs_nk = ( Ia_n - s_k*c_k*R_k*Gm_n )/( s_k*c_k*R_k - dEs_nk );

            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs_nk ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        

        % Implement a function to compute the maximum synaptic conductance of subtraction subnetwork synapses.
        function gs_nk = compute_subtraction_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c_k = parameters{ 1 };
                s_k = parameters{ 2 };
                R_k = parameters{ 3 };
                Gm_n = parameters{ 4 };
                dEs_nk = parameters{ 5 };
                Ia_n = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs_nk = self.compute_absolute_subtraction_gs( c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c_k = parameters{ 1 };
                s_k = parameters{ 2 };
                R_k = parameters{ 3 };
                Gm_n = parameters{ 4 };
                dEs_nk = parameters{ 5 };
                Ia_n = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs_nk = self.compute_relative_subtraction_gs( c_k, s_k, R_k, Gm_n, dEs_nk, Ia_n, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance of absolute inversion subnetwork synapses.
        function gs21 = compute_absolute_inversion_gs21( self, delta, Gm2, dEs21, Ia2, validation_flag )
            
            % Define the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia2 = self.Ia2_absolute_inversion_DEFAULT; end                               % [A] Applied Current.
            if nargin < 4, dEs21 = self.dEs_absolute_inversion_DEFAULT; end                         	% [V] Synaptic Reversal Potential.
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs21 = ( delta*Gm2 - Ia2 )/( dEs21 - delta );                                             	% [S] Maximum Synaptic Conductance.

            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs21 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of relative inversion subnetwork synapses.
        function gs21 = compute_relative_inversion_gs21( self, delta, Gm2, dEs21, Ia2, validation_flag )
            
            % Define the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia2 = self.Ia2_absolute_inversion_DEFAULT; end                               % [A] Applied Current.
            if nargin < 4, dEs21 = self.dEs_absolute_inversion_DEFAULT; end                         	% [V] Synaptic Reversal Potential.
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs21 = ( delta*Gm2 - Ia2 )/( dEs21 - delta );                                             	% [S] Maximum Synaptic Conductance.

            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs21 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of inversion subnetwork synapses.
        function gs21 = compute_inversion_gs21( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm2 = parameters{ 2 };
                dEs21 = parameters{ 3 };
                Ia2 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs21 = self.compute_absolute_inversion_gs21( delta, Gm2, dEs21, Ia2, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm2 = parameters{ 2 };
                dEs21 = parameters{ 3 };
                Ia2 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs21 = self.compute_relative_inversion_gs21( delta, Gm2, dEs21, Ia2, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to compute the maximum synaptic conductance of reduced absolute inversion subnetwork synapses.
        function gs21 = compute_reduced_absolute_inversion_gs21( self, delta, Gm2, dEs21, Ia2, validation_flag )
        
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia2 = self.Ia_DEFAULT; end
            if nargin < 4, dEs21 = self.dEs_DEFAULT; end
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_reduced_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs21 = ( Ia2 - delta*Gm2 )/( delta - dEs21 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs21 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced relative inversion subnetwork synapses.
        function gs21 = compute_reduced_relative_inversion_gs21( self, delta, Gm2, dEs21, Ia2, validation_flag )
        
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia2 = self.Ia_DEFAULT; end
            if nargin < 4, dEs21 = self.dEs_DEFAULT; end
            if nargin < 3, Gm2 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_reduced_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs21 = ( Ia2 - delta*Gm2 )/( delta - dEs21 );   
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs21 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced inversion subnetwork synapses.
        function gs21 = compute_reduced_inversion_gs21( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm2 = parameters{ 2 };
                dEs21 = parameters{ 3 };
                Ia2 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs21 = self.compute_reduced_absolute_inversion_gs21( delta, Gm2, dEs21, Ia2, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm2 = parameters{ 2 };
                dEs21 = parameters{ 3 };
                Ia2 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs21 = self.compute_reduced_relative_inversion_gs21( delta, Gm2, dEs21, Ia2, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to compute the maximum synaptic conductance of numerator absolute division subnetwork synapses.
        function gs31 = compute_absolute_division_gs31( self, R3, Gm3, dEs31, Ia3, validation_flag )
            
            % Set the default input arugments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia_DEFAULT; end
            if nargin < 4, dEs31 = self.dEs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs31 = ( Ia3 - R3*Gm3 )/( R3 - dEs31 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of numerator relative division subnetwork synapses.
        function gs31 = compute_relative_division_gs31( self, R3, Gm3, dEs31, Ia3, validation_flag )
            
            % Set the default input arugments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia_DEFAULT; end
            if nargin < 4, dEs31 = self.dEs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs31 = ( Ia3 - R3*Gm3 )/( R3 - dEs31 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of numerator reduced inversion subnetwork synapses.
        function gs31 = compute_division_gs31( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                R3 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs31 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs31 = self.compute_absolute_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                R3 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs31 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs31 = self.compute_relative_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division Subnetwork Functions (Synapse 32) ----------
        
        % Implement a function to compute the maximum synaptic conductance of denominator absolute division subnetwork synapses.
        function gs32 = compute_absolute_division_gs32( self, delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, gs31 = self.gs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = ( ( dEs31 - delta )*gs31 + ( Ia3 - delta*Gm3 ) )/( delta - dEs32 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of denominator relative division subnetwork synapses.
        function gs32 = compute_relative_division_gs32( self, delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag )
            
            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, gs31 = self.gs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = ( ( dEs31 - delta )*gs31 + ( Ia3 - delta*Gm3 ) )/( delta - dEs32 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of denominator reduced inversion subnetwork synapses.
        function gs32 = compute_division_gs32( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm3 = parameters{ 2 };
                gs31 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs32 = self.compute_absolute_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm3 = parameters{ 2 };
                gs31 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs32 = self.compute_relative_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division Subnetwork Functions (Combined) ----------

        % Implement a function to compute the maximum synaptic conductance of combined absolute division subnetwork synapses.
        function [ gs31, gs32 ] = compute_absolute_division_gs( self, delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, Gm3 = self.Gm_DEFAULT; end
            if nargin < 3, R3 = self.R_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end
            
            % Compute the maximum synaptic conductance for synapse 31.
            gs31 = self.compute_absolute_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.
            gs32 = self.compute_absolute_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of combined relative division subnetwork synapses.
        function [ gs31, gs32 ] = compute_relative_division_gs( self, delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, Gm3 = self.Gm_DEFAULT; end
            if nargin < 3, R3 = self.R_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end    
            
            % Compute the maximum synaptic conductance for synapse 31.            
            gs31 = self.compute_relative_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.            
            gs32 = self.compute_relative_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced division subnetwork synapses.
        function [ gs31, gs32 ] = compute_division_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta = parameters{ 1 };
                R3 = parameters{ 2 };
                Gm3 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                [ gs31, gs32 ] = self.compute_absolute_division_gs( delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta = parameters{ 1 };
                R3 = parameters{ 2 };
                Gm3 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                [ gs31, gs32 ] = self.compute_relative_division_gs( delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Synapse 31) ----------

        % Implement a function to compute the maximum synaptic conductance of numerator reduced absolute division subnetwork synapses.
        function gs31 = compute_reduced_absolute_division_gs31( self, R3, Gm3, dEs31, Ia3, validation_flag )
            
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia_DEFAULT; end
            if nargin < 4, dEs31 = self.dEs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs31 = ( Ia3 - R3*Gm3 )/( R3 - dEs31 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of numerator reduced relative division subnetwork synapses.
        function gs31 = compute_reduced_relative_division_gs31( self, R3, Gm3, dEs31, Ia3, validation_flag )
            
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia_DEFAULT; end
            if nargin < 4, dEs31 = self.dEs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, R3 = self.R_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs31 = ( Ia3 - R3*Gm3 )/( R3 - dEs31 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of numerator reduced division subnetwork synapses.
        function gs31 = compute_reduced_division_gs31( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                R3 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs31 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs31 = self.compute_reduced_absolute_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                R3 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs31 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs31 = self.compute_reduced_relative_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Synapse 32) ----------        
        
        % Implement a function to compute the maximum synaptic conductance of denominator reduced relative division subnetwork synapses.
        function gs32 = compute_reduced_absolute_division_gs32( self, delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, gs31 = self.gs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = ( ( dEs31 - delta )*gs31 + ( Ia3 - delta*Gm3 ) )/( delta - dEs32 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of denominator reduced relative division subnetwork synapses.
        function gs32 = compute_reduced_relative_division_gs32( self, delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, gs31 = self.gs_DEFAULT; end
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = ( ( dEs31 - delta )*gs31 + ( Ia3 - delta*Gm3 ) )/( delta - dEs32 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of numerator reduced division subnetwork synapses.
        function gs32 = compute_reduced_division_gs32( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm3 = parameters{ 2 };
                gs31 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs32 = self.compute_reduced_absolute_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta = parameters{ 1 };
                Gm3 = parameters{ 2 };
                gs31 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs32 = self.compute_reduced_relative_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions (Combined) ----------

        % Implement a function to compute the maximum synaptic conductance of combined reduced absolute division subnetwork synapses.
        function [ gs31, gs32 ] = compute_reduced_absolute_division_gs( self, delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, Gm3 = self.Gm_DEFAULT; end
            if nargin < 3, R3 = self.R_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end
            
            % Compute the maximum synaptic conductance for synapse 31.
            gs31 = self.compute_reduced_absolute_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.
            gs32 = self.compute_reduced_absolute_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of combined reduced relative division subnetwork synapses.
        function [ gs31, gs32 ] = compute_reduced_relative_division_gs( self, delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, Ia3 = self.Ia_DEFAULT; end
            if nargin < 6, dEs32 = self.dEs_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, Gm3 = self.Gm_DEFAULT; end
            if nargin < 3, R3 = self.R_DEFAULT; end
            if nargin < 2, delta = self.delta_absolute_division_DEFAULT; end
            
            % Compute the maximum synaptic conductance for synapse 31.            
            gs31 = self.compute_reduced_relative_division_gs31( R3, Gm3, dEs31, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.            
            gs32 = self.compute_reduced_relative_division_gs32( delta, Gm3, gs31, dEs31, dEs32, Ia3, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of combined reduced division subnetwork synapses.
        function [ gs31, gs32 ] = compute_reduced_division_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta = parameters{ 1 };
                R3 = parameters{ 2 };
                Gm3 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                [ gs31, gs32 ] = self.compute_reduced_absolute_division_gs( delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta = parameters{ 1 };
                R3 = parameters{ 2 };
                Gm3 = parameters{ 3 };
                dEs31 = parameters{ 4 };
                dEs32 = parameters{ 5 };
                Ia3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                [ gs31, gs32 ] = self.compute_reduced_relative_division_gs( delta, R3, Gm3, dEs31, dEs32, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
                
        
        % ---------- Division After Inversion Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to compute the maximum synaptic conductance of numerator absolute division after inversion subnetwork synapse 31.
        function gs31 = compute_absolute_dai_gs31( self, c1, c3, delta1, delta2, R1, R2, validation_flag )
            
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, R2 = self.R_DEFAULT; end
            if nargin < 6, R1 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_relative_dai_DEFAULT; end
            if nargin < 3, c3 = self.c3_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs31 = ( c1*c3 )/( ( c1*R1*delta1 + c3*R2*delta2 - c3*delta1*delta2 )*R2 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of numerator relative division subnetwork synapse 31.
        function gs31 = compute_relative_dai_gs31( self, c1, c3, delta1, delta2, R2, dEs31, validation_flag )
            
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, R2 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_absolute_inversion_DEFAULT; end
            if nargin < 3, c3 = self.c3_absolute_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_absolute_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs31 = ( ( c3^2 )*delta1*delta2 + ( c1 - c3 )*R2*c3*delta2 )/( -c3*delta1*delta2 + c3*dEs31*delta1 + ( c3 - c1 )*R2*delta2 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of numerator reduced inversion subnetwork synapse 31.
        function gs31 = compute_dai_gs31( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c1 = parameters{ 1 };
                c3 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R1 = parameters{ 5 };
                R2 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs31 = self.compute_absolute_dai_gs31( c1, c3, delta1, delta2, R1, R2, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c1 = parameters{ 1 };
                c3 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R2 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs31 = self.compute_relative_dai_gs31( c1, c3, delta1, delta2, R2, dEs31, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
    
        % ---------- Division After Inversion Subnetwork Functions (Synapse 32) ----------

        % Implement a function to compute the maximum synaptic conductance of absolute division after inversion subnetwork synapse 32.
        function gs32 = compute_absolute_dai_gs32( self, c1, c3, delta2, R1, R2, dEs31, validation_flag )
            
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, R2 = self.R_DEFAULT; end
            if nargin < 5, R1 = self.R_DEFAULT; end
            if nargin < 4, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 3, c3 = self.c3_absolute_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_absolute_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = ( ( delta2*c3 - R1*c1 )*dEs31*c3 )/( ( R1*c1 - dEs31*c3 )*R1*R2*delta2 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of relative division after inversion subnetwork synapse 32.
        function gs32 = compute_relative_dai_gs32( self, c1, c3, delta1, delta2, R2, dEs31, validation_flag )
            
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, R2 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_relative_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_relative_inversion_DEFAULT; end
            if nargin < 3, c3 = self.c3_relative_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_relative_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = ( ( c1 - c3 )*c3*R2*dEs31 )/( -c3*delta1*delta2 + c3*dEs31*delta1 + ( c3 - c1 )*R2*delta2 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of division after inversion subnetwork synapse 32.
        function gs32 = compute_dai_gs32( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c1 = parameters{ 1 };
                c3 = parameters{ 2 };
                delta2 = parameters{ 3 };
                R1 = parameters{ 4 };
                R2 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.                
                gs32 = self.compute_absolute_dai_gs32( c1, c3, delta2, R1, R2, dEs31, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c1 = parameters{ 1 };
                c3 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R2 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.            
                gs32 = self.compute_relative_dai_gs32( c1, c3, delta1, delta2, R2, dEs31, validation_flag );
                
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions (Combined) ----------

        % Implement a function to compute the maximum synaptic conductance of combined absolute division subnetwork synapses.
        function [ gs31, gs32 ] = compute_absolute_dai_gs( self, c1, c3, delta1, delta2, R1, R2, dEs31, validation_flag )

            % Set the default input arguments.
            if nargin < 9, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 8, dEs31 = self.dEs_DEFAULT; end
            if nargin < 7, R2 = self.R_DEFAULT; end
            if nargin < 6, R1 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_absolute_inversion_DEFAULT; end
            if nargin < 3, c3 = self.c3_absolute_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_absolute_dai_DEFAULT; end
                        
            % Compute the maximum synaptic conductance for synapse 31.
            gs31 = self.compute_absolute_dai_gs31( c1, c3, delta1, delta2, R1, R2, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.
            gs32 = self.compute_absolute_dai_gs32( c1, c3, delta2, R1, R2, dEs31, validation_flag );
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of combined relative division subnetwork synapses.
        function [ gs31, gs32 ] = compute_relative_dai_gs( self, c1, c3, delta1, delta2, R2, dEs31, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, R2 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_relative_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_relative_inversion_DEFAULT; end
            if nargin < 3, c3 = self.c3_relative_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_relative_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance for synapse 31.            
            gs31 = self.compute_relative_dai_gs31( c1, c3, delta1, delta2, R2, dEs31, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.            
            gs32 = self.compute_relative_dai_gs32( c1, c3, delta1, delta2, R2, dEs31, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced division subnetwork synapses.
        function [ gs31, gs32 ] = compute_dai_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c1 = parameters{ 1 };
                c3 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R1 = parameters{ 5 };
                R2 = parameters{ 6 };
                dEs31 = parameters{ 7 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                [ gs31, gs32 ] = self.compute_absolute_dai_gs( c1, c3, delta1, delta2, R1, R2, dEs31, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c1 = parameters{ 1 };
                c3 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R2 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                [ gs31, gs32 ] = self.compute_relative_dai_gs( c1, c3, delta1, delta2, R2, dEs31, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Synapse 31) ----------
        
        % Implement a function to compute the maximum synaptic conductance of reduced absolute division after inversion subnetwork synapse 31.
        function gs31 = compute_reduced_absolute_dai_gs31( self, delta1, delta2, R2, R3, Gm3, dEs31, validation_flag )
            
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, Gm3 = self.Gm_DEFAULT; end
            if nargin < 5, R3 = self.R_DEFAULT; end
            if nargin < 4, R2 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs31 = ( ( delta1 - R2 )*delta2*R3*Gm3 )/( ( R2 - delta1 )*delta2*R3 + ( delta1*R3 - delta2*R2 )*dEs31 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of reduced relative division subnetwork synapse 31.
        function gs31 = compute_reduced_relative_dai_gs31( self, delta1, delta2, R2, R3, dEs31, validation_flag )
            
            % Set the default input arugments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, dEs31 = self.dEs_DEFAULT; end
            if nargin < 4, R2 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_absolute_inversion_DEFAULT; end
                        
            % Compute the maximum synaptic conductance.
            gs31 = ( ( delta1 - R2 )*delta2*R3*Gm3 )/( ( R2 - delta1 )*delta2*R3 + ( delta1*R3 - delta2*R2 )*dEs31 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs31 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced division after inversion subnetwork synapse 31.
        function gs31 = compute_reduced_dai_gs31( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R2 = parameters{ 3 };
                R3 = parameters{ 4 };
                Gm3 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs31 = self.compute_reduced_absolute_dai_gs31( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R2 = parameters{ 3 };
                R3 = parameters{ 4 };
                dEs31 = parameters{ 5 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs31 = self.compute_reduced_relative_dai_gs31( delta1, delta2, R2, R3, dEs31, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Synapse 32) ----------

        % Implement a function to compute the maximum synaptic conductance of reduced absolute division after inversion subnetwork synapse 32.
        function gs32 = compute_reduced_absolute_dai_gs32( self, delta1, delta2, R2, R3, Gm3, dEs31, validation_flag )
            
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, Gm3 = self.Gm_DEFAULT; end
            if nargin < 5, R3 = self.R_DEFAULT; end
            if nargin < 4, R2 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = ( ( delta2 - R3 )*R2*Gm3*dEs31 )/( ( R2 - delta1 )*delta2*R3 + ( delta1*R3- delta2*R2 )*dEs31 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of reduced relative division after inversion subnetwork synapse 32.
        function gs32 = compute_reduced_relative_dai_gs32( self, delta1, delta2, R2, R3, Gm3, dEs31, validation_flag )
            
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, Gm3 = self.Gm_DEFAULT; end
            if nargin < 5, R3 = self.R_DEFAULT; end
            if nargin < 4, R2 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
                        
            % Compute the maximum synaptic conductance.
            gs32 = ( ( delta2 - R3 )*R3*Gm3*dEs31 )/( ( R2 - delta1 )*delta2*R3 + ( delta1*R3 - delta2*R2 )*dEs31 );                                                         % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs32 ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced division after inversion subnetwork synapse 32.
        function gs32 = compute_reduced_dai_gs32( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R2 = parameters{ 3 };
                R3 = parameters{ 4 };
                Gm3 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.                
                gs32 = self.compute_reduced_absolute_dai_gs32( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R2 = parameters{ 3 };
                R3 = parameters{ 4 };
                Gm3 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.            
                gs32 = self.compute_reduced_relative_dai_gs32( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
                
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions (Combined) ----------

        % Implement a function to compute the maximum synaptic conductance of combined reduced absolute division subnetwork synapses.
        function [ gs31, gs32 ] = compute_reduced_absolute_dai_gs( self, delta1, delta2, R2, R3, Gm3, dEs31, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, Gm3 = self.Gm_DEFAULT; end
            if nargin < 5, R3 = self.R_DEFAULT; end
            if nargin < 4, R2 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance for synapse 31.            
            gs31 = self.compute_reduced_absolute_dai_gs31( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.            
            gs32 = self.compute_reduced_absolute_dai_gs32( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of combined reduced relative division subnetwork synapses.
        function [ gs31, gs32 ] = compute_reduced_relative_dai_gs( self, delta1, delta2, R2, R3, Gm3, dEs31, validation_flag )

            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs31 = self.dEs_DEFAULT; end
            if nargin < 6, Gm3 = self.Gm_DEFAULT; end
            if nargin < 5, R3 = self.R_DEFAULT; end
            if nargin < 4, R2 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_relative_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_relative_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance for synapse 31.                        
            gs31 = self.compute_reduced_relative_dai_gs31( delta1, delta2, R2, R3, dEs31, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.                        
            gs32 = self.compute_reduced_relative_dai_gs32( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced division subnetwork synapses.
        function [ gs31, gs32 ] = compute_reduced_dai_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R2 = parameters{ 3 };
                R3 = parameters{ 4 };
                Gm3 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                [ gs31, gs32 ] = self.compute_reduced_absolute_dai_gs( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R2 = parameters{ 3 };
                R3 = parameters{ 4 };
                Gm3 = parameters{ 5 };
                dEs31 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                [ gs31, gs32 ] = self.compute_reduced_relative_dai_gs( delta1, delta2, R2, R3, Gm3, dEs31, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 41) ----------
        
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 41.
        function gs41 = compute_absolute_multiplication_gs41( self, c4, c6, delta1, delta2, R1, R3, validation_flag )
        
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, R3 = self.R_DEFAULT; end
            if nargin < 6, R1 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_absolute_inversion_DEFAULT; end
            if nargin < 3, c6 = self.c3_absolute_dai_DEFAULT; end
            if nargin < 2, c4 = self.c1_absolute_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs41 = self.compute_absolute_dai_gs31( c4, c6, delta1, delta2, R1, R3, validation_flag );
            
        end
        
            
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 41.
        function gs41 = compute_relative_multiplication_gs41( self, c4, c6, delta1, delta2, R3, dEs41, validation_flag )
           
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs41 = self.dEs_DEFAULT; end
            if nargin < 6, R3 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_relative_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_relative_inversion_DEFAULT; end
            if nargin < 3, c6 = self.c3_relative_dai_DEFAULT; end
            if nargin < 2, c4 = self.c1_relative_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs41 = self.compute_relative_dai_gs31( c4, c6, delta1, delta2, R3, dEs41, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 41.
        function gs41 = compute_multiplication_gs41( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c4 = parameters{ 1 };
                c6 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R1 = parameters{ 5 };
                R3 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs41 = self.compute_absolute_multiplication_gs41( c4, c6, delta1, delta2, R1, R3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c4 = parameters{ 1 };
                c6 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R3 = parameters{ 5 };
                dEs41 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs41 = self.compute_relative_multiplication_gs41( c4, c6, delta1, delta2, R3, dEs41, validation_flag );
                            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 32) ----------

        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 32.
        function gs32 = compute_absolute_multiplication_gs32( self, delta1, Gm3, dEs32, Ia3, validation_flag )
        
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia2_absolute_inversion_DEFAULT; end                               % [A] Applied Current.
            if nargin < 4, dEs32 = self.dEs_absolute_inversion_DEFAULT; end                         	% [V] Synaptic Reversal Potential.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta1 = self.delta_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = self.compute_absolute_inversion_gs21( delta1, Gm3, dEs32, Ia3, validation_flag );
            
        end
        
            
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 32.
        function gs32 = compute_relative_multiplication_gs32( self, delta1, Gm3, dEs32, Ia3, validation_flag )
           
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia2_absolute_inversion_DEFAULT; end                               % [A] Applied Current.
            if nargin < 4, dEs32 = self.dEs_absolute_inversion_DEFAULT; end                         	% [V] Synaptic Reversal Potential.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta1 = self.delta_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = self.compute_relative_inversion_gs21( delta1, Gm3, dEs32, Ia3, validation_flag );            
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 32.
        function gs32 = compute_multiplication_gs32( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs32 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs32 = self.compute_absolute_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs32 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs32 = self.compute_relative_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Synapse 43) ----------

        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 43.
        function gs43 = compute_absolute_multiplication_gs43( self, c4, c6, delta2, R1, R3, dEs41, validation_flag )
        
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs41 = self.dEs_DEFAULT; end
            if nargin < 6, R3 = self.R_DEFAULT; end
            if nargin < 5, R1 = self.R_DEFAULT; end
            if nargin < 4, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 3, c6 = self.c3_absolute_dai_DEFAULT; end
            if nargin < 2, c4 = self.c1_absolute_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.            
            gs43 = self.compute_absolute_dai_gs32( c4, c6, delta2, R1, R3, dEs41, validation_flag );
            
        end
        
            
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 43.
        function gs43 = compute_relative_multiplication_gs43( self, c4, c6, delta1, delta2, R3, dEs41, validation_flag )
           
            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs41 = self.dEs_DEFAULT; end
            if nargin < 6, R3 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_relative_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_relative_inversion_DEFAULT; end
            if nargin < 3, c6 = self.c3_relative_dai_DEFAULT; end
            if nargin < 2, c4 = self.c1_relative_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance.            
            gs43 = compute_relative_dai_gs32( c4, c6, delta1, delta2, R3, dEs41, validation_flag );
                        
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 43.
        function gs43 = compute_multiplication_gs43( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c4 = parameters{ 1 };
                c6 = parameters{ 2 };
                delta2 = parameters{ 3 };
                R1 = parameters{ 4 };
                R3 = parameters{ 5 };
                dEs41 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.                                
                gs43 = self.compute_absolute_multiplication_gs43( c4, c6, delta2, R1, R3, dEs41, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c4 = parameters{ 1 };
                c6 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R3 = parameters{ 5 };
                dEs41 = parameters{ 6 };
                                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs43 = self.compute_relative_multiplication_gs43( c4, c6, delta1, delta2, R3, dEs41, validation_flag );
                
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions (Combined) ----------

        % Implement a function to compute the maximum synaptic conductance of combined absolute multiplication subnetwork synapses.
        function [ gs41, gs32, gs43 ] = compute_absolute_multiplication_gs( self, c4, c6, delta1, delta2, R1, R3, Gm3, dEs41, dEs32, Ia3, validation_flag )
            
            % Set the default input arguments.
            if nargin < 12, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 11, Ia3 = self.Ia_DEFAULT; end
            if nargin < 10, dEs32 = self.dEs_DEFAULT; end
            if nargin < 9, dEs41 = self.dEs_DEFAULT; end
            if nargin < 8, Gm3 = self.Gm_DEFAULT; end
            if nargin < 7, R3 = self.R_DEFAULT; end
            if nargin < 6, R1 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_absolute_inversion_DEFAULT; end
            if nargin < 3, c6 = self.c3_absolute_dai_DEFAULT; end
            if nargin < 2, c4 = self.c1_absolute_dai_DEFAULT; end
                        
            % Compute the maximum synaptic conductance for synapse 41.                        
            gs41 = self.compute_absolute_multiplication_gs41( c4, c6, delta1, delta2, R1, R3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.            
            gs32 = self.compute_absolute_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 43.
            gs43 = self.compute_absolute_multiplication_gs43( c4, c6, delta2, R1, R3, dEs41, validation_flag );
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of combined relative multiplication subnetwork synapses.
        function [ gs41, gs32, gs43 ] = compute_relative_multiplication_gs( self, c4, c6, delta1, delta2, R3, Gm3, dEs41, dEs32, Ia3, validation_flag )
            
            % Set the default input arguments.
            if nargin < 11, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 10, Ia3 = self.Ia_DEFAULT; end
            if nargin < 9, dEs32 = self.dEs_DEFAULT; end
            if nargin < 8, dEs41 = self.dEs_DEFAULT; end
            if nargin < 7, Gm3 = self.Gm_DEFAULT; end
            if nargin < 6, R3 = self.R_DEFAULT; end
            if nargin < 5, delta2 = self.delta_relative_dai_DEFAULT; end
            if nargin < 4, delta1 = self.delta_relative_inversion_DEFAULT; end
            if nargin < 3, c6 = self.c3_relative_dai_DEFAULT; end
            if nargin < 2, c4 = self.c1_relative_dai_DEFAULT; end
            
            % Compute the maximum synaptic conductance for synapse 41.            
            gs41 = self.compute_relative_multiplication_gs41( c4, c6, delta1, delta2, R3, dEs41, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.            
            gs32 = self.compute_relative_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 43.
            gs43 = self.compute_relative_multiplication_gs43( c4, c6, delta1, delta2, R3, dEs41, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced multiplication subnetwork synapses.
        function [ gs41, gs32, gs43 ] = compute_multiplication_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                c4 = parameters{ 1 };
                c6 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R1 = parameters{ 5 };
                R3 = parameters{ 6 };
                Gm3 = parameters{ 7 };
                dEs41 = parameters{ 8 };
                dEs32 = parameters{ 9 };
                Ia3 = parameters{ 10 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.                
                [ gs41, gs32, gs43 ] = self.compute_absolute_multiplication_gs( c4, c6, delta1, delta2, R1, R3, Gm3, dEs41, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                c4 = parameters{ 1 };
                c6 = parameters{ 2 };
                delta1 = parameters{ 3 };
                delta2 = parameters{ 4 };
                R3 = parameters{ 5 };
                Gm3 = parameters{ 6 };
                dEs41 = parameters{ 7 };
                dEs32 = parameters{ 8 };
                Ia3 = parameters{ 9 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.                
                [ gs41, gs32, gs43 ] = self.compute_relative_multiplication_gs( c4, c6, delta1, delta2, R3, Gm3, dEs41, dEs32, Ia3, validation_flag );
                
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 41) ----------
        
        % Implement a function to compute the maximum synaptic conductance of reduced absolute multiplication subnetwork synapse 41.
        function gs41 = compute_reduced_absolute_multiplication_gs41( self, delta1, delta2, R3, R4, Gm4, dEs41, validation_flag )
        
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs41 = self.dEs_DEFAULT; end
            if nargin < 6, Gm4 = self.Gm_DEFAULT; end
            if nargin < 5, R4 = self.R_DEFAULT; end
            if nargin < 4, R3 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.                        
            gs41 = self.compute_reduced_absolute_dai_gs31( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
            
        end
        
            
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 41.
        function gs41 = compute_reduced_relative_multiplication_gs41( self, delta1, delta2, R3, R4, dEs41, validation_flag )
           
            % Set the default input arugments.
            if nargin < 7, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 6, dEs41 = self.dEs_DEFAULT; end
            if nargin < 5, R4 = self.R_DEFAULT; end
            if nargin < 4, R3 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_relative_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_relative_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs41 = self.compute_reduced_relative_dai_gs31( delta1, delta2, R3, R4, dEs41, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 41.
        function gs41 = compute_reduced_multiplication_gs41( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R3 = parameters{ 3 };
                R4 = parameters{ 4 };
                Gm4 = parameters{ 5 };
                dEs41 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.                
                gs41 = self.compute_reduced_absolute_multiplication_gs41( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R3 = parameters{ 3 };
                R4 = parameters{ 4 };
                dEs41 = parameters{ 5 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs41 = self.compute_reduced_relative_multiplication_gs41( delta1, delta2, R3, R4, dEs41, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 32) ----------

        % Implement a function to compute the maximum synaptic conductance of reduced absolute multiplication subnetwork synapse 32.
        function gs32 = compute_reduced_absolute_multiplication_gs32( self, delta1, Gm3, dEs32, Ia3, validation_flag )
        
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia2_absolute_inversion_DEFAULT; end                               % [A] Applied Current.
            if nargin < 4, dEs32 = self.dEs_absolute_inversion_DEFAULT; end                         	% [V] Synaptic Reversal Potential.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta1 = self.delta_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = self.compute_reduced_absolute_inversion_gs21( delta1, Gm3, dEs32, Ia3, validation_flag );
            
        end
        
            
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 32.
        function gs32 = compute_reduced_relative_multiplication_gs32( self, delta1, Gm3, dEs32, Ia3, validation_flag )
           
            % Set the default input arguments.
            if nargin < 6, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 5, Ia3 = self.Ia2_absolute_inversion_DEFAULT; end                               % [A] Applied Current.
            if nargin < 4, dEs32 = self.dEs_absolute_inversion_DEFAULT; end                         	% [V] Synaptic Reversal Potential.
            if nargin < 3, Gm3 = self.Gm_DEFAULT; end
            if nargin < 2, delta1 = self.delta_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.
            gs32 = self.compute_reduced_relative_inversion_gs21( delta1, Gm3, dEs32, Ia3, validation_flag );            
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 32.
        function gs32 = compute_reduced_multiplication_gs32( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs32 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.
                gs32 = self.compute_reduced_absolute_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                Gm3 = parameters{ 2 };
                dEs32 = parameters{ 3 };
                Ia3 = parameters{ 4 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs32 = self.compute_reduced_relative_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Synapse 43) ----------

        % Implement a function to compute the maximum synaptic conductance of reduced absolute multiplication subnetwork synapse 43.
        function gs43 = compute_reduced_absolute_multiplication_gs43( self, delta1, delta2, R3, R4, Gm4, dEs41, validation_flag )
        
            % Set the default input arugments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs41 = self.dEs_DEFAULT; end
            if nargin < 6, Gm4 = self.Gm_DEFAULT; end
            if nargin < 5, R4 = self.R_DEFAULT; end
            if nargin < 4, R3 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.            
            gs43 = self.compute_reduced_absolute_dai_gs32( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
                        
        end
        
            
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 43.
        function gs43 = compute_reduced_relative_multiplication_gs43( self, delta1, delta2, R3, R4, Gm4, dEs41, validation_flag )
           
            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 7, dEs41 = self.dEs_DEFAULT; end
            if nargin < 6, Gm4 = self.Gm_DEFAULT; end
            if nargin < 5, R4 = self.R_DEFAULT; end
            if nargin < 4, R3 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_reduced_relative_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_reduced_relative_inversion_DEFAULT; end
            
            % Compute the maximum synaptic conductance.                        
            gs43 = self.compute_reduced_relative_dai_gs32( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of absolute multiplication subnetwork synapse 43.
        function gs43 = compute_reduced_multiplication_gs43( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R3 = parameters{ 3 };
                R4 = parameters{ 4 };
                Gm4 = parameters{ 5 };
                dEs41 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.                
                gs43 = self.compute_reduced_absolute_multiplication_gs43( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R3 = parameters{ 3 };
                R4 = parameters{ 4 };
                Gm4 = parameters{ 5 };
                dEs41 = parameters{ 6 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.
                gs43 = self.compute_reduced_relative_multiplication_gs43( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
                
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions (Combined) ----------
        
        % Implement a function to compute the maximum synaptic conductance of combined reduced absolute multiplication subnetwork synapses.
        function [ gs41, gs32, gs43 ] = compute_reduced_absolute_multiplication_gs( self, delta1, delta2, R3, R4, Gm3, Gm4, dEs41, dEs32, Ia3, validation_flag )
            
            % Set the default input arguments.
            if nargin < 11, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 10, Ia3 = self.Ia_DEFAULT; end
            if nargin < 9, dEs32 = self.dEs_DEFAULT; end
            if nargin < 8, dEs41 = self.dEs_DEFAULT; end
            if nargin < 7, Gm4 = self.Gm_DEFAULT; end
            if nargin < 6, Gm3 = self.Gm_DEFAULT; end
            if nargin < 5, R4 = self.R_DEFAULT; end
            if nargin < 4, R3 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_absolute_dai_DEFAULT; end
            if nargin < 2, delta1 = self.delta_absolute_inversion_DEFAULT; end
                        
            % Compute the maximum synaptic conductance for synapse 41.            
            gs41 = self.compute_reduced_absolute_multiplication_gs41( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.
            gs32 = self.compute_reduced_absolute_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 43.
            gs43 = self.compute_reduced_absolute_multiplication_gs43( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
            
        end

        
        % Implement a function to compute the maximum synaptic conductance of combined relative multiplication subnetwork synapses.
        function [ gs41, gs32, gs43 ] = compute_reduced_relative_multiplication_gs( self, delta1, delta2, R3, R4, Gm3, Gm4, dEs41, dEs32, Ia3, validation_flag )
            
            % Set the default input arguments.
            if nargin < 11, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 10, Ia3 = self.Ia_DEFAULT; end
            if nargin < 9, dEs32 = self.dEs_DEFAULT; end
            if nargin < 8, dEs41 = self.dEs_DEFAULT; end
            if nargin < 7, Gm4 = self.Gm_DEFAULT; end
            if nargin < 6, Gm3 = self.Gm_DEFAULT; end
            if nargin < 5, R4 = self.R_DEFAULT; end
            if nargin < 4, R3 = self.R_DEFAULT; end
            if nargin < 3, delta2 = self.delta_division_DEFAULT; end
            if nargin < 2, delta1 = self.delta_inversion_DEFAULT; end
                        
            % Compute the maximum synaptic conductance for synapse 41.            
            gs41 = self.compute_reduced_relative_multiplication_gs41( delta1, delta2, R3, R4, dEs41, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 32.            
            gs32 = self.compute_reduced_relative_multiplication_gs32( delta1, Gm3, dEs32, Ia3, validation_flag );
            
            % Compute the maximum synaptic conductance for synapse 43.
            gs43 = self.compute_reduced_relative_multiplication_gs43( delta1, delta2, R3, R4, Gm4, dEs41, validation_flag );
            
        end
        
        
        % Implement a function to compute the maximum synaptic conductance of reduced multiplication subnetwork synapses.
        function [ gs41, gs32, gs43 ] = compute_reduced_multiplication_gs( self, parameters, encoding_scheme, validation_flag )
        
            % Set the default input arguments.
            if nargin < 4, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the synaptic reversal potential.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R3 = parameters{ 3 };
                R4 = parameters{ 4 };
                Gm3 = parameters{ 5 };
                Gm4 = parameters{ 6 };
                dEs41 = parameters{ 7 };
                dEs32 = parameters{ 8 };
                Ia3 = parameters{ 9 };
                
                % Compute the synaptic reversal potential using an absolute encoding scheme.                
                [ gs41, gs32, gs43 ] = self.compute_reduced_absolute_multiplication_gs( delta1, delta2, R3, R4, Gm3, Gm4, dEs41, dEs32, Ia3, validation_flag );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Unpack the parameters.
                delta1 = parameters{ 1 };
                delta2 = parameters{ 2 };
                R3 = parameters{ 3 };
                R4 = parameters{ 4 };
                Gm3 = parameters{ 5 };
                Gm4 = parameters{ 6 };
                dEs41 = parameters{ 7 };
                dEs32 = parameters{ 8 };
                Ia3 = parameters{ 9 };
                
                % Compute the synaptic reversal potential using a relative encoding scheme.                
                [ gs41, gs32, gs43 ] = self.compute_reduced_relative_multiplication_gs( delta1, delta2, R3, R4, Gm3, Gm4, dEs41, dEs32, Ia3, validation_flag );
                
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Derivation Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductances for a derivative subnetwork.
        function gs = compute_derivation_gs( self, Gm3, R1, dEs13, dEs23, Ia3, k, validation_flag )
            
            %{
            Input(s):
                Gm3             =   [S] Membrane Conductance (Neuron 3).
                R1              =   [V] Maximum Membrane Voltage (Neuron 1).
                dE_syn13        =   [V] Synaptic Reversal Potential (Synapse 13).
                dE_syn23        =   [V] Synaptic Reversal Potential (Synapse 23).
                I_app3          =   [A] Applied Current (Neuron 3).
                k               =   [-] Derivation Subnetwork Gain.
            
            Output(s):
                g_syn_max13     =   [S] Maximum Synaptic Conductance (Synapse 13).
                g_syn_max23     =   [S] Maximum Synaptic Conductnace (Synapse 23).
            %}
            
            % Set the default input arguments.
            if nargin < 8, validation_flag = self.validation_flag; end
            
            % Compute the maximum synaptic conductances for a derivative subnetwork in the same way as for a subtraction subnetwork.
            gs = self.compute_subtraction_gs( Gm3, R1, dEs13, dEs23, Ia3, k, validation_flag );
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
                
        % Implement a function to compute the maximum synaptic conductances for an integration subnetwork.
        function gs = compute_integration_gs( self, Gm, Cm, ki_range, validation_flag )
        
            %{
            Input(s):
                Gm          =   [S] Membrane Conductance.
                Cm          =   [F] Membrane Capacitance.
                ki_range    =   [-] Integration Subnetwork Gain.
            
            Output(s):
                gs          =   [S] Maximum Synaptic Conductance.
            %}
            
            % Set the default input arguments.
            if nargin < 5, validation_flag = self.validation_flag_DEFAULT; end
            if nargin < 4, ki_range = self.c_integration_range_DEFAULT; end
            
            % Compute the integration subnetwork maximum synaptic conductances.
            gs = ( -2*Gm.*Cm.*ki_range )./( Cm.*ki_range - 1 );
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
        %{
        
        % Implement a function to compute the maximum synaptic conductances for a voltage based integration subnetwork.
        function gs12 = compute_vbi_gs( ~, R2, dEs12, Is12 )
                    
            %{
            Input(s):
                R2              =   [V] Maximum Membrane Voltage.
                dEs12        =   [V] Synaptic Reversal Potential (Synapse 12).
                Is12         =   [A] Synaptic Current (Synapse 12).
            
            Output(s):
                gs12     =   [S] Maximum Synaptic Conductance (Synapse 12).
            %}
            
            % Compute the maximum synaptic conductance for a voltage based integration subnetwork.
            gs12 = Is12./( dEs12 - ( R2/2 ) );
            
        end
        
        %}
        
        
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to compute the maximum synaptic conductance for a driven multistate cpg subnetwork.
        function gs = compute_dmcpg_gs( self, dEs, delta_oscillatory, Id_max, validation_flag )
            
            % Define the default input arguments.
            if nargin < 5, validation_flag = self.validation_flag; end
            if nargin < 4, Id_max = self.Id_max_DEFAULT; end                                            % [A] Maximum Drive Current.
            if nargin < 3, delta_oscillatory = self.delta_oscillatory_DEFAULT; end                      % [-] Oscillatory Delta.
            if nargin < 2, dEs = self.dEs_max_DEFAULT; end                                              % [V] Synaptic Reversal Potential.
            
            % Compute the maximum synaptic conductance.
            gs = Id_max./( dEs - delta_oscillatory );                                                   % [S] Maximum Synaptic Conductance.
            
            % Determine whether to validate the synaptic conductance.
            if validation_flag              % If we want to validate the synaptic conductances...
            
                % Ensure that the synaptic conductance is valid.
                assert( self.validate_gs( gs ), 'Invalid synaptic conductance detected.' )
            
            end
            
        end
        
        
    end
end