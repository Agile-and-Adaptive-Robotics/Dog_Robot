classdef neuron_class
    
    % This class contains properties and methods related to neurons.
    
    %% NEURON PROPERTIES
    
    % Define the class properties.
    properties
        
        ID                                                      	% [#] Neuron ID.
        name                                                     	% [-] Neuron Name.
        
        U                                                        	% [V] Membrane Voltage.
        h                                                           % [-] Sodium Channel Deactivation Parameter.
        
        Cm                                                       	% [C] Membrance Capacitance.
        Gm                                                        	% [S] Membrane Conductance.
        Er                                                       	% [V] Membrane Resting Potential.
        R                                                        	% [V] Activation Domain.
        
        Am                                                        	% [V] Sodium Channel Activation Amplitude.
        Sm                                                       	% [V] Sodium Channel Activation Slope.
        dEm                                                        	% [V] Sodium Channel Activation Reversal Potential.
        
        Ah                                                       	% [-] Sodium Channel Deactivation Amplitude.
        Sh                                                        	% [-] Sodium Channel Deactivation Slope.
        dEh                                                     	% [V] Sodium Channel Deactivation Reversal Potential.
        
        dEna                                                       	% [V] Sodium Channel Reversal Potential.
        tauh_max                                                 	% [s] Maximum Sodium Channel Deactivation Time Constant.
        tauh                                                    	% [s] Sodium Channel Deactivation Time Constant.
        Gna                                                         % [S] Sodium Channel Conductance.
        
        minf                                                      	% [-] Steady State Sodium Channel Activation Parameter.
        hinf                                                      	% [-] Steady State Sodium Channel Deactivation Parameter.
        
        Ileak                                                      	% [A] Leak Current.
        Isyn                                                    	% [A] Synaptic Current.
        Ina                                                      	% [A] Sodium Channel Current.
        Itonic                                                   	% [A] Tonic Current.
        Iapp                                                     	% [A] Applied Current.
        Itotal                                                    	% [A] Total Current.
        
        enabled_flag                                                % [-] [T/F] Enable Flag.
        
        neuron_utilities                                            % [-] Neuron Utilities Class.
        
    end
    
    
    % Define private, constant class properties.
    properties ( Access = private, Constant = true )
        
        % ---------- Neuron Properties ----------
        
        % Define the neuron parameters.
        ID_DEFAULT = 0;                                          	% [#] Default Neuron ID.
        name_DEFAULT = '';                                         	% [-] Default Neuron Name.
        U_DEFUALT = 0;                                             	% [V] Default Membrane Voltage.
        h_DEFAULT = NaN;                                          	% [-] Default Sodium Channel Deactivation Parameter.
        Cm_DEFAULT = 5e-9;                                        	% [C] Default Membrane Capacitance.
        Gm_DEFAULT = 1e-6;                                        	% [S] Default Membrane Conductance.
        Er_DEFAULT = -60e-3;                                      	% [V] Default Equilibrium Voltage.
        R_DEFAULT = 20e-3;                                         	% [V] Default Activation Domain.
        Am_DEFAULT = 1;                                          	% [-] Default Sodium Channel Activation Parameter Amplitude.
        Sm_DEFAULT = -50;                                        	% [-] Default Sodium Channel Activation Parameter Slope.
        dEm_DEFAULT = 40e-3;                                     	% [V] Default Sodium Channel Activation Reversal Potential.
        Ah_DEFAULT = 0.5;                                        	% [-] Default Sodium Channel Deactivation Parameter Amplitude.
        Sh_DEFAULT = 50;                                          	% [-] Default Sodium Channel Deactivation Parameter Slope.
        dEh_DEFAULT = 0;                                          	% [V] Default Sodium Channel Deactivation Reversal Potential.
        dEna_DEFAULT = 110e-3;                                    	% [V] Default Sodium Channel Reversal Potential.
        tauh_max_DEFAULT = 0.25;                                   	% [s] Default Maximum Sodium Channel Steady State Time Constant.
        Gna_DEFAULT = 1e-6;                                      	% [S] Default Sodium Channel Conductance.
        Ileak_DEFAULT = 0;                                         	% [A] Default Leak Current.
        Isyn_DEFAULT = 0;                                        	% [A] Default Synaptic Current.
        Ina_DEFAULT = 0;                                        	% [A] Default Sodium Channel Current.
        Itonic_DEFAULT = 0;                                       	% [A] Default Tonic Current.
        Iapp_DEFAULT = 0;                                         	% [A] Default Applied Current.
        Itotal_DEFAULT = 0;                                      	% [A] Default Total Current.
        
        
        % ---------- Transmission Properties ----------
        
        % Define transmission subnetwork gains.
        c_absolute_transmission_DEFAULT = 1.0;                      % [-] Absolute Transmission Gain.
        c_relative_transmission_DEFAULT = 1.0;                      % [-] Relative Transmission Gain.

        
        % ---------- Addition Properties ----------
        
        % Define addition subnetwork gains.
        c_absolute_addition_DEFAULT = 1.0;                          % [-] Absolute Addition Gain.
        c_relative_addition_DEFAULT = 1.0;                          % [-] Relative Addition Gain.
        
        
        % ---------- Subtraction Properties ----------
        
        % Define subtraction subnetwork gains.
        c_absolute_subtraction_DEFAULT = 1.0;                       % [-] Absolute Subtraction Gain.
        c_relative_subtraction_DEFAULT = 1.0;                       % [-] Relative Subtraction Gain.
        
        % Define subtraction subnetwork parameters.
        s_ks_DEFAULT = [ 1, -1 ];                               	% [-] Default Subtraction Input Signature.
        
        % ---------- Inversion Properties ----------

        % Define absolute inversion subnetwork gains.
        c1_absolute_inversion_DEFAULT = 1.0;                        % [-] Absolute Inversion Gain 1.
        c2_absolute_inversion_DEFAULT = 1.0;                        % [-] Absolute Inversion Gain 2.
        c3_absolute_inversion_DEFAULT = 1.0;                        % [-] Absolute Inversion Gain 3.
        
        % Define relative inversion subnetwork gains.
        c1_relative_inversion_DEFAULT = 1.0;                        % [-] Relative Inversion Gain 1.
        c2_relative_inversion_DEFAULT = 1.0;                        % [-] Relative Inversion Gain 2.
        c3_relative_inversion_DEFAULT = 1.0;                        % [-] Relative Inversion Gain 3.
        
        % Define inversion subnetwork offsets.
        delta_absolute_inversion_DEFAULT = 1e-3;                    % [V] Absolute Inversion Offset.
        delta_relative_inversion_DEFAULT = 1e-3;                    % [V] Relative Inversion Offset.

        
        % ---------- Reduced Inversion Properties ----------
        
        % Define the reduced absolute inversion subnetwork gain.
        c1_reduced_absolute_inversion_DEFAULT = 1.0;                % [-] Reduced Absolute Inversion Gain 1.
        c2_reduced_absolute_inversion_DEFAULT = 1.0;                % [-] Reduced Absolute Inversion Gain 2.
        
        % Define the reduced relative inversion subnetwork gain.
        c1_reduced_relative_inversion_DEFAULT = 1.0;                % [-] Reduced Relative Inversion Gain 1.
        c2_reduced_relative_inversion_DEFAULT = 1.0;                % [-] Reduced Relative Inversion Gain 2.
        
        % Define reduced inversion subnetwork offsets.
        delta_reduced_absolute_inversion_DEFAULT = 1e-3;            % [V] Reduced Absolute Inversion Offset.
        delta_reduced_relative_inversion_DEFAULT = 1e-3;            % [V] Reduced Relative Inversion Offset.
        
        
        % ---------- Division Properties ----------
        
        % Define the absolute division subnetwork gains.
        c1_absolute_division_DEFAULT = 1.0;                         % [-] Absolute Division Gain 1.
        c2_absolute_division_DEFAULT = 1.0;                         % [-] Absolute Division Gain 2.
        c3_absolute_division_DEFAULT = 1.0;                         % [-] Absolute Division Gain 3.
        
        % Define the relative division subnetwork gains.
        c1_relative_division_DEFAULT = 1.0;                         % [-] Relative Division Gain 1.
        c2_relative_division_DEFAULT = 1.0;                         % [-] Relative Division Gain 2.
        c3_relative_division_DEFAULT = 1.0;                         % [-] Relative Division Gain 3.
        
        % Define division subnetwork offsets.
        delta_absolute_division_DEFAULT = 1e-3;                     % [V] Absolute Division Offset.
        delta_relative_division_DEFAULT = 1e-3;                     % [V] Relative Division Offset.
        
        
        % ---------- Reduced Division Properties ----------

        % Define the reduced absolute division subnetwork gains.
        c1_reduced_absolute_division_DEFAULT = 1.0;                 % [-] Reduced Absolute Division Gain 1.
        c2_reduced_absolute_division_DEFAULT = 1.0;                 % [-] Reduced Absolute Division Gain 2.
        
        % Define the reduced relative division subnetwork gains.
        c1_reduced_relative_division_DEFAULT = 1.0;                 % [-] Reduced Relative Division Gain 1.
        c2_reduced_relative_division_DEFAULT = 1.0;                 % [-] Reduced Relative Division Gain 2.
        
        % Define reduced division subnetwork offsets.
        delta_reduced_absolute_division_DEFAULT = 1e-3;             % [V] Reduced Absolute Division Offset.
        delta_reduced_relative_division_DEFAULT = 1e-3;          	% [V] Reduced Relative Division Offset.
        
        
        % ---------- Division After Inversion Properties ----------
        
        % Define the absolute division after inversion subnetwork gains.
        c1_absolute_dai_DEFAULT = 1.0;                              % [-] Absolute Division Gain 1.
        c2_absolute_dai_DEFAULT = 1.0;                              % [-] Absolute Division Gain 2.
        c3_absolute_dai_DEFAULT = 1.0;                              % [-] Absolute Division Gain 3.
        
        % Define the relative division after inversion subnetwork gains.
        c1_relative_dai_DEFAULT = 1.0;                              % [-] Relative Division Gain 1.
        c2_relative_dai_DEFAULT = 1.0;                              % [-] Relative Division Gain 2.
        c3_relative_dai_DEFAULT = 1.0;                              % [-] Relative Division Gain 3.
        
        % Define division after inversion subnetwork offsets.
        delta_absolute_dai_DEFAULT = 2e-3;                          % [V] Absolute Division After Inversion Offset.
        delta_relative_dai_DEFAULT = 2e-3;                          % [V] Relative Division After Inversion Offset.
        
        
        % ---------- Reduced Division After Inversion Properties ----------
        
        % Define the reduced absolute division after inversion subnetwork gains.
        c1_reduced_absolute_dai_DEFAULT = 1.0;                      % [-] Reduced Absolute Division Gain 1.
        c2_reduced_absolute_dai_DEFAULT = 1.0;                      % [-] Reduced Absolute Division Gain 2.
        
        % Define the reduced relative division after inversion subnetwork gains.
        c1_reduced_relative_dai_DEFAULT = 1.0;                      % [-] Reduced Relative Division Gain 1.
        c2_reduced_relative_dai_DEFAULT = 1.0;                      % [-] Reduced Relative Division Gain 2.
        
        % Define reduced division after inversion subnetwork offsets.
        delta_reduced_absolute_dai_DEFAULT = 2e-3;                 	% [V] Reduced Absolute Division After Inversion Offset.
        delta_reduced_relative_dai_DEFAULT = 2e-3;                 	% [V] Reduced Relative Division After Inversion Offset.
        

        % ---------- Derivation Properties ----------

        % Define derivative subnetwork parameters.
        c_derivation_DEFAULT = 1e6;                                 % [-] Default Derivative Gain.
        w_derivation_DEFAULT = 1;                                 	% [Hz?] Default Derivative Cutoff Frequency?
        sf_derivation_DEFAULT = 0.05;                            	% [-] Default Derivative safety Factor.
        
        
        % ---------- Integration Properties ----------
        
        % Define integration subnetwork parameters.
        c_integration_mean_DEFAULT = 0.01e9;                       	% [-] Default Average Integration Gain.
        
        
        % ---------- Central Pattern Generator Properties ----------

        % Define centeral pattern generator subnetwork parameters.
        T_oscillation_DEFAULT = 2;                              	% [s] Default Oscillation Period.
        r_oscillation_DEFAULT = 0.90;                           	% [-] Default Oscillation Decay.
        num_cpg_neurons_DEFAULT = 2;                              	% [#} Default Number of CPG Neurons.
       
        
        % ---------- Design Properties ----------

        % Define the default encoding scheme.
        encoding_scheme_DEFAULT = 'Absolute';                   	% [str] Default Encoding Scheme.
        
        % Define the default flags.
        enabled_flag_DEFAULT = true;                            	% [T/F] Default Enabled Flag.
        set_flag_DEFAULT = true;                                  	% [T/F] Default Set Flag (Determines whether to update the neuron object.)
        
    end
    
    
    %% NEURON METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = neuron_class( ID, name, U, h, Cm, Gm, Er, R, Am, Sm, dEm, Ah, Sh, dEh, dEna, tauh_max, Gna, Ileak, Isyn, Ina, Itonic, Iapp, Itotal, enabled_flag, neuron_utilities )
            
            % Set the default neuron properties.
            if nargin < 25, neuron_utilities = neuron_utilities_class(  ); end      % [class] Neuron Utilities Class.
            if nargin < 24, enabled_flag = self.enabled_flag_DEFAULT; end           % [T/F] Enable Flag.
            if nargin < 23, Itotal = self.Itotal_DEFAULT; end                       % [A] Total Current.
            if nargin < 22, Iapp = self.Iapp_DEFAULT; end                          	% [A] Applied Current.
            if nargin < 21, Itonic = self.Itonic_DEFAULT; end                       % [A] Tonic Current.
            if nargin < 20, Ina = self.Ina_DEFAULT; end                             % [A] Sodium Channel Current.
            if nargin < 19, Isyn = self.Isyn_DEFAULT; end                         	% [A] Synaptic Current.
            if nargin < 18, Ileak = self.Ileak_DEFAULT; end                         % [A] Leak Current.
            if nargin < 17, Gna = self.Gna_DEFAULT; end                             % [S] Sodium Channel Conductance.
            if nargin < 16, tauh_max = self.tauh_max_DEFAULT; end                   % [s] Maximum Sodium Channel Deactivation Time Constant.
            if nargin < 15, dEna = self.dEna_DEFAULT; end                           % [V] Sodium Channel Reveral Potential.
            if nargin < 14, dEh = self.dEh_DEFAULT; end                             % [V] Sodium Channel Deactivation Reversal Potential.
            if nargin < 13, Sh = self.Sh_DEFAULT; end                               % [V] Sodium Channel Deacitvation Slope.
            if nargin < 12, Ah = self.Ah_DEFAULT; end                               % [V] Sodium Channel Deactivation Amplitude.
            if nargin < 11, dEm = self.dEm_DEFAULT; end                             % [V] Sodium Channel Activation Reversal Potential.
            if nargin < 10, Sm = self.Sm_DEFAULT; end                               % [V] Sodium Channel Activation Slope.
            if nargin < 9, Am = self.Am_DEFAULT; end                                % [V] Sodium Channel Activation Amplitude.
            if nargin < 8, R = self.R_DEFAULT; end                                  % [V] Activation Domain.
            if nargin < 7, Er = self.Er_DEFAULT; end                                % [V] Membrane Reversal Potential.
            if nargin < 6, Gm = self.Gm_DEFAULT; end                                % [S] Membrane Conductance.
            if nargin < 5, Cm = self.Cm_DEFAULT; end                                % [C] Membrane Capacitance.
            if nargin < 4, h = self.h_DEFAULT; end                                	% [-] Sodium Channel Deactivation.
            if nargin < 3, U = self.U_DEFUALT; end                                	% [V] Membrane Voltage.
            if nargin < 2, name = self.name_DEFAULT; end                          	% [-] Neuron Name.
            if nargin < 1, ID = self.ID_DEFAULT; end                               	% [#] ID Number.
            
            % Store an instance of the neuron utilities class.
            self.neuron_utilities = neuron_utilities;
            
            % Store whether this neuron is active.
            self.enabled_flag = enabled_flag;
            
            % Store the current properties.
            self.Itotal = Itotal;
            self.Iapp = Iapp;
            self.Itonic = Itonic;
            self.Ina = Ina;
            self.Isyn = Isyn;
            self.Ileak = Ileak;
            
            % Store the ion channel properties.
            self.Gna = Gna;
            self.tauh_max = tauh_max;
            self.dEna = dEna;
            self.dEh = dEh;
            self.Sh = Sh;
            self.Ah = Ah;
            self.dEm = dEm;
            self.Sm = Sm;
            self.Am = Am;
            
            % Store the membrane properties.
            self.R = R;
            self.Er = Er;
            self.Gm = Gm;
            self.Cm = Cm;
            
            % Store the neuron states.
            self.h = h;
            self.U = U;
            
            % Store the neuron identification information.
            self.name = name;
            self.ID = ID;
            
            % Compute the steady state sodium channel activation and deactivation parameters.
            [ ~, self ] = self.compute_minf( U, Am, Sm, dEm, true, neuron_utilities );                      % [-] Steady State Sodium Channel Activation Parameter
            [ hinf, self ] = self.compute_hinf( U, Ah, Sh, dEh, true, neuron_utilities );                   % [-] Steady State Sodium Channel Deactivation Parameter
            
            % Determine whether to set the sodium channel activation parameter to its steady state value.
            if any( isnan( self.h ) ), self.h = hinf; end                                                  	% [-] Steady State Sodium Channel Deactivation Parameter
            
            % Compute and set the sodium channel deactivation time constant.
            [ ~, self ] = self.compute_tauh( U, tauh_max, hinf, Ah, Sh, dEh, true, neuron_utilities );      % [-] Sodium Channel Deactivation Time Constant                                                                         
            
        end
        
        
        %% Name Functions.
        
        % Implement a function to generate a name for this neuron.
        function [ name, self ] = generate_name( self, ID, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end
            if nargin < 2, ID = self.ID; end
            
            % Generate a name for the neuron.
            name = neuron_utilities.ID2name( ID );
            
            % Determine whether to update the name.
            if set_flag, self.name = name; end
            
        end
                
        
        %% Sodium Channel Activation & Deactivation Compute Functions.
        
        % Implement a function to compute the steady state sodium channel activation parameter.
        function [ minf, self ] = compute_minf( self, U, Am, Sm, dEm, set_flag, neuron_utilities )
            
            % Define the default input arguments.
            if nargin < 7, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 5, dEm = self.dEm; end                                                          % [V] Sodium Channel Activation Reversal Potential.
            if nargin < 4, Sm = self.Sm; end                                                         	% [-] Sodium Channel Activation Slope.
            if nargin < 3, Am = self.Am; end                                                         	% [-] Sodium Channel Activation Amplitude.
            if nargin < 2, U = self.U; end                                                           	% [V] Membrane Voltage.
            
            % Compute the steady state sodium channel activation parameter.
            minf = neuron_utilities.compute_mhinf( U, Am, Sm, dEm );                                	% [-] Sodium Channel Activation Parameter.
            
            % Determine whether to update the neuron object.
            if set_flag, self.minf = minf; end
            
        end
        
        
        % Implement a function to compute the steady state sodium channel deactivation parameter.
        function [ hinf, self ] = compute_hinf( self, U, Ah, Sh, dEh, set_flag, neuron_utilities )
            
            % Define the default input arguments.
            if nargin < 7, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 6, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)                 
            if nargin < 5, dEh = self.dEh; end                                                                          % [V] Sodium Channel Deactivation Reversal Potential
            if nargin < 4, Sh = self.Sh; end                                                                            % [-] Sodium Channel Deactivation Slope
            if nargin < 3, Ah = self.Ah; end                                                                            % [-] Sodium Channel Deactivation Amplitude
            if nargin < 2, U = self.U; end                                                                              % [V] Membrane Voltage
            
            % Compute the steady state sodium channel deactivaiton parameter.
            hinf = neuron_utilities.compute_mhinf( U, Ah, Sh, dEh );                                              % [-] Sodium Channel Deactivation Parameter
            
            % Determine whether to update the neuron object.
            if set_flag, self.hinf = hinf; end
            
        end
        
        
        % Implement a function to compute the sodium channel deactivation time constant.
        function [ tauh, self ] = compute_tauh( self, U, tauh_max, hinf, Ah, Sh, dEh, set_flag, neuron_utilities )
            
            % Define the default input arguments.
            if nargin < 9, neuron_utilities = self.neuron_utilities; end                     	% [class] Neuron Utilities.
            if nargin < 8, set_flag = self.set_flag_DEFAULT; end                              	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 7, dEh = self.dEh; end                                              	% [V] Sodium Channel Deactivation Reversal Potential.
            if nargin < 6, Sh = self.Sh; end                                                	% [-] Sodium Channel Deactivation Slope.
            if nargin < 5, Ah = self.Ah; end                                                    % [-] Sodium Channel Deactivation Amplitude.
            if nargin < 4, hinf = self.hinf; end                                              % [-] Steady State Sodium Channel Deactivation Parameter.
            if nargin < 3, tauh_max = self.tauh_max; end                                        % [s] Maximum Sodium Channel Deactivation Time Constant.
            if nargin < 2, U = self.U; end                                                      % [V] Membrane Voltage.
            
            % Compute the sodium channel deactivation time constant.
            tauh = neuron_utilities.compute_tauh( U, tauh_max, hinf, Ah, Sh, dEh );            % [s] Sodium Channel Deactivation Time Constant.
            
            % Determine whether to update the neuron object.
            if set_flag, self.tauh = tauh; end
            
        end
        
        
        %% Sodium Channel Conductance Compute Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a transmission subnetwork neuron.
        function [ Gna, self ] = compute_transmission_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end             	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                      	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end          % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this transmission subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue transmission subnetwork.
                Gna = neuron_utilities.compute_absolute_transmission_Gna(  );           % [S] Sodium Channel Conductance.
            
            elseif strcmpi( encoding_scheme, 'relative' )                               % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative transmission subnetwork.
                Gna = neuron_utilities.compute_relative_transmission_Gna(  );           % [S] Sodium Channel Conductance.

            else                                                                        % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for an addition subnetwork neuron.
        function [ Gna, self ] = compute_addition_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                 	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this addition subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue addition subnetwork.
                Gna = neuron_utilities.compute_absolute_addition_Gna(  );                   % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative addition subnetwork.
                Gna = neuron_utilities.compute_relative_addition_Gna(  );                   % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a subtraction subnetwork neuron.
        function [ Gna, self ] = compute_subtraction_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                 	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                           	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this subtraction subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue subtraction subnetwork.
                Gna = neuron_utilities.compute_absolute_subtraction_Gna(  );                % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative subtraction subnetwork.
                Gna = neuron_utilities.compute_relative_subtraction_Gna(  );                % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        

        % Implement a function to compute the required sodium channel conductance for a double subtraction subnetwork neuron.
        function [ Gna, self ] = compute_double_subtraction_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                 	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this double subtraction subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue double subtraction subnetwork.
                Gna = neuron_utilities.compute_absolute_double_subtraction_Gna(  );         % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative double subtraction subnetwork.
                Gna = neuron_utilities.compute_relative_double_subtraction_Gna(  );         % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
      
        % Implement a function to compute the required sodium channel conductance for an inversion subnetwork neuron.
        function [ Gna, self ] = compute_inversion_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                           	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue inversion subnetwork.
                Gna = neuron_utilities.compute_absolute_inversion_Gna(  );                  % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative inversion subnetwork.
                Gna = neuron_utilities.compute_relative_inversion_Gna(  );                  % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a reduced inversion subnetwork neuron.
        function [ Gna, self ] = compute_reduced_inversion_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                           	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue inversion subnetwork.
                Gna = neuron_utilities.compute_reduced_absolute_inversion_Gna(  );                  % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative inversion subnetwork.
                Gna = neuron_utilities.compute_reduced_relative_inversion_Gna(  );                  % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Division Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a division subnetwork neuron.
        function [ Gna, self ] = compute_division_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue division subnetwork.
                Gna = neuron_utilities.compute_absolute_division_Gna(  );                   % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative division subnetwork.
                Gna = neuron_utilities.compute_relative_division_Gna(  );                   % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a reduced division subnetwork neuron.
        function [ Gna, self ] = compute_reduced_division_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue division subnetwork.
                Gna = neuron_utilities.compute_reduced_absolute_division_Gna(  );                   % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative division subnetwork.
                Gna = neuron_utilities.compute_reduced_relative_division_Gna(  );                   % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a division after inversion subnetwork neuron.
        function [ Gna, self ] = compute_dai_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue division subnetwork.
                Gna = neuron_utilities.compute_absolute_dai_Gna(  );                   % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative division subnetwork.
                Gna = neuron_utilities.compute_relative_dai_Gna(  );                   % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
                
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a reduced division after inversion subnetwork neuron.
        function [ Gna, self ] = compute_reduced_dai_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue division subnetwork.
                Gna = neuron_utilities.compute_reduced_absolute_dai_Gna(  );                   % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative division subnetwork.
                Gna = neuron_utilities.compute_reduced_relative_dai_Gna(  );                   % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a multiplication subnetwork neuron.
        function [ Gna, self ] = compute_multiplication_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                   	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                         	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this multiplication subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue multiplication subnetwork.
                Gna = neuron_utilities.compute_absolute_multiplication_Gna(  );             % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative multiplication subnetwork.
                Gna = neuron_utilities.compute_relative_multiplication_Gna(  );             % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for a reduced multiplication subnetwork neuron.
        function [ Gna, self ] = compute_reduced_multiplication_Gna( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                   	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                         	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the sodium channel conductance for this multiplication subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the sodium channel conductance for this neuron assuming that it belongs to an absolue multiplication subnetwork.
                Gna = neuron_utilities.compute_reduced_absolute_multiplication_Gna(  );             % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the sodium channel conductance for this neuron assuming that it belongs to a relative multiplication subnetwork.
                Gna = neuron_utilities.compute_reduced_relative_multiplication_Gna(  );             % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Derivation Subnetwork Functions ----------
                
        % Implement a function to compute the required sodium channel conductance for a derivation subnetwork neuron.
        function [ Gna, self ] = compute_derivation_Gna( self, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 3, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Compute the sodium channel conductance for a derivation subnetwork neuron.
            Gna = neuron_utilities.compute_derivation_Gna(  );                                      	% [S] Sodium Channel Conductance
           
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        % Implement a function to compute the required sodium channel conductance for an integration subnetwork neuron.
        function [ Gna, self ] = compute_integration_Gna( self, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 3, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Compute the sodium channel conductance for an integration subnetwork neuron.
            Gna = neuron_utilities.compute_integration_Gna(  );                                      	% [S] Sodium Channel Conductance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % Implement a function to compute the required sodium channel conductance for a voltage based integration subnetwork neuron.
        function [ Gna, self ] = compute_vbi_Gna( self, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 3, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Compute the sodium channel conductance for a voltage based integration subnetwork neuron.
            Gna = neuron_utilities.compute_vbi_Gna(  );                                                 % [S] Sodium Channel Conductance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % Implement a function to compute the required sodium channel conductance for a split voltage based integration subnetwork neuron.
        function [ Gna, self ] = compute_svbi_Gna( self, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 3, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Compute the sodium channel conductance for a split voltage based integration subnetwork neuron.
            Gna = neuron_utilities.compute_svbi_Gna(  );                                                % [S] Sodium Channel Conductance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        % ---------- Central Pattern Generator Subnetwork Functions ----------

        % Implement a function to compute the required sodium channel conductance to create oscillation in a CPG subnetwork.
        function [ Gna, self ] = compute_cpg_Gna( self, R, Gm, Am, Sm, dEm, Ah, Sh, dEh, dEna, set_flag, neuron_utilities )
            
            % Define the default input arguments.
            if nargin < 12, neuron_utilities = self.neuron_utilities; end                                               % [class] Neuron Utilities.
            if nargin < 11, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 10, dEna = self.dEna; end                                                                       % [V] Sodium Channel Reversal Potential
            if nargin < 9, dEh = self.dEh; end                                                                          % [V] Sodium Channel Deactivation Reversal Potential
            if nargin < 8, Sh = self.Sh; end                                                                            % [-] Sodium Channel Deactivation Slope
            if nargin < 7, Ah = self.Ah; end                                                                            % [-] Sodium Channel Deactivation Amplitude
            if nargin < 6, dEm = self.dEm; end                                                                          % [V] Sodium Channel Activation Reversal Potential
            if nargin < 5, Sm = self.Sm; end                                                                            % [-] Sodium Channel Activation Slope
            if nargin < 4, Am = self.Am; end                                                                            % [-] Sodium Channel Activation Amplitude
            if nargin < 3, Gm = self.Gm; end                                                                            % [S] Membrane Conductance
            if nargin < 2, R = self.R; end                                                                              % [V] Activation Domain
            
            % Compute the required sodium channel conductance to create oscillation in a two neuron CPG subnetwork.
            Gna = neuron_utilities.compute_cpg_Gna( R, Gm, Am, Sm, dEm, Ah, Sh, dEh, dEna );                            % [S] Sodium Channel Conductance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end

        end
        
        
        % Implement a function to compute the required sodium channel conductance for a driven multistate cpg subnetwork neuron.
        function [ Gna, self ] = compute_dmcpg_Gna( self, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 3, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)

            % Compute the sodium channel conductance for a driven multistate cpg subnetwork neuron.
            Gna = neuron_utilities.compute_dmcpg_Gna(  );                                               % [S] Sodium Channel Conductance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gna = Gna; end
            
        end
        
        
        %% Membrane Conductance Compute Functions.
        
        %{
        
        % Implement a function to compute the membrance conductance for addition subnetwork input neurons.
        function [ Gm, self ] = compute_addition_Gm_input( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                        	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this addition subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an input to an absolue addition subnetwork.
                Gm = neuron_utilities.compute_absolute_addition_Gm_input(  );               % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an input to a relative addition subnetwork.
                Gm = neuron_utilities.compute_relative_addition_Gm_input(  );               % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        
        
        % Implement a function to compute the membrance conductance for addition subnetwork output neurons.
        function [ Gm, self ] = compute_addition_Gm_output( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                   	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this addition subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an output to an absolue addition subnetwork.
                Gm = neuron_utilities.compute_absolute_addition_Gm_output(  );              % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an output to a relative addition subnetwork.
                Gm = neuron_utilities.compute_relative_addition_Gm_output(  );              % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        
        
        % Implement a function to compute the membrance conductance for subtraction subnetwork input neurons.
        function [ Gm, self ] = compute_subtraction_Gm_input( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                   	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                           	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this subtraction subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an input to an absolue subtraction subnetwork.
                Gm = neuron_utilities.compute_absolute_subtraction_Gm_input(  );            % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an input to a relative subtraction subnetwork.
                Gm = neuron_utilities.compute_relative_subtraction_Gm_input(  );            % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        
        
        % Implement a function to compute the membrance conductance for subtraction subnetwork output neurons.
        function [ Gm, self ] = compute_subtraction_Gm_output( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                        	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this subtraction subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an input to an absolue subtraction subnetwork.
                Gm = neuron_utilities.compute_absolute_subtraction_Gm_output(  );           % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an input to a relative subtraction subnetwork.
                Gm = neuron_utilities.compute_relative_subtraction_Gm_output(  );           % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        

        % Implement a function to compute the membrance conductance for inversion subnetwork input neurons.
        function [ Gm, self ] = compute_inversion_Gm_input( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an input to an absolue inversion subnetwork.
                Gm = neuron_utilities.compute_absolute_inversion_Gm_input(  );              % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an input to a relative inversion subnetwork.
                Gm = neuron_utilities.compute_relative_inversion_Gm_input(  );              % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        
        
        % Implement a function to compute the membrance conductance for inversion subnetwork output neurons.
        function [ Gm, self ] = compute_inversion_Gm_output( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an output to an absolue inversion subnetwork.
                Gm = neuron_utilities.compute_absolute_inversion_Gm_output(  );             % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an output to a relative inversion subnetwork.
                Gm = neuron_utilities.compute_relative_inversion_Gm_output(  );             % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        
        
        % Implement a function to compute the membrance conductance for division subnetwork input neurons.
        function [ Gm, self ] = compute_division_Gm_input( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                           	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an input to an absolue division subnetwork.
                Gm = neuron_utilities.compute_absolute_division_Gm_input(  );               % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an input to a relative division subnetwork.
                Gm = neuron_utilities.compute_relative_division_Gm_input(  );               % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        
        
        % Implement a function to compute the membrance conductance for division subnetwork output neurons.
        function [ Gm, self ] = compute_division_Gm_output( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                            % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane conductance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane conductance for this neuron assuming that it is an output to an absolue division subnetwork.
                Gm = neuron_utilities.compute_absolute_division_Gm_output(  );              % [S] Sodium Channel Conductance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                 	% If the encoding scheme is set to relative...
            
                % Compute the membrane conductance for this neuron assuming that it is an output to a relative division subnetwork.
                Gm = neuron_utilities.compute_relative_division_Gm_output(  );              % [S] Sodium Channel Conductance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        

        % Implement a function to compute the required membrane conductance for a derivation neuron.
        function [ Gm, self ] = compute_derivation_Gm( self, k, w, safety_factor, set_flag, neuron_utilities )
            
            % Set the default input arugments.
            if nargin < 6, neuron_utilities = self.neuron_utilities; end                 	% [class] Neuron Utilities.
            if nargin < 5, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 4, safety_factor = self.sf_derivation_DEFAULT; end               	% [-] Derivation Subnetwork Safety Factor.
            if nargin < 3, w = self.w_derivation_DEFAULT; end                             	% [Hz?] Derivation Subnetwork Cutoff Frequency?
            if nargin < 2, k = self.c_derivation_DEFAULT; end                             	% [-] Derivation Subnetwork Gain.
            
            % Compute the membrane conductance for this derivation neuron.
            Gm = neuron_utilities.compute_derivation_Gm( k, w, safety_factor );             % [S] Membrane Conductance.
            
            % Determine whether to update the neuron object.
            if set_flag, self.Gm = Gm; end
            
        end
        
        %}
        
        
        %% Membrane Capacitance Compute Functions.
        
        %{
        
        % Implement a function to compute the required membrane capacitance for a transmission subnetwork neuron.
        function [ Cm, self ] = compute_transmission_Cm( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end         	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                 	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end     	% [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane capacitance for this transmission subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                               % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue transmission subnetwork.
                Cm = neuron_utilities.compute_absolute_transmission_Cm(  );         % [C] Membrane Capacitance
            
            elseif strcmpi( encoding_scheme, 'relative' )                           % If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative transmission subnetwork.
                Cm = neuron_utilities.compute_relative_transmission_Cm(  );         % [C] Membrane Capacitance

            else                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for a transmission subnetwork neuron.
        function [ Cm, self ] = compute_slow_transmission_Cm( self, Gm, num_cpg_neurons, T, r, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 8, neuron_utilities = self.neuron_utilities; end                                                        % [class] Neuron Utilities.
            if nargin < 7, set_flag = self.set_flag_DEFAULT; end                                                                % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 6, encoding_scheme = self.encoding_scheme_DEFAULT; end                                                  % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 5, r = self.r_oscillation_DEFAULT; end                                                                  % [-] Oscillation Decay
            if nargin < 4, T = self.T_oscillation_DEFAULT; end                                                                  % [s] Oscillation Period
            if nargin < 3, num_cpg_neurons = self.num_cpg_neurons_DEFAULT; end                                                  % [#] Number of CPG Neurons
            if nargin < 2, Gm = self.Gm; end                                                                                    % [S] Membrane Conductance
            
            % Determine how to compute the membrane capacitance for this transmission subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                                           % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue transmission subnetwork.
                Cm = neuron_utilities.compute_slow_absolute_transmission_Cm( Gm, num_cpg_neurons, T, r );                       % [C] Membrane Capacitance

            elseif strcmpi( encoding_scheme, 'relative' )                                                                   	% If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative transmission subnetwork.
                Cm = neuron_utilities.compute_slow_relative_transmission_Cm( Gm, num_cpg_neurons, T, r );                       % [C] Membrane Capacitance

            else                                                                                                                % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for a modulation subnetwork neuron.
        function [ Cm, self ] = compute_modulation_Cm( self, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 3, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                                    % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Compute the membrane capacitance for a modulation subnetwork neuron.
            Cm = neuron_utilities.compute_modulation_Cm(  );                                                        % [C] Membrane Capacitance
           
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for an addition subnetwork neuron.
        function [ Cm, self ] = compute_addition_Cm( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane capacitance for this addition subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue addition subnetwork.
                Cm = neuron_utilities.compute_absolute_addition_Cm(  );                     % [C] Membrane Capacitance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative addition subnetwork.
                Cm = neuron_utilities.compute_relative_addition_Cm(  );                     % [C] Membrane Capacitance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for a subtraction subnetwork neuron.
        function [ Cm, self ] = compute_subtraction_Cm( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                   	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                         	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane capacitance for this subtraction subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue subtraction subnetwork.
                Cm = neuron_utilities.compute_absolute_subtraction_Cm(  );                  % [C] Membrane Capacitance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                 	% If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative subtraction subnetwork.
                Cm = neuron_utilities.compute_relative_subtraction_Cm(  );                  % [C] Membrane Capacitance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        

        % Implement a function to compute the membrane capacitance for a double subtraction subnetwork neuron.
        function [ Cm, self ] = compute_double_subtraction_Cm( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                   	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                           	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane capacitance for this double subtraction subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue double subtraction subnetwork.
                Cm = neuron_utilities.compute_absolute_double_subtraction_Cm(  );           % [C] Membrane Capacitance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative double subtraction subnetwork.
                Cm = neuron_utilities.compute_relative_double_subtraction_Cm(  );           % [C] Membrane Capacitance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for an inversion subnetwork neuron.
        function [ Cm, self ] = compute_inversion_Cm( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                          	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane capacitance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue inversion subnetwork.
                Cm = neuron_utilities.compute_absolute_inversion_Cm(  );                    % [C] Membrane Capacitance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                  	% If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative inversion subnetwork.
                Cm = neuron_utilities.compute_relative_inversion_Cm(  );                    % [C] Membrane Capacitance

            else                                                                          	% Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        

        % Implement a function to compute the membrane capacitance for a division subnetwork neuron.
        function [ Cm, self ] = compute_division_Cm( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                         	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue division subnetwork.
                Cm = neuron_utilities.compute_absolute_division_Cm(  );                     % [C] Membrane Capacitance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                   % If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative division subnetwork.
                Cm = neuron_utilities.compute_relative_division_Cm(  );                     % [C] Membrane Capacitance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for a multiplication subnetwork neuron.
        function [ Cm, self ] = compute_multiplication_Cm( self, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                  	% [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                           	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end              % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            
            % Determine how to compute the membrane capacitance for this multiplication subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                       % If the encoding scheme is set to absolute...

                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue multiplication subnetwork.
                Cm = neuron_utilities.compute_absolute_multiplication_Cm(  );               % [C] Membrane Capacitance
            
            elseif strcmpi( encoding_scheme, 'relative' )                                 	% If the encoding scheme is set to relative...
            
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative multiplication subnetwork.
                Cm = neuron_utilities.compute_relative_multiplication_Cm(  );               % [C] Membrane Capacitance

            else                                                                            % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the first membrane capacitance for a derivation subnetwork neuron.
        function [ Cm1, self ] = compute_derivation_Cm1( self, Gm, Cm2, k, set_flag, neuron_utilities  )
            
            % Set the default input arguments.
            if nargin < 6, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 5, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 4, k = self.c_derivation_DEFAULT; end                                          	% [-] Derivative Subnetwork Gain.
            if nargin < 3, Cm2 = 1e-9; end                                                            	% [C] Membrane Capacitance.
            if nargin < 2, Gm = self.Gm; end                                                         	% [S] Membrance Conductance.
            
            % Compute the first membrane capacitance for a derivation subnetwork neuron.
            Cm1 = neuron_utilities.compute_derivation_Cm1( Gm, Cm2, k );                              	% [C] Membrane Capacitance.
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm1; end
            
        end
        
        
        % Implement a function to compute the second membrane capacitance for a derivation subnetwork neuron.
        function [ Cm2, self ] = compute_derivation_Cm2( self, Gm, w, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, w = self.w_derivation_DEFAULT; end                                          	% [Hz?] Derivative Subnetwork Cutoff Frequency?
            if nargin < 2, Gm = self.Gm; end                                                          	% [S] Membrane Conductance.
            
            % Compute the second membrane capacitance for a derivation subnetwork neuron.
            Cm2 = neuron_utilities.compute_derivation_Cm2( Gm, w );                                  	% [C] Membrane Capacitance.
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm2; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for an integration neuron.
        function [ Cm, self ] = compute_integration_Cm( self, ki_mean, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, ki_mean = self.c_integration_mean_DEFAULT; end                             	% [-] Average Integration Gain
            
            % Compute the membrane capacitance for this integration neuron.
            Cm = neuron_utilities.compute_integration_Cm( ki_mean );                                  	% [C] Membrane Capacitance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to compute the membrane capacitance for a voltage based integration neuron.
        function [ Cm, self ] = compute_vbi_Cm( self, ki_mean, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, ki_mean = self.c_integration_mean_DEFAULT; end                             	% [-] Average Integration Gain
            
            % Compute the membrane capacitance for this voltage based integration neuron.
            Cm = neuron_utilities.compute_vbi_Cm( ki_mean );                                            % [C] Membrane Capacitance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to first compute the membrane capacitance for a split voltage based integration neuron.
        function [ Cm, self ] = compute_svbi_Cm1( self, ki_mean, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 4, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 2, ki_mean = self.c_integration_mean_DEFAULT; end                              	% [-] Average Integration Gain
            
            % Compute the first membrane capacitance for this split voltage based integration neuron.
            Cm = neuron_utilities.compute_svbi_Cm1( ki_mean );                                          % [C] Membrane Capacitance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        
        % Implement a function to second compute the membrane capacitance for a split voltage based integration neuron.
        function [ Cm, self ] = compute_svbi_Cm2( self, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 3, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Compute the second membrane capacitance for this split voltage based integration neuron.
            Cm = neuron_utilities.compute_svbi_Cm2(  );                                                 % [C] Membrane Capacitance
            
            % Determine whether to update the neuron object.
            if set_flag, self.Cm = Cm; end
            
        end
        
        %}
        
        
        %% Parameter Unpacking Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute transmission output activation domain.
        function [ c, R1 ] = unpack_absolute_transmission_Rn_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                c = self.c_absolute_transmission_DEFAULT;               % [-] Absolute Transmission Gain.
                R1 = self.R_DEFAULT;                                    % [V] Activation Domain
                
            elseif length( parameters ) == 2                          	% If there are a specific number of parameters...
                
                % Unpack the parameters.
                c = parameters{ 1 };                                    % [-] Absolute Transmission Gain.
                R1 = parameters{ 2 };                                	% [V] Activation Domain
                
            else                                                     	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute addition output activation domain.
        function [ cs, Rs ] = unpack_absolute_addition_Rn_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                               % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                            % If the parameters are empty...
                
                % Set the parameters to default values.
                cs = self.c_absolute_addition_DEFAULT*ones( 1, 2 );             % [-] Absolute Addition Gain.
                Rs = self.R_DEFAULT*ones( 1, 2 );                               % [V] Activation Domain
                
            elseif length( parameters ) == 2                                    % If there are a specific number of parameters...
                
                % Unpack the parameters.
                cs = parameters{ 1 };
                Rs = parameters{ 2 };                                           % [V] Activation Domain.
                
            else                                                                % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end
            
        end
        
    
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute subtraction output activation domain.
        function [ cs, s_ks, Rs ] = unpack_absolute_subtraction_Rn_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                       % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                    % If the parameters are empty...
                
                % Set the parameters to default values.
                cs = self.c_absolute_subtraction_DEFAULT;               % [-] Absolute Subtraction Gain.
                s_ks = self.s_ks_DEFAULT;                               % [-] Subtraction Signature.
                Rs = self.R_DEFAULT*ones( 1, 2 );                       % [V] Activation Domain.

            elseif length( parameters ) == 2                            % If there are a specific number of parameters...
                
                % Unpack the parameters.
                cs = parameters{ 1 };                                   % [-] Absolute Subtraction Gain.
                s_ks = parameters{ 2 };                                 % [-] Subtraction Signature.
                Rs = parameters{ 3 };                                   % [V] Activation Domain.
                                    
            else                                                        % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
                
        % Implement a function to unpack the parameters required to compute the absolute inversion output activation domain.
        function [ c1, c3 ] = unpack_absolute_inversion_R2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                   % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_absolute_inversion_DEFAULT;          	% [-] General Subnetwork Gain 1.
                c3 = self.c3_absolute_inversion_DEFAULT;         	% [-] General Subnetwork Gain 3.

            elseif length( parameters ) == 2                     	% If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                            	% [-] General Subnetwork Gain 1.
                c3 = parameters{ 2 };                               % [-] General Subnetwork Gain 3.
                
            else                                                    % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the relative inversion output activation domain.
        function R2 = unpack_relative_inversion_R2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                   % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...

                % Set the default parameters.
                R2 = self.R2_relative_inversion_DEFAULT;          	% [V] Maximum Membrane Voltage.

            elseif length( parameters ) == 1                     	% If there are a specific number of parameters...

                % Retrieve the parameters.
                R2 = parameters{ 1 };                            	% [V] Maximum Membrane Voltage.
                
            else                                                    % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the reduced absolute inversion output activation domain.
        function [ c1, c2 ] = unpack_reduced_absolute_inversion_R2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_reduced_absolute_inversion_DEFAULT;          	% [-] General Subnetwork Gain 1.
                c2 = self.c2_reduced_absolute_inversion_DEFAULT;            % [-] General Subnetwork Gain 3.

            elseif length( parameters ) == 2                                % If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                                       % [-] General Subnetwork Gain 1.
                c2 = parameters{ 2 };                                       % [-] General Subnetwork Gain 3.
                
            else                                                            % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the reduced absolute inversion output activation domain.
        function R2 = unpack_reduced_relative_inversion_R2_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...

                % Set the default parameters.
                R2 = self.R2_relative_inversion_DEFAULT;          	% [V] Maximum Membrane Voltage.

            elseif length( parameters ) == 1                     	% If there are a specific number of parameters...

                % Retrieve the parameters.
                R2 = parameters{ 1 };                            	% [V] Maximum Membrane Voltage.
                
            else                                                    % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end               
            
        end
        
        
        % ---------- Division Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute division output activation domain.
        function [ c1, c3, R1 ] = unpack_absolute_division_R3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end               	% [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_absolute_division_DEFAULT;             % [-] Subnetwork Gain 1.
                c3 = self.c3_absolute_division_DEFAULT;           	% [-] Subnetwork Gain 3.
                R1 = self.R_DEFAULT;                                % [V] Activation Domain 1.

            elseif length( parameters ) == 3                      	% If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                              	% [-] Subnetwork Gain 1.
                c3 = parameters{ 2 };                               % [-] Subnetwork Gain 3.
                R1 = parameters{ 3 };                               % [V] Activation Domain 1.
                
            else                                                   	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end

        
        % Implement a function to unpack the parameters required to compute the absolute division after inversion output activation domain.
        function [ c1, c2, c3, delta1, R1 ] = unpack_absolute_dai_R3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                   % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_absolute_dai_DEFAULT;                                  % [-] Absolute Division After Inversion Subnetwork Gain 1.
                c2 = self.c2_absolute_dai_DEFAULT;                                  % [-] Absolute Division After Inversion Subnetwork Gain 2.
                c3 = self.c3_absolute_dai_DEFAULT;                                  % [-] Absolute Division After Inversion Subnetwork Gain 3.
                delta1 = self.delta_absolute_inversion_DEFAULT;                     % [V] Absolute Inversion Subnetwork Offest.
                R1 = self.R_DEFAULT;                                                % [V] Activation Domain 1.

            elseif length( parameters ) == 5                                        % If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                                               % [-] Absolute Division After Inversion Subnetwork Gain 1.
                c2 = parameters{ 2 };                                               % [-] Absolute Division After Inversion Subnetwork Gain 2.
                c3 = parameters{ 3 };                                               % [-] Absolute Division After Inversion Subnetwork Gain 3.
                delta1 = parameters{ 4 };                                           % [V] Absolute Inversion Subnetwork Offset.
                R1 = parameters{ 5 };                                               % [V] Activation Domain 1.
                
            else                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
                
        % ---------- Reduced Division Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the reduced absolute division output activation domain.
        function [ c1, c2, R1 ] = unpack_reduced_absolute_division_R3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end               	% [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_absolute_division_DEFAULT;           	% [-] Absolute Division Subnetwork Gain 1.
                c2 = self.c2_absolute_division_DEFAULT;         	% [-] Absolute Division Subnetwork Gain 2.
                R1 = self.R_DEFAULT;                                % [V] Activation Domain 1.

            elseif length( parameters ) == 3                      	% If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                              	% [-] Absolute Division Subnetwork Gain 1.
                c2 = parameters{ 2 };                               % [-] Absolute Division Subnetwork Gain 2.
                R1 = parameters{ 3 };                               % [V] Activation Domain 1.
                
            else                                                   	% Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the reduced absolute division after inversion output activation domain.
        function [ c1, c2, delta1, R1 ] = unpack_reduced_absolute_dai_R3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                                   % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                                % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_absolute_dai_DEFAULT;                                  % [-] Absolute Division After Inversion Subnetwork Gain 1.
                c2 = self.c2_absolute_dai_DEFAULT;                                  % [-] Absolute Division After Inversion Subnetwork Gain 2.
                delta1 = self.delta_absolute_inversion_DEFAULT;                     % [V] Absolute Inversion Subnetwork Offest.
                R1 = self.R_DEFAULT;                                                % [V] Activation Domain 1.

            elseif length( parameters ) == 4                                        % If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                                               % [-] Absolute Division After Inversion Subnetwork Gain 1.
                c2 = parameters{ 2 };                                               % [-] Absolute Division After Inversion Subnetwork Gain 2.
                delta1 = parameters{ 3 };                                           % [V] Absolute Inversion Subnetwork Offset.
                R1 = parameters{ 4 };                                               % [V] Activation Domain 1.
                
            else                                                                    % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to unpack the parameters required to compute the absolute multiplication neuron 3 activation domain.
        function [ c1, c3 ] = unpack_absolute_multiplication_R3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_absolute_inversion_DEFAULT;                    % [-] Absolute Inversion Gain 1.
                c3 = self.c3_absolute_inversion_DEFAULT;                    % [-] Absolute Inversion Gain 3.

            elseif length( parameters ) == 7                                % If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                                       % [-] Absolute Inversion Gain 1.
                c3 = parameters{ 2 };                                       % [-] Absolute Inversion Gain 3.
                
            else                                                            % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the absolute multiplication neuron 4 activation domain.
        function [ c4, c5, c6, delta1, R1 ] = unpack_absolute_multiplication_R4_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...

                % Set the default parameters.
                c4 = self.c1_absolute_division_DEFAULT;                     % [-] Absolute Division Gain 1.
                c5 = self.c2_absolute_division_DEFAULT;                     % [-] Absolute Division Gain 2.
                c6 = self.c3_absolute_division_DEFAULT;                     % [-] Absolute Division Gain 3.
                delta1 = self.delta_absolute_inversion_DEFAULT;             % [V] Absolute Inversion Offset.
                R1 = self.R_DEFAULT;                                        % [V] Activation Domain.

            elseif length( parameters ) == 5                                % If there are a specific number of parameters...

                % Retrieve the parameters.
                c4 = parameters{ 1 };                                       % [-] Absolute Division Gain 1.
                c5 = parameters{ 2 };                                       % [-] Absolute Division Gain 2.
                c6 = parameters{ 3 };                                       % [-] Absolute Division Gain 3.
                delta1 = parameters{ 4 };                                   % [V] Absolute Inversion Offset.
                R1 = parameters{ 5 };                                       % [V] Activation Domain.
                
            else                                                            % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to unpack the parameters required to compute the reduced absolute multiplication neuron 3 activation domain.
        function [ c1, c2 ] = unpack_reduced_absolute_multiplication_R3_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...

                % Set the default parameters.
                c1 = self.c1_reduced_absolute_inversion_DEFAULT;            % [-] Reduced Absolute Inversion Gain 1.
                c2 = self.c2_reduced_absolute_inversion_DEFAULT;          	% [-] Reduced Absolute Inversion Gain 3.

            elseif length( parameters ) == 7                                % If there are a specific number of parameters...

                % Retrieve the parameters.
                c1 = parameters{ 1 };                                       % [-] Reduced Absolute Inversion Gain 1.
                c2 = parameters{ 2 };                                       % [-] Reduced Absolute Inversion Gain 3.
                
            else                                                            % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        % Implement a function to unpack the parameters required to compute the reduced absolute multiplication neuron 4 activation domain.
        function [ c3, c4, delta1, R1 ] = unpack_reduced_absolute_multiplication_R4_parameters( self, parameters )
        
            % Set the default input arguments.
            if nargin < 2, parameters = {  }; end                           % [-] Parameters Cell.
            
            % Determine how to set the parameters.
            if isempty( parameters )                                        % If the parameters are empty...

                % Set the default parameters.
                c3 = self.c1_reduced_absolute_division_DEFAULT;             % [-] Reduced Absolute Division Gain 1.
                c4 = self.c2_reduced_absolute_division_DEFAULT;             % [-] Reduced Absolute Division Gain 2.
                delta1 = self.delta_absolute_inversion_DEFAULT;             % [V] Absolute Inversion Offset.
                R1 = self.R_DEFAULT;                                        % [V] Activation Domain.

            elseif length( parameters ) == 4                                % If there are a specific number of parameters...

                % Retrieve the parameters.
                c3 = parameters{ 1 };                                       % [-] Reduced Absolute Division Gain 1.
                c4 = parameters{ 2 };                                       % [-] Reduced Absolute Division Gain 2.
                delta1 = parameters{ 3 };                                   % [V] Absolute Inversion Offset.
                R1 = parameters{ 4 };                                       % [V] Activation Domain.
                
            else                                                            % Otherwise...
                
                % Throw an error.
                error( 'Unable to unpack parameters.' )
                
            end            
            
        end
        
        
        %% Activation Domain Compute Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the transmission output neuron.
        function [ R2, self ] = compute_transmission_R2( self, parameters, encoding_scheme, set_flag, neuron_utilities )
        
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end           	% [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                  	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end
            
            % Determine how to compute the membrane capacitance for this addition subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                               % If the encoding scheme is set to absolute...

                % Unpack the absolute transmission parameters.
                [ c, R1 ] = self.unpack_absolute_transmission_Rn_parameters( parameters );
                
                % Compute the activation domain for this neuron assuming that it belongs to an absolue addition subnetwork.
                R2 = neuron_utilities.compute_absolute_transmission_R2( c, R1 );         % [V] Activation Domain.
            
            elseif strcmpi( encoding_scheme, 'relative' )                           % If the encoding scheme is set to relative...
            
                % Retrieve the maximum membrane voltage.
                R2 = parameters{ 1 };

                % % Throw an error.
                % error( 'R2 is a free parameter for relative transmission subnetworks.' )

            else                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R2; end
                    
        end
        
        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the addition subnetwork output neuron.
        function [ Rn, self ] = compute_addition_Rn( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end              	% [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                       	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end          % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                       % [-] Parameters Cell.                                                                             

            % Determine how to compute the membrane capacitance for this addition subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                   % If the encoding scheme is set to absolute...

                % Unpack the parameters required to compute the absolute addition subnetwork output activation domain.                    
                [ cs, Rs ] = self.unpack_absolute_addition_Rn_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue addition subnetwork.
                Rn = neuron_utilities.compute_absolute_addition_Rn( cs, Rs );          % [V] Activation Domain.
            
            elseif strcmpi( encoding_scheme, 'relative' )                               % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'Rn is a free parameter for relative addition subnetworks.' )

            else                                                                        % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = Rn; end
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the subtraction subnetwork output neurons.
        function [ Rn, self ] = compute_subtraction_Rn( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                    % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this subtraction subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                               % If the encoding scheme is set to absolute...

                % Unpack the parameters required to compute the absolute subtraction subnetwork output activation domain.
                [ cs, s_ks, Rs ] = self.unpack_absolute_subtraction_Rn_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue subtraction subnetwork.            
                Rn = neuron_utilities.compute_absolute_subtraction_Rn( cs, s_ks, Rs );              % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'Rn is a free parameter for relative subtraction subnetworks.' )
                
            else                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = Rn; end
            
        end

        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the inversion subnetwork output neuron.
        function [ R2, self ] = compute_inversion_R2( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                         	% [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                       % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                   % If the encoding scheme is set to absolute...

                % Unpack the parameters required to compute the absolute inversion subnetwork output activation domain.
                [ c1, c3 ] = self.unpack_absolute_inversion_R2_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue inversion subnetwork.            
                R2 = neuron_utilities.compute_absolute_inversion_R2( c1, c3 );                          % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                               % If the encoding scheme is set to relative...
            
                % Unpack the parameters required to compute the relative inversion subnetwork output activation domain.
                R2 = self.unpack_relative_inversion_R2_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to a relative inversion subnetwork.            
                R2 = neuron_utilities.compute_relative_inversion_R2( R2 );                          % [V] Activation Domain.

            else                                                                                        % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R2; end

        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the reduced inversion subnetwork output neuron.
        function [ R2, self ] = compute_reduced_inversion_R2( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                         	% [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                       % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this inversion subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                   % If the encoding scheme is set to absolute...

                % Unpack the parameters required to compute the absolute inversion subnetwork output activation domain.
                [ c1, c2 ] = self.unpack_reduced_absolute_inversion_R2_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue inversion subnetwork.            
                R2 = neuron_utilities.compute_reduced_absolute_inversion_R2( c1, c2 );                          % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                               % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R2 is a free parameter for reduced relative inversion subnetworks.' )

            else                                                                                        % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R2; end

        end


        % ---------- Division Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the division subnetwork output neurons.
        function [ R3, self ] = compute_division_R3( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute division subnetwork output activation domain.
                [ c1, c3, R1 ] = self.unpack_absolute_division_R3_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue division subnetwork.            
                R3 = neuron_utilities.compute_absolute_division_R3( c1, c3, R1 );                                   % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R3 is a free parameter for relative division subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R3; end
                        
        end
        

        % ---------- Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the division after inversion subnetwork output neurons.
        function [ R3, self ] = compute_dai_R3( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute division subnetwork output activation domain.
                [ c1, c2, c3, delta1, R1 ] = self.unpack_absolute_dai_R3_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue division subnetwork.            
                R3 = neuron_utilities.compute_absolute_dai_R3( c1, c2, c3, delta1, R1 );                                   % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R3 is a free parameter for relative division subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R3; end
                        
        end
        
        
        % ---------- Reduced Division Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the reduced division subnetwork output neurons.
        function [ R3, self ] = compute_reduced_division_R3( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute division subnetwork output activation domain.
                [ c1, c2, R1 ] = self.unpack_reduced_absolute_division_R3_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue division subnetwork.            
                R3 = neuron_utilities.compute_reduced_absolute_division_R3( c1, c2, R1 );                                   % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R3 is a free parameter for relative division subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R3; end
                        
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the reduced division after inversion subnetwork output neurons.
        function [ R3, self ] = compute_reduced_dai_R3( self, parameters, encoding_scheme, set_flag, neuron_utilities )
            
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute division subnetwork output activation domain.
                [ c1, c2, delta1, R1 ] = self.unpack_reduced_absolute_dai_R3_parameters( parameters );
                
                % Compute the membrane capacitance for this neuron assuming that it belongs to an absolue division subnetwork.            
                R3 = neuron_utilities.compute_reduced_absolute_dai_R3( c1, c2, delta1, R1 );                                   % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R3 is a free parameter for reduced relative division subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R3; end
                        
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the absolute multiplication subnetwork neuron 3.
        function [ R3, self ] = compute_multiplication_R3( self, parameters, encoding_scheme, set_flag, neuron_utilities )
        
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute multiplication subnetwork output activation domain.                
                [ c1, c3 ] = self.unpack_absolute_multiplication_R3_parameters( parameters );
                
                % Compute the activation domain for this neuron assuming that it belongs to an absolue multiplicationb subnetwork.                            
                R3 = neuron_utilities.compute_absolute_multiplication_R3( c1, c3 );                                 % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R3 is a free parameter for relative multiplication subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R3; end
        
        end
        
            
        % Implement a function to compute the operational domain of the absolute multiplication subnetwork neuron 4.
        function [ R4, self ] = compute_multiplication_R4( self, parameters, encoding_scheme, set_flag, neuron_utilities )
        
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute multiplication subnetwork output activation domain.                
                [ c4, c5, c6, delta1, R1 ] = self.unpack_absolute_multiplication_R4_parameters( parameters );
                
                % Compute the activation domain for this neuron assuming that it belongs to an absolue multiplicationb subnetwork.                            
                R4 = neuron_utilities.compute_absolute_multiplication_R4( c4, c5, c6, delta1, R1 );                                 % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R4 is a free parameter for relative multiplication subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R4; end
        
        end
        
                
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to compute the operational domain of the reduced absolute multiplication subnetwork neuron 3.
        function [ R3, self ] = compute_reduced_multiplication_R3( self, parameters, encoding_scheme, set_flag, neuron_utilities )
        
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute multiplication subnetwork output activation domain.                
                [ c1, c2 ] = self.unpack_reduced_absolute_multiplication_R3_parameters( parameters );
                
                % Compute the activation domain for this neuron assuming that it belongs to an absolue multiplicationb subnetwork.                            
                R3 = neuron_utilities.compute_reduced_absolute_multiplication_R3( c1, c2 );                                 % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R3 is a free parameter for reduced relative multiplication subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R3; end
        
        end
        
        
        % Implement a function to compute the operational domain of the reduced absolute multiplication subnetwork neuron 4.
        function [ R4, self ] = compute_reduced_multiplication_R4( self, parameters, encoding_scheme, set_flag, neuron_utilities )
        
            % Set the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end                                            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                                                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, encoding_scheme = self.encoding_scheme_DEFAULT; end                                      % [str] Encoding Scheme (Either 'Absolute' or 'Relative'.)
            if nargin < 2, parameters = {  }; end                                                                   % [-] Parameters Cell.
            
            % Determine how to compute the membrane capacitance for this division subnetwork neuron.
            if strcmpi( encoding_scheme, 'absolute' )                                                               % If the encoding scheme is set to absolute...
                
                % Unpack the parameters required to compute the absolute multiplication subnetwork output activation domain.                
                [ c3, c4, delta1, R1 ] = self.unpack_reduced_absolute_multiplication_R4_parameters( parameters );
                
                % Compute the activation domain for this neuron assuming that it belongs to an absolue multiplicationb subnetwork.                            
                R4 = neuron_utilities.compute_reduced_absolute_multiplication_R4( c3, c4, delta1, R1 );                                 % [V] Activation Domain.
                
            elseif strcmpi( encoding_scheme, 'relative' )                                                           % If the encoding scheme is set to relative...
            
                % Throw an error.
                error( 'R4 is a free parameter for relative multiplication subnetworks.' )

            else                                                                                                    % Otherwise...

                % Throw an error.
                error( 'Invalid encoding scheme %s.  Encoding scheme must be one of: ''absolute'', ''relative''', encoding_scheme )
                
            end
            
            % Determine whether to update the neuron object.
            if set_flag, self.R = R4; end
        
        end
        
        
        %% Current Compute Functions.
        
        % Implement a function to compute the leak current associated with this neuron.
        function [ Ileak, self ] = compute_Ileak( self, U, Gm, set_flag, neuron_utilities )
            
            % Define the default input arguments.
            if nargin < 5, neuron_utilities = self.neuron_utilities; end            % [class] Neuron Utilities.
            if nargin < 4, set_flag = self.set_flag_DEFAULT; end                   	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 3, Gm = self.Gm; end                                     	% [S] Membrane Conductance.
            if nargin < 2, U = self.U; end                                         	% [V] Membrane Voltage.
            
            % Compute the leak current associated with this neuron.
            Ileak = neuron_utilities.compute_Ileak( U, Gm );                     	% [A] Leak Current.
            
            % Determine whether to update the neuron object.
            if set_flag, self.Ileak = Ileak; end
            
        end
        
        
        % Implement a function to compute the sodium channel current associated with this neuron.
        function [ Ina, self ] = compute_Ina( self, U, Gna, Am, Sm, dEm, Ah, Sh, dEh, dEna, set_flag, neuron_utilities )
            
            % Define the default input arguments.
            if nargin < 12, neuron_utilities = self.neuron_utilities; end                       	% [class] Neuron Utilities.
            if nargin < 11, set_flag = self.set_flag_DEFAULT; end                               	% [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 10, dEna = self.dEna; end                                                 	% [V] Sodium Channel Reversal Potential.
            if nargin < 9, dEh = self.dEh; end                                                      % [V] Sodium Channel Deactivation Reversal Potential.
            if nargin < 8, Sh = self.Sh; end                                                      	% [-] Sodium Channel Deactivation Slope.
            if nargin < 7, Ah = self.Ah; end                                                     	% [-] Sodium Channel Deactivation Amplitude.
            if nargin < 6, dEm = self.dEm; end                                                    	% [V] Sodium Channel Activation Reversal Potential.
            if nargin < 5, Sm = self.Sm; end                                                        % [-] Sodium Channel Activation Slope.
            if nargin < 4, Am = self.Am; end                                                     	% [-] Sodium Channel Activation Amplitude.
            if nargin < 3, Gna = self.Gna; end                                                  	% [S] Sodium Channel Conductance.
            if nargin < 2, U = self.U; end                                                       	% [V] Membrane Voltage.
            
            % Compute the sodium channel current associated with this neuron.
            Ina = neuron_utilities.compute_Ina( U, Gna, Am, Sm, dEm, Ah, Sh, dEh, dEna );          % [A] Sodium Channel Current.
            
            % Determine whether to update the neuron object.
            if set_flag, self.Ina = Ina; end
            
        end
        
        
        % Implement a function to compute the total current associated with this neuron.
        function [ Itotal, self ] = compute_Itotal( self, Ileak, Isyn, Ina, Itonic, Iapp, set_flag, neuron_utilities )
            
            % Define the default input arguments.
            if narign < 8, neuron_utilities = self.neuron_utilities; end                                % [class] Neuron Utilities.
            if nargin < 7, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if nargin < 6, Iapp = self.Iapp; end                                                    	% [A] Applied Currents
            if nargin < 5, Itonic = self.Itonic; end                                               	% [A] Tonic Current
            if nargin < 4, Ina = self.Ina; end                                                     	% [A] Sodium Channel Current
            if nargin < 3, Isyn = self.Isyn; end                                                     	% [A] Synaptic Current
            if nargin < 2, Ileak = self.Ileak; end                                                 	% [A] Leak Current
            
            % Compute the total current.
            Itotal = neuron_utilities.compute_Itotal( Ileak, Isyn, Ina, Itonic, Iapp );           % [A] Total Current
            
            % Determine whether to update the neuron object.
            if set_flag, self.Itotal = Itotal; end
            
        end
  
        
        %% Enable & Disable Functions.
        
        % Implement a function to toogle whether this neuron is enabled.
        function [ enabled_flag, self ] = toggle_enabled( self, enabled_flag, set_flag )
            
            % Set the default input arguments.
            if nargin < 3, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            if narign < 2, enabled_flag = self.enabled_flag; end                                        % [T/F] Enabled Flag.
            
            % Toggle whether the neuron is enabled.
            enabled_flag = ~enabled_flag;                                                            	% [T/F] Neuron Enabled Flag.
            
            % Determine whether to update the neuron object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        % Implement a function to enable this neuron.
        function [ enabled_flag, self ] = enable( self, set_flag )
            
            % Set the default input arguments.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Enable this neuron.
            enabled_flag = true;                                                                        % [T/F] Neuron Enabled Flag
            
            % Determine whether to update the neuron object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end
        
        
        % Implement a function to disable this neuron.
        function [ enabled_flag, self ] = disable( self, set_flag )
            
            % Set the default input arguments.
            if nargin < 2, set_flag = self.set_flag_DEFAULT; end                                        % [T/F] Set Flag (Determines whether to update the neuron object.)
            
            % Disable this neuron.
            enabled_flag = false;                                                                   	% [T/F] Neuron Enabled Flag
            
            % Determine wehther to update the neuron object.
            if set_flag, self.enabled_flag = enabled_flag; end
            
        end

        
        %% Print Functions.
        
        % Implement a function to print the properties of this neuron.
        function print( self, verbose_flag )
                
            % Set the default input arguments.
            if nargin < 2, verbose_flag = self.verbose_flag_DEFAULT; end
            
            % Print the network information.
            self.neuron_utilities.print( self.ID, self.name, self.U, self.h, self.Cm, self.Gm, self.Er, self.R, self.Am, self.Sm, self.dEm, self.Ah, self.Sh, self.dEh, self.dEna, self.tauh_max, self.tauh, self.Gna, self.minf, self.hinf, self.Ileak, self.Isyn, self.Ina, self.Itonic, self.Iapp, self.Itotal, self.enabled_flag, verbose_flag )
            
        end
        
        
        %% Save & Load Functions.
        
        % Implement a function to save neuron data as a matlab object.
        function save( self, directory, file_name, neuron )
            
            % Set the default input arguments.
            if nargin < 4, neuron = self; end
            if nargin < 3, file_name = 'Neuron.mat'; end                                                                % [-] File Name
            if nargin < 2, directory = '.'; end                                                                         % [-] Directory Path
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];                                                                  % [-] Full Directory Path
            
            % Save the neuron data.
            save( full_path, neuron )
            
        end
        
        
        % Implement a function to load neuron data as a matlab object.
        function neuron = load( ~, directory, file_name )
            
            % Set the default input arguments.
            if nargin < 3, file_name = 'Neuron.mat'; end                                                                % [-] File Name
            if nargin < 2, directory = '.'; end                                                                         % [-] Directory Name
            
            % Create the full path to the file of interest.
            full_path = [ directory, '\', file_name ];                                                                  % [-] Full Directory Path
            
            % Load the data.
            data = load( full_path );
            
            % Retrieve the desired variable from the loaded data structure.
            neuron = data.neuron;
            
        end
        
        
    end
end

