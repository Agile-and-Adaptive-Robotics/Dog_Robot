classdef neuron_utilities_class
    
    % This class contains properties and methods related to neuron utilities.
    
    
    %% NEURON UTILITIES PROPERTIES
    
    % Define the class properties.
    properties
        
        
    end
    
    
    % Define private, constant class properties.
    properties ( Access = private, Constant = true )
        
        % ---------- Neuron Properties ----------
        
        % Define default neuron parameters.
        Cm_DEFAULT = 5e-9;                                          % [C] Membrane Capacitance.
        Gm_DEFAULT = 1e-6;                                          % [S] Membrane Conductance.
        Er_DEFAULT = -60e-3;                                        % [V] Equilibrium Voltage.
        R_DEFAULT = 20e-3;                                          % [V] Activation Domain.
        Am_DEFAULT = 1;                                             % [-] Sodium Channel Activation Parameter Amplitude.
        Sm_DEFAULT = -50;                                           % [-] Sodium Channel Activation Parameter Slope.
        dEm_DEFAULT = 40e-3;                                        % [V] Sodium Channel Activation Reversal Potential.
        Ah_DEFAULT = 0.5;                                           % [-] Sodium Channel Deactivation Parameter Amplitude.
        Sh_DEFAULT = 50;                                            % [-] Sodium Channel Deactivation Parameter Slope.
        dEh_DEFAULT = 0;                                            % [V] Sodium Channel Deactivation Reversal Potential.
        dEna_DEFAULT = 110e-3;                                      % [V] Sodium Channel Reversal Potential.
        tauh_max_DEFAULT = 0.25;                                    % [s] Maximum Sodium Channel Steady State Time Constant.
        Gna_DEFAULT = 0;                                            % [S] Sodium Channel Conductance.
        Ileak_DEFAULT = 0;                                          % [A] Leak Current.
        Isyn_DEFAULT = 0;                                           % [A] Synaptic Current.
        Ina_DEFAULT = 0;                                            % [A] Sodium Channel Current.
        Itonic_DEFAULT = 0;                                         % [A] Tonic Current.
        Iapp_DEFAULT = 0;                                           % [A] Applied Current.
        Itotal_DEFAULT = 0;                                         % [A] Total Current.
                
        
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
        
        % Define the subtraction subnetwork signature.
        signature_DEFAULT = 1;                                      % [-1/+1] Subtraction Signature.
        
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
        
    end
    
    
    %% NEURON UTILITIES METHODS SETUP
    
    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = neuron_utilities_class(  )
            
            
            
        end
        
        
        %% Name Functions.
        
        % Implement a function to generate a name from an ID.
        function name = ID2name( ~, ID )
            
            % Generate a name for the neuron.
            name = sprintf( 'Neuron %0.0f', ID );
            
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
        
        
        %% Sodium Channel Functions.
        
        % Implement a function to compute the steady state sodium channel activation and deactivation parameters.
        function mhinf = compute_mhinf( self, U, Amh, Smh, dEmh )
            
            % Set the default input arguments.
            if nargin < 5, dEmh = self.dEh_DEFAULT; end                 % [V] Sodium Channel Deactivation Reversal Potential.
            if nargin < 4, Smh = self.Sh_DEFAULT; end                   % [-] Sodium Channel Deactivation Slope.
            if nargin < 3, Amh = self.Ah_DEFAULT; end                   % [-] Sodium Channel Deactivation Amplitude.
            if nargin < 2, U = 0; end                                   % [V] Membrane Voltage
            
            % Compute the steady state sodium channel activation / deactivation parameter.
            mhinf = 1./( 1 + Amh.*exp( -Smh.*( dEmh - U ) ) );          % [-] Steady State Sodium Channel Activation/Deactivation Parameter (Defaults to Deactivation Parameter)
            
        end
        
        
        % Compute the sodium channel deactivation time constant.
        function tauhs = compute_tauh( self, Us, tauh_maxs, hinfs, Ahs, Shs, dEhs )
            
            % This function computes the sodium channel deactivation time constant associated with each neuron in a network.
            
            % Inputs:
            % Us = num_neurons x 1 vector of neuron membrane voltages w.r.t. their resting potential.
            % tauh_maxs = num_neurons x 1 vector of maximum sodium channel deactivation parameter time constants.
            % hinfs = num_neurons x 1 vector of steady state sodium channel deactivation parameter values.
            % Ahs = num_neurons x 1 vector of sodium channel deactivation A parameter values.
            % Shs = num_neurons x 1 vector of sodium channel deactivation S parameter values.
            % dEhs = num_neurons x 1 vector of sodium channel deactivation parameter reversal potentials.
            
            % Outputs:
            % tauhs = num_neurons x 1 vector of sodium channel deactivation parameter time constants.
            
            % Set the default input arguments.
            if nargin < 7, dEhs = self.dEh_DEFAULT; end                                     % [V] Sodium Channel Deactivation Reversal Potential
            if nargin < 6, Shs = self.Sh_DEFAULT; end                                       % [-] Sodium Channel Deactivation Slope
            if nargin < 5, Ahs = self.Ah_DEFAULT; end                                       % [-] Sodium Channel Deactivation Amplitude
            if nargin < 4, hinfs = self.compute_mhinf( Us, Ahs, Shs, dEhs ); end            % [-] Steady State Sodium Channel Deactivation Parameter
            if nargin < 3, tauh_maxs = self.tauh_max_DEFAULT; end                           % [s] Maximum Sodium Channel Deactivation Time Constant
            
            % Compute the sodium channel deactivation time constant.
            tauhs = tauh_maxs.*hinfs.*sqrt( Ahs.*exp( -Shs.*( dEhs - Us ) ) );              % [s] Sodium Channel Deactivation Time Constant
            
        end
        
        
        % Implement a function to perform a sodium channel time constant step.
        function [ tauhs, hinfs ] = tauh_step( self, Us, tauh_maxs, Ahs, Shs, dEhs )
            
            % This function computes the sodium channel current for each neuron in a network.
            
            % Inputs:
            % Us = num_neurons x 1 vector of neuron membrane voltages w.r.t. their resting potential.
            % hs = num_neurons x 1 vector of neuron sodium channel deactivation parameters.
            % Ams = num_neurons x 1 vector of sodium channel activation A parameter values.
            % Sms = num_neurons x 1 vector of sodium channel activation S parameter values.
            % dEms = num_neurons x 1 vector of sodium channel activation parameter reversal potentials.
            % Ahs = num_neurons x 1 vector of sodium channel deactivation A parameter values.
            % Shs = num_neurons x 1 vector of sodium channel deactivation S parameter values.
            % dEhs = num_neurons x 1 vector of sodium channel deactivation parameter reversal potentials.
            % tauh_maxs = num_neurons x 1 vector of maximum sodium channel deactivation parameter time constants.
            % Gnas = num_neurons x 1 vector of sodium channel conductances for each neuron.
            % dEnas = num_neurons x 1 vector of sodium channel reversal potentials for each neuron.
            
            % Outputs:
            % Inas = num_neurons x 1 vector of sodium channel currents for each neuron.
            % minfs = num_neurons x 1 vector of neuron steady state sodium channel activation values.
            % hinfs = num_neurons x 1 vector of neuron steady state sodium channel deactivation values.
            % tauhs = num_neurons x 1 vector of sodium channel deactivation parameter time constants.
            
            % Set the default input arguments.
            if nargin < 6, dEhs = self.dEh_DEFAULT; end                                     % [V] Sodium Channel Deactivation Reversal Potential
            if nargin < 5, Shs = self.Sh_DEFAULT; end                                       % [-] Sodium Channel Deactivation Slope
            if nargin < 4, Ahs = self.Ah_DEFAULT; end                                       % [-] Sodium Channel Deactivation Amplitude
            if nargin < 3, tauh_maxs = self.tauh_max_DEFAULT; end                           % [s] Maximum Sodium Channel Deactivation Time Constant
            
            % Compute the steady state sodium channel deactivation parameters.
            hinfs = self.compute_mhinf( Us, Ahs, Shs, dEhs );                               % [-] Steady State Sodium Channel Deactivation Parameter
            
            % Compute the sodium channel deactivation time constants.
            tauhs = self.compute_tauh( Us, tauh_maxs, hinfs, Ahs, Shs, dEhs );              % [s] Sodium Channel Deactivation Time Constant
            
        end
        
        
        %% Current Functions.
        
        % Implement a function to compute leak currents.
        function Ileak = compute_Ileak( self, U, Gm )
            
            % Set the default input arguments.
            if nargin < 3, Gm = self.Gm_DEFAULT; end                % [S] Membrane Conductance.
            if nargin < 2, U = 0; end                               % [V] Membrane Voltage.
            
            % Compute the leak current.
            Ileak = -Gm.*U;                                         % [A] Leak Current.
            
        end
        
        
        % Implement a function to compute a sodium current.
        function Ina = compute_Ina( self, U, h, minf, Gna, dEna, Am, Sm, dEm )
            
            % Set the default input arguments.
            if nargin < 9, dEm = self.dEm_DEFAULT; end                              % [V] Sodium Channel Activation Reversal Potential.
            if nargin < 8, Sm = self.Sm_DEFAULT; end                                % [-] Sodium Channel Activation Slope.
            if nargin < 7, Am = self.Am_DEFAULT; end                                % [-] Sodium Channel Activation Amplitude.
            if nargin < 6, dEna = self.dEna_DEFAULT; end                            % [V] Sodium Channel Reversal Potential.
            if nargin < 5, Gna = self.Gna_DEFAULT; end                              % [S] Sodium Channel Conductance.
            if nargin < 4, minf = self.compute_mhinf( U, Am, Sm, dEm ); end         % [-] Steady State Sodium Channel Activation Parameter.
            
            % Compute the sodium current.
            Ina = Gna.*minf.*h.*( dEna - U );                                       % [A] Sodium Channel Current.
            
        end
        
        
        % Implement a function to compute sodium channel currents.
        function [ Ina, minf ] = Ina_step( self, U, h, Gna, Am, Sm, dEm, dEna )
            
            % Set the default input arguments.
            if nargin < 8, dEna = self.dEna_DEFAULT; end                            % [V] Sodium Channel Reversal Potential
            if nargin < 7, dEm = self.dEm_DEFAULT; end                              % [-] Sodium Channel Activation Reversal Potential
            if nargin < 6, Sm = self.Sm_DEFAULT; end                                % [-] Sodium Channel Activation Slope
            if nargin < 5, Am = self.Am_DEFAULT; end                                % [-] Sodium Channel Activation Amplitude
            if nargin < 4, Gna = self.Gna_DEFAULT; end                              % [S] Sodium Channel Conductance
            if nargin < 3, h = 0; end                                               % [-] Sodium Channel Deactivation Parameter
            if nargin < 2, U = 0; end                                               % [V] Membrane Voltage
            
            % Compute the steady state sodium channel activation parameter.
            minf = self.compute_mhinf( U, Am, Sm, dEm );                           % [-] Steady State Sodium Channel Activation Parameter
            
            % Compute the sodium channel current.
            Ina = self.compute_Ina( U, h, minf, Gna, dEna );                      % [A] Sodium Channel Current
            
        end
        
        
        % Implement a function to compute the total current.
        function Itotal = compute_Itotal( self, Ileak, Isyn, Ina, Itonic, Iapp )
            
            % Set the default input arguments.
            if nargin < 6, Iapp = self.Iapp_DEFAULT; end                            % [A] Applied Current.
            if nargin < 5, Itonic = self.Itonic_DEFAULT; end                        % [A] Tonic Current.
            if nargin < 4, Ina = self.Ina_DEFAULT; end                              % [A] Sodium Channel Current.
            if nargin < 3, Isyn = self.Isyn_DEFAULT; end                            % [A] Synaptic Current.
            if nargin < 2, Ileak = self.Ileak_DEFAULT; end                          % [A] Leak Current.
            
            % Compute the the total current.
            Itotal = Ileak + Isyn + Ina + Itonic + Iapp;                         	% [A] Total Current.
            
        end
        
        
        %% Neuron State Flow Functions.
        
        % Implement a function to compute the derivative of the membrane voltage with respect to time.
        function dUs = compute_dU( self, Itotals, Cms )
            
            % Set the default input arguments.
            if nargin < 3, Cms = self.Cm_DEFAULT; end               % [C] Membrane Capacitance
            if nargin < 2, Itotals = self.Itotal_DEFAULT; end       % [A] Total Current
            
            % Compute the membrane voltage derivative with respect to time.
            dUs = Itotals./Cms;                                     % [V/s] Voltage Derivative With Respect To Time
            
        end
        
        
        % Implement a function to compute the derivative of the sodium channel deactivation parameter with respect to time.
        function dhs = compute_dh( self, hs, hinfs, tauhs )
            
            % Set the default input arguments.
            if nargin < 4, tauhs = self.tauh_max_DEFAULT; end       % [s] Sodium Channel Deactivation Time Constant
            if nargin < 3, hinfs = 0; end                           % [-] Steady State Sodium Channel Deactivation Parameter
            if nargin < 2, hs = 0; end                              % [-] Sodium Channel Deactivation Parameter
            
            % Compute the sodium channel deactivation parameter derivative with respect to time.
            dhs = ( hinfs - hs )./tauhs;                            % [-/s] Sodium Channel Deactivation Parameter Derivative With Respect To Time.
            
        end
        
        
        %% Membrane Conductance Functions.
        
        % NOTE: Membrane conductance, Gm, is typically a free variable for subnetwork's who encode their output as their steady state value.  This is because the membrane conductance affects the temporal dynamics of the network, but not its steady state output.
        
        % ---------- Derivation Subnetwork Functions ----------
        
        %{
        
        % Implement a function to compute membrane conductance for a derivative subnetwork.
        function Gm = compute_derivation_Gm( ~, k, w, safety_factor )
            
            %{
            Input(s):
                k               =   [-] Derivation Subnetwork Gain.
                w               =   [-] Derivation Subnetwork Cutoff.
                safety_factor   =   [-] Derivation Subnetwork Safety Factor.
            
            Output(s):
                Gm              =   [S] Membrane Conductance.
            %}
            
            % Set the default input arugments.
            if nargin < 4, safety_factor = self.sf_derivation_DEFAULT; end
            if nargin < 3, w = self.w_derivation_DEFAULT; end
            if nargin < 2, k = self.c_derivation_DEFAULT; end
            
            % Compute the required membrance conductance.
            Gm = ( 1 - safety_factor )./( k.*w );    
            
        end
        
        %}
        
        
        %% Membrane Capacitance Functions.
        
        % NOTE: Membrane capacitance, Cm, is typically a free variable for subnetwork's who encode their output as their steady state value.  This is because the membrane capacitance affects the temporal dynamics of the network, but not its steady state output.

        
        % ---------- Derivation Subnetwork Functions ----------
        
        %{
        
        % Implement a function to compute membrane capacitances for a derivative subnetwork.
        function [ Cm1, Cm2 ] = compute_derivation_Cms( ~, Gm, k, w )
            
            %{
            Input(s):
                Gm = [S] Membrane Conductance.
                k = [-] Derivation Subnetwork Gain.
            `   w = [-] Derivation Subnetwork Cutoff.
            
            Output(s):
                Cm1 = [F] Membrane Conductance (Neuron 1).
                Cm2 = [F] Membrane Conductnace (Neuron 2).
            %}
            
            % Set the default input arugments.
            if nargin < 4, w = self.w_derivation_DEFAULT; end
            if nargin < 3, k = self.c_derivation_DEFAULT; end
            if nargin < 2, Gm = 1e-6; end
            
           % Compute the required time constant.
            tau = 1./w;
            
            % Compute the required membrane capacitance of the second neuron.
            Cm2 = Gm.*tau;
            
            % Compute the required membrane capacitance of the first neuron.
            Cm1 = Cm2 - ( Gm.^2 ).*k; 
            
        end
        
        %}
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        %{
        
        % Implement a function to compute the membrane capacitances for an integration subnetwork.
        function Cm = compute_integration_Cm( ~, ki_mean )
        
            %{
            Input(s):
                ki_mean     =   [-] Intergration Subnetwork Gain.
            
            Output(s):
                Cm          =   [F] Membrane Conductance.
            %}
            
            % Set the default input arguments.
            if nargin < 2, ki_mean = self.c_integration_mean_DEFAULT; end
            
            % Compute the integration subnetwork membrane capacitance.
            Cm = 1./( 2*ki_mean );
            
        end
        
        %}
        
        
        %% Sodium Channel Conductance Functions.
        
        % ---------- Transmission Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of an absolute transmission subnetwork neuron.
        function Gna = compute_absolute_transmission_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end

        
        % Implement a function to compute the sodium channel conductance of a relative transmission subnetwork neuron.
        function Gna = compute_relative_transmission_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function compute the sodium channel conductance of a transmission subnetwork neuron.
        function Gna = compute_transmission_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_absolute_transmission_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_relative_transmission_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Addition Subnetwork Functions ----------

        % Implement a function to compute the sodium channel conductance of an absolute addition subnetwork neuron.
        function Gna = compute_absolute_addition_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a relative addition subnetwork neuron.
        function Gna = compute_relative_addition_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function compute the sodium channel conductance of a addition subnetwork neuron.
        function Gna = compute_addition_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_absolute_addition_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_relative_addition_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of an absolute subtraction subnetwork neuron.
        function Gna = compute_absolute_subtraction_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to comptue the sodium channel conductance of a relative subtraction subnetwork neuron.
        function Gna = compute_relative_subtraction_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        

        % Implement a function to comute the sodium channel conductance of a subtraction subnetwork neuron.
        function Gna = compute_subtraction_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_absolute_subtraction_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_relative_subtraction_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Double Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of an absolute double subtraction subnetwork neuron.
        function Gna = compute_absolute_double_subtraction_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a relative double subtraction subnetwork neuron.
        function Gna = compute_relative_double_subtraction_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to comute the sodium channel conductance of a double subtraction subnetwork neuron.
        function Gna = compute_double_subtraction_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_absolute_double_subtraction_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_relative_double_subtraction_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of an absolute inversion subnetwork neuron.
        function Gna = compute_absolute_inversion_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a relative inversion subnetwork neuron.
        function Gna = compute_relative_inversion_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to comute the sodium channel conductance of an inversion subnetwork neuron.
        function Gna = compute_inversion_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_absolute_inversion_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_relative_inversion_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------

        % Implement a function to compute the sodium channel conductance of a reduced absolute inversion subnetwork neuron.
        function Gna = compute_reduced_absolute_inversion_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a reduced relative inversion subnetwork neuron.
        function Gna = compute_reduced_relative_inversion_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to comute the sodium channel conductance of a reduced inversion subnetwork neuron.
        function Gna = compute_reduced_inversion_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_reduced_absolute_inversion_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_reduced_relative_inversion_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division Subnetwork Functions ----------

        % Implement a function to compute the sodium channel conductance of an absolute division subnetwork neuron.
        function Gna = compute_absolute_division_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a relative division subnetwork neuron.
        function Gna = compute_relative_division_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to comute the sodium channel conductance of a division subnetwork neuron.
        function Gna = compute_division_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_absolute_division_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_relative_division_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division Subnetwork Functions ----------

        % Implement a function to compute the sodium channel conductance of a reduced absolute division subnetwork neuron.
        function Gna = compute_reduced_absolute_division_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a reduced relative division subnetwork neuron.
        function Gna = compute_reduced_relative_division_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a reduced division subnetwork neuron.
        function Gna = compute_reduced_division_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_reduced_absolute_division_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_reduced_relative_division_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Division After Inversion Subnetwork Functions ----------

        % Implement a function to compute the sodium channel conductance of an absolute division after inversion subnetwork neuron.
        function Gna = compute_absolute_dai_Gna( ~ )
        
            % Compute the sodium channel conductance.
            Gna = 0;
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a relative division after inversion subnetwork neuron.
        function Gna = compute_relative_dai_Gna( ~ )
        
            % Compute the sodium channel conductance.
            Gna = 0;
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a division after inversion subnetwork neuron.
        function Gna = compute_dai_Gna( ~ )
           
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_reduced_absolute_dai_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_reduced_relative_dai_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------

        % Implement a function to compute the sodium channel conductance of a reduced absolute division after inversion subnetwork neuron.
        function Gna = compute_reduced_absolute_dai_Gna( ~ )
        
            % Compute the sodium channel conductance.
            Gna = 0;
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a reduced relative division after inversion subnetwork neuron.
        function Gna = compute_reduced_relative_dai_Gna( ~ )
        
            % Compute the sodium channel conductance.
            Gna = 0;
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a reduced division after inversion subnetwork neuron.
        function Gna = compute_reduced_dai_Gna( ~ )
           
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                Gna = self.compute_reduced_absolute_dai_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                Gna = self.compute_reduced_relative_dai_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of an absolute multiplication subnetwork neuron.
        function [ Gna3, Gna4 ] = compute_absolute_multiplication_Gna( self )
            
            % Compute the absolute inversion subnetwork sodium channel conductance.
            Gna3 = compute_absolute_inversion_Gna(  );                                      % [S] Sodium Channel Conductance.
            
            % Compute the absolute division subnetwork sodium channel conductance.
            Gna4 = self.compute_absolute_dai_Gna(  );                  % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a relative multiplication subnetwork neuron.
        function [ Gna3, Gna4 ] = compute_relative_multiplication_Gna( ~ )
            
            % Compute the relative inversion subnetwork sodium channel conductance.
            Gna3 = compute_relative_inversion_Gna(  );                                      % [S] Sodium Channel Conductance.
            
            % Compute the relative division subnetwork sodium channel conductance.
            Gna4 = self.compute_relative_dai_Gna(  );                  % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to comute the sodium channel conductance of a multiplication subnetwork neuron.
        function [ Gna3, Gna4 ] = compute_multiplication_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                [ Gna3, Gna4 ] = self.compute_absolute_multiplication_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                [ Gna3, Gna4 ] = self.compute_relative_multiplication_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of a reduced absolute multiplication subnetwork neuron.
        function [ Gna3, Gna4 ] = compute_reduced_absolute_multiplication_Gna( ~ )
            
            % Compute the absolute inversion subnetwork sodium channel conductance.
            Gna3 = compute_reduced_absolute_inversion_Gna(  );                                      % [S] Sodium Channel Conductance.
            
            % Compute the absolute division subnetwork sodium channel conductance.
            Gna4 = self.compute_reduced_absolute_dai_Gna(  );                  % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a reduced relative multiplication subnetwork neuron.
        function [ Gna3, Gna4 ] = compute_reduced_relative_multiplication_Gna( ~ )
            
            % Compute the relative inversion subnetwork sodium channel conductance.
            Gna3 = compute_reduced_relative_inversion_Gna(  );                                      % [S] Sodium Channel Conductance.
            
            % Compute the relative division subnetwork sodium channel conductance.
            Gna4 = self.compute_reduced_relative_dai_Gna(  );                  % [S] Sodium Channel Conductance.
            
        end
        
        
        % Implement a function to comute the sodium channel conductance of a reduced multiplication subnetwork neuron.
        function [ Gna3, Gna4 ] = compute_reduced_multiplication_Gna( self, encoding_scheme )
            
            % Set the default input arguments.
            if nargin < 2, encoding_scheme = self.encoding_scheme_DEFAULT; end
            
            % Determine how to compute the sodium channel conductance.
            if strcmpi( encoding_scheme, 'absolute' )
               
                % Compute the sodium channel conductance using an absolute encoding scheme.
                [ Gna3, Gna4 ] = self.compute_reduced_absolute_multiplication_Gna(  );
                
            elseif strcmpi( encoding_scheme, 'relative' )
                
                % Compute the sodium channel conductance using a relative encoding scheme.
                [ Gna3, Gna4 ] = self.compute_reduced_relative_multiplication_Gna(  );
            
            else
            
                % Throw an error.
                error( 'Invalid encoding scheme.  Must be either: ''absolute'' or ''relative''.' )
                
            end
            
        end
        
        
        % ---------- Derivation Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of a derivation subnetwork neuron.
        function Gna = compute_derivation_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        % ---------- Integration Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductance of a integration subnetwork neuron.
        function Gna = compute_integration_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        

        % Implement a function to compute the sodium channel conductance of a voltage based integration subnetwork neuron.
        function Gna = compute_vbi_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end


        % Implement a function to compute the sodium channel conductance of a split voltage based integration subnetwork neuron.
        function Gna = compute_svbi_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
       
        % ---------- Central Pattern Generator Subnetwork Functions ----------
        
        % Implement a function to compute the sodium channel conductances for a CPG subnetwork.
        function Gna = compute_cpg_Gna( self, R, Gm, Am, Sm, dEm, Ah, Sh, dEh, dEna )
            
            % Set the default input arguments.
            if nargin < 10, dEna = self.dEna_DEFAULT; end                   % [V] Sodium Channel Reversal Potential
            if nargin < 9, dEh = self.dEh_DEFAULT; end                      % [V] Sodium Channel Deactivation Reversal Potential
            if nargin < 8, Sh = self.Sh_DEFAULT; end                        % [-] Sodium Channel Deactivation Slope
            if nargin < 7, Ah = self.Ah_DEFAULT; end                        % [-] Sodium Channel Deactivation Amplitude
            if nargin < 6, dEm = self.dEm_DEFAULT; end                      % [-] Sodium Channel Activation Reversal Potential
            if nargin < 5, Sm = self.Sm_DEFAULT; end                        % [-] Sodium Channel Activation Slope
            if nargin < 4, Am = self.Am_DEFAULT; end                        % [-] Sodium Channel Activation Amplitude
            if nargin < 3, Gm = self.Gm_DEFAULT; end                        % [S] Membrane Conductance
            if nargin < 2, R = self.R_DEFAULT; end                          % [V] Activation Domain
            
            % Compute the steady state sodium channel activation & devactivation parameters at the upper equilibrium.
            minf_upper = self.compute_mhinf( R, Am, Sm, dEm );              % [-] Steady State Sodium Channel Activation Parameter
            hinf_upper = self.compute_mhinf( R, Ah, Sh, dEh );              % [-] Steady State Sodium Channel Deactivation Parameter.
            
            % Compute the sodium channel conductance for each half-center neuron.
            Gna = ( Gm.*R )./( minf_upper.*hinf_upper.*( dEna - R ) );      % [S] Sodium Channel Conductance.  Equation straight from Szczecinski's CPG example.
            
        end
        
        
        % Implement a function to compute the sodium channel conductance of a driven multistate cpg subnetwork neuron.
        function Gna = compute_dmcpg_Gna( ~ )
            
            % Compute the sodium channel conductance.
            Gna = 0;              % [S] Sodium Channel Conductance.
            
        end
        
        
        %% Activation Domain Functions.
        
        % ---------- Transmission Subnetwork Functions ----------

        % Implement a function to compute the operational domain of the absolute transmission subnetwork output neuron.
        function R2 = compute_absolute_transmission_R2( self, c, R1 )
        
            % Set the default input arguments.
            if nargin < 3, R1 = self.R_DEFAULT; end
            if nargin < 2, c = self.c_absolute_transmission_DEFAULT; end
            
            % Compute the operational domain.
            R2 = c*R1;
            
        end

        
        % ---------- Addition Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the absolute addition subnetwork output neurons.
        function Rn = compute_absolute_addition_Rn( self, cs, Rs_input )
            
            % Set the default input arguments.
            if nargin < 3, Rs_input = self.R_DEFAULT; end                        % [V] Activation Domain
            if nargin < 2, cs = self.c_absolute_addition_DEFAULT*ones( 1, length( Rs_input ) ); end
            
            % Compute the operational domain.
            Rn = sum( cs.*Rs_input );                                              % [V] Activation Domain
            
        end
        

        % ---------- Subtraction Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the absolute subtraction subnetwork output neurons.
        function Rn = compute_absolute_subtraction_Rn( self, cs, s_ks, Rs_input )
            
            % Set the default input arguments.
            if nargin < 4, Rs_input = self.R_DEFAULT; end                                                           % [V] Activation Domain.
            if nargin < 3, s_ks = self.signature_DEFAULT*ones( 1, length( Rs_input ) ); end                       	% [-1, +1] Input Signature.
            if nargin < 2, cs = self.c_absolute_subtraction_DEFAULT*ones( 1, length( Rs_input ) ); end
            
            % Compute the excitatory and inhibitory input indexes.
            excitatory_indexes = s_ks == 1;
            inhibitory_indexes = s_ks == -1;
            
            % Compute the excitatory and inhibitory gains.
            cs_excitatory = cs( excitatory_indexes );
            cs_inhibitory = cs( inhibitory_indexes );
            
            % Retrieve the excitatory and inhibitory domains.       
            R_excitatory = sum( Rs_input( excitatory_indexes ) );           % [V] Excitatory Activation Domain.
            R_inhibitory = sum( Rs_input( inhibitory_indexes ) );          	% [V] Inhibitory Activation Domain.
            
            % Compute the operational domain.
            Rn = max( cs_excitatory.*R_excitatory, cs_inhibitory.*R_inhibitory );                       	% [V] Activation Domain.
            
        end
        
        
        % ---------- Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the absolute inversion subnetwork output neuron.
        function R2 = compute_absolute_inversion_R2( self, c1, c3 )

            % Set the default input arguments.
            if nargin < 3, c3 = self.c3_absolute_inversion_DEFAULT; end        	% [-] Inversion Subnetwork Gain 3.
            if nargin < 2, c1 = self.c1_absolute_inversion_DEFAULT; end          % [-] Inversion Subnetwork Gain 1.

            % Compute the operational domain.
            R2 = c1/c3;                                                  % [V] Activation Domain.

        end
        
        
        % Implement a function to compute the operational domain of the relative inversion subnetwork output neuron.
        function R2 = compute_relative_inversion_R2( self, R2 )

            % Set the default input arguments.
            if nargin < 2, R2 = self.R_relative_inversion_DEFAULT; end          % [V] Maximum Membrane Voltage.

        end
        
        
        % ---------- Reduced Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the reduced absolute inversion subnetwork output neuron.
        function R2 = compute_reduced_absolute_inversion_R2( self, c1, c2 )
        
            % Set the default input arguments.
            if nargin < 3, c2 = self.c2_reduced_absolute_inversion_DEFAULT; end        	% [-] Inversion Subnetwork Gain 2.
            if nargin < 2, c1 = self.c1_reduced_absolute_inversion_DEFAULT; end          % [-] Inversion Subnetwork Gain 1.
            
            % Compute the operational domain.
            R2 = c1/c2;
            
        end
        
        
        % ---------- Division Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the absolute division subnetwork output neuron.
        function R3 = compute_absolute_division_R3( self, c1, c3, R1 )
            
            % Set the default input arguments.
            if nargin < 4, R1 = self.R_DEFAULT; end
            if nargin < 3, c3 = self.c3_absolute_division_DEFAULT; end
            if nargin < 2, c1 = self.c1_absolute_division_DEFAULT; end
            
            % Compute the operational domain.
            R3 = ( c1*R1 )/c3;

        end

                
        % ---------- Reduced Division Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the reduced absolute division subnetwork output neuron.
        function R3 = compute_reduced_absolute_division_R3( self, c1, c2, R1 )
        
            % Set the default input arguments.
            if nargin < 4, R1 = self.R_DEFAULT; end
            if nargin < 3, c2 = self.c2_reduced_absolute_division_DEFAULT; end
            if nargin < 2, c1 = self.c1_reduced_absolute_division_DEFAULT; end
            
            % Compute the operational domain.
            R3 = ( c1*R1 )/c2;
            
        end

        
        % ---------- Division After Inversion Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the absolute division after inversion subnetwork output neuron.
        function R3 = compute_absolute_dai_R3( self, c1, c2, c3, delta1, R1 )
            
            % Set the default input arguments.
            if nargin < 5, R1 = self.R_DEFAULT; end
            if nargin < 4, c3 = self.c3_absolute_dai_DEFAULT; end
            if nargin < 3, c2 = self.c2_absolute_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_absolute_dai_DEFAULT; end
            
            % Compute the operational domain.
            R3 = ( c1*R1 )/( c2*delta1 + c3 );

        end
        
        
        % ---------- Reduced Division After Inversion Subnetwork Functions ----------

        % Implement a function to compute the operational domain of the reduced absolute division after inversion subnetwork output neuron.
        function R3 = compute_reduced_absolute_dai_R3( self, c1, c2, delta1, R1 )
            
            % Set the default input arguments.
            if nargin < 5, R1 = self.R_DEFAULT; end
            if nargin < 4, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
            if nargin < 3, c2 = self.c2_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, c1 = self.c1_reduced_absolute_dai_DEFAULT; end

            % Compute the opertional domain.
            R3 = ( c1*R1 )/( delta1 + c2 );
            
        end
        
        
        % ---------- Multiplication Subnetwork Functions ----------
        
        % Implement a function to compute the operational domain of the absolute multiplication subnetwork neuron 3.
        function R3 = compute_absolute_multiplication_R3( self, c1, c3 )
            
            % Set the default input arguments.
            if nargin < 3, c3 = self.c3_absolute_inversion_DEFAULT; end
            if nargin < 2, c1 = self.c1_absolute_inversion_DEFAULT; end
            
            % Compute the operational domain of the absolute inversion subnetwork.
            R3 = self.compute_absolute_inversion_R2( c1, c3 );
            
        end
        
        
        % Implement a function to compute the operational domain of the absolute multiplication subnetwork neuron 4.
        function R4 = compute_absolute_multiplication_R4( self, c4, c5, c6, delta1, R1 )
        
            % Set the default input arguments.
            if nargin < 6, R1 = self.R_DEFAULT; end
            if nargin < 5, delta1 = self.delta_absolute_inversion_DEFAULT; end
            if nargin < 4, c6 = self.c3_absolute_division_DEFAULT; end
            if nargin < 3, c5 = self.c2_absolute_division_DEFAULT; end
            if nargin < 2, c4 = self.c1_absolute_division_DEFAULT; end
            
            % Compute the operational domain of the absolute division subnetwork.
            R4 = self.compute_absolute_dai_R3( c4, c5, c6, delta1, R1 );
                
        end
        
            
        % Implement a function to compute the operational domain of the absolute multiplication subnetwork neurons.
        function [ R3, R4 ] = compute_absolute_multiplication_Rs( self, c1, c3, c4, c5, c6, delta1, R1 )
        
            % Set the default input arguments.
            if nargin < 7, R1 = self.R_DEFAULT; end
            if nargin < 6, c6 = self.c3_absolute_division_DEFAULT; end
            if nargin < 5, c5 = self.c2_absolute_division_DEFAULT; end
            if nargin < 4, c4 = self.c1_absolute_division_DEFAULT; end
            if nargin < 3, c3 = self.c3_absolute_inversion_DEFAULT; end
            if nargin < 2, c1 = self.c1_absolute_inversion_DEFAULT; end
            
            % Compute the operational domain of the absolute inversion subnetwork.
            R3 = self.compute_absolute_multiplication_R3( c1, c3 );
            
            % Compute the operational domain of the absolute division subnetwork.
            R4 = self.compute_absolute_multiplication_R4( c4, c5, c6, delta1, R1 );
            
        end
        
        
        % ---------- Reduced Multiplication Subnetwork Functions ----------

        % Implement a function to compute the operational domain of the reduced absolute multiplication subnetwork neuron 3.
        function R3 = compute_reduced_absolute_multiplication_R3( self, c1, c2 )
            
            % Set the default input arguments.
            if nargin < 3, c2 = self.c2_reduced_absolute_inversion_DEFAULT; end
            if nargin < 2, c1 = self.c1_reduced_absolute_inversion_DEFAULT; end
            
            % Compute the operational domain of the absolute inversion subnetwork.
            R3 = self.compute_reduced_absolute_inversion_R2( c1, c2 );
            
        end
        
        
        % Implement a function to compute the operational domain of the reduced absolute multiplication subnetwork neuron 4.
        function R4 = compute_reduced_absolute_multiplication_R4( self, c3, c4, delta1, R1 )
        
            % Set the default input arguments.
            if nargin < 5, R1 = self.R_DEFAULT; end
            if nargin < 4, delta1 = self.delta_reduced_absolute_inversion_DEFAULT; end
            if nargin < 3, c4 = self.c2_reduced_absolute_dai_DEFAULT; end
            if nargin < 2, c3 = self.c1_reduced_absolute_dai_DEFAULT; end
            
            % Compute the operational domain of the absolute division subnetwork.
            R4 = self.compute_reduced_absolute_dai_R3( c3, c4, delta1, R1 );

        end
        
        
        % Implement a function to compute the operational domain of the reduced multiplication subnetwork neurons.
        function [ R3, R4 ] = compute_reduced_absolute_multiplication_Rs( self, c1, c2, c3, c4, delta1, R1 )
        
            % Set the default input arguments.
            if nargin < 7, R1 = self.R_DEFAULT; end
            if nargin < 6, delta1 = self.delta_absolute_inversion_DEFAULT; end
            if nargin < 5, c4 = self.c2_reduced_absolute_dai_DEFAULT; end
            if nargin < 4, c3 = self.c1_reduced_absolute_dai_DEFAULT; end
            if nargin < 3, c2 = self.c2_reduced_absolute_inversion_DEFAULT; end
            if anrgin < 2, c1 = self.c1_reduced_absolute_inversion_DEFAULT; end
            
            % Compute the reduced absolute inversion output operational domain.
            R3 = self.compute_reduced_absolute_multiplication_R3( c1, c2 );
            
            % Compute the reduced absolute division output operational domain.
            R4 = self.compute_reduced_absolute_multiplication_R4( c3, c4, delta1, R1 );
            
        end
        
        
        %% Print Functions.
        
        % Implement a function to print the properties of this neuron.
        function print( ~, ID, name, U, h, Cm, Gm, Er, R, Am, Sm, dEm, Ah, Sh, dEh, dEna, tauh_max, tauh, Gna, minf, hinf, Ileak, Isyn, Ina, Itonic, Iapp, Itotal, enabled_flag, verbose_flag )
        
            % Define the default input arguments.
            if nargin < 29, verbose_flag = false; end
            
            % Print out a header for this neuron.
            fprintf( '----------------------------------- NEURON %0.0f: %s -----------------------------------\n', ID, name )
            
            % Determine which information to print about this neuron.
            if verbose_flag         % If we want to print all of the information...
                
                fprintf( 'Membrane Voltage: \t\t\t\t\t\t\t\t\t\tU \t\t\t= \t%0.2f \t\t[mV]\n', U*( 10^3 ) )
                fprintf( 'Sodium Channel Deactivation: \t\t\t\t\t\t\th \t\t\t= \t%0.2f \t\t[-]\n', h )

                fprintf( 'Membrane Capacitance: \t\t\t\t\t\t\t\t\tCm \t\t\t= \t%0.2f \t\t[nF]\n', Cm*( 10^9 ) )
                fprintf( 'Membrane Conductance: \t\t\t\t\t\t\t\t\tGm \t\t\t= \t%0.2f \t\t[muS]\n', Gm*( 10^6 ) )
                fprintf( 'Resting Membrane Voltage: \t\t\t\t\t\t\t\tEr \t\t\t= \t%0.2f \t\t[mV]\n', Er*( 10^3 ) )
                fprintf( 'Maximum Membrane Voltage: \t\t\t\t\t\t\t\tR \t\t\t= \t%0.2f \t\t[mV]\n', R*( 10^3 ) )
                
                fprintf( 'Sodium Channel Activation Amplitude: \t\t\t\t\tAm \t\t\t= \t%0.2f \t[mV]\n', Am*( 10^3 ) )
                fprintf( 'Sodium Channel Activation Slope: \t\t\t\t\t\tSm \t\t\t= \t%0.2f \t[mV]\n', Sm*( 10^3 ) )
                fprintf( 'Sodium Channel Activation Reversal Potential: \t\t\tdEm \t\t= \t%0.2f \t\t[mV]\n', dEm*( 10^3 ) )
                
                fprintf( 'Sodium Channel Deactivation Amplitude: \t\t\t\t\tAh \t\t\t= \t%0.2f \t\t[mV]\n', Ah*( 10^3 ) )
                fprintf( 'Sodium Channel Deactivation Slope: \t\t\t\t\t\tSh \t\t\t= \t%0.2f \t[mV]\n', Sh*( 10^3 ) )
                fprintf( 'Sodium Channel Deactivation Reversal Potential: \t\tdEh \t\t= \t%0.2f \t\t[mV]\n', dEh*( 10^3 ) )

                fprintf( 'Sodium Channel Reversal Potential: \t\t\t\t\t\tdEna \t\t= \t%0.2f \t\t[mV]\n', dEna*( 10^3 ) )
                fprintf( 'Maximum Sodium Channel Deactivation Time Constant: \t\ttauh_max \t= \t%0.2f \t\t[ms]\n', tauh_max*( 10^3 ) )
                fprintf( 'Sodium Channel Deactivation Time Constant: \t\t\t\ttauh \t\t= \t%0.2f \t\t[ms]\n', tauh*( 10^3 ) )
                fprintf( 'Sodium Channel Conductance: \t\t\t\t\t\t\tGna \t\t= \t%0.2f \t\t[muS]\n', Gna*( 10^6 ) )
                
                fprintf( 'Steady State Sodium Channel Activation Parameter: \t\tminf \t\t= \t%0.2f \t\t[-]\n', minf )
                fprintf( 'Steady State Sodium Channel Deactivation Parameter: \thinf \t\t= \t%0.2f \t\t[-]\n', hinf )

                fprintf( 'Leak Current: \t\t\t\t\t\t\t\t\t\t\tIleak \t\t= \t%0.2f \t\t[nA]\n', Ileak*( 10^9 ) )
                fprintf( 'Synaptic Current: \t\t\t\t\t\t\t\t\t\tIsyn \t\t= \t%0.2f \t\t[nA]\n', Isyn*( 10^9 ) )
                fprintf( 'Sodium Channel Current: \t\t\t\t\t\t\t\tIna \t\t= \t%0.2f \t\t[nA]\n', Ina*( 10^9 ) )
                fprintf( 'Tonic Current: \t\t\t\t\t\t\t\t\t\t\tItonic \t\t= \t%0.2f \t\t[nA]\n', Itonic*( 10^9 ) )
                fprintf( 'Applied Current: \t\t\t\t\t\t\t\t\t\tIapp \t\t= \t%0.2f \t\t[nA]\n', Iapp*( 10^9 ) )
                fprintf( 'Total Current: \t\t\t\t\t\t\t\t\t\t\tItotal \t\t= \t%0.2f \t\t[nA]\n', Itotal*( 10^9 ) )

                fprintf( 'Enabled Flag: \t\t\t\t\t\t\t\t\t\t\tEnabled \t= \t%0.2f \t\t[T/F]\n', enabled_flag )

            else                    % Otherwise...
                
                fprintf( 'Membrane Capacitance: \t\t\t\t\t\t\t\t\tCm \t\t\t= \t%0.2f \t\t[nF]\n', Cm*( 10^9 ) )
                fprintf( 'Membrane Conductance: \t\t\t\t\t\t\t\t\tGm \t\t\t= \t%0.2f \t\t[muS]\n', Gm*( 10^6 ) )
                fprintf( 'Resting Membrane Voltage: \t\t\t\t\t\t\t\tEr \t\t\t= \t%0.2f \t\t[mV]\n', Er*( 10^3 ) )
                fprintf( 'Maximum Membrane Voltage: \t\t\t\t\t\t\t\tR \t\t\t= \t%0.2f \t\t[mV]\n', R*( 10^3 ) )
                
                fprintf( 'Sodium Channel Activation Amplitude: \t\t\t\t\tAm \t\t\t= \t%0.2f \t[mV]\n', Am*( 10^3 ) )
                fprintf( 'Sodium Channel Activation Slope: \t\t\t\t\t\tSm \t\t\t= \t%0.2f \t[mV]\n', Sm*( 10^3 ) )
                fprintf( 'Sodium Channel Activation Reversal Potential: \t\t\tdEm \t\t= \t%0.2f \t\t[mV]\n', dEm*( 10^3 ) )
                
                fprintf( 'Sodium Channel Deactivation Amplitude: \t\t\t\t\tAh \t\t\t= \t%0.2f \t\t[mV]\n', Ah*( 10^3 ) )
                fprintf( 'Sodium Channel Deactivation Slope: \t\t\t\t\t\tSh \t\t\t= \t%0.2f \t[mV]\n', Sh*( 10^3 ) )
                fprintf( 'Sodium Channel Deactivation Reversal Potential: \t\tdEh \t\t= \t%0.2f \t\t[mV]\n', dEh*( 10^3 ) )

                fprintf( 'Sodium Channel Reversal Potential: \t\t\t\t\t\tdEna \t\t= \t%0.2f \t\t[mV]\n', dEna*( 10^3 ) )
                fprintf( 'Maximum Sodium Channel Deactivation Time Constant: \t\ttauh_max \t= \t%0.2f \t\t[ms]\n', tauh_max*( 10^3 ) )
                fprintf( 'Sodium Channel Conductance: \t\t\t\t\t\t\tGna \t\t= \t%0.2f \t\t[muS]\n', Gna*( 10^6 ) )
                
            end
            
            % Print out a footer.
            fprintf( '-------------------------------------------------------------------------------------------\n\n' )
            
        end
        
    end
    
    
end