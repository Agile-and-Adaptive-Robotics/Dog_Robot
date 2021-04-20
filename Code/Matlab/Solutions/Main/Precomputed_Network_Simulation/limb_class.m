classdef limb_class

    % This class contains properties and methods related to limbs (limbs are kinematic chains comprised of links and joints).

    %% LIMB PROPERTIES
    
    % Define the class properties.
    properties
        ID
        name
        origin
        links
        num_links
        joints
        num_joints
        Ss
        Ms_joints
        Js_joints
    end
    
    %% LIMB METHODS SETUP

    % Define the class methods.
    methods
        
        % Implement the class constructor.
        function self = limb_class( ID, name, origin, links, joints )

            % Define the class properties.
            if nargin < 5, self.joints = []; else, self.joints = joints; end
            if nargin < 4, self.links = []; else, self.links = links; end
            if nargin < 3, self.origin = []; else, self.origin = origin; end
            if nargin < 2, self.name = ''; else, self.name = name; end
            if nargin < 1, self.ID = []; else, self.ID = ID; end

            % Compute the number of links.
            self.num_links = length(self.links);
            
            % Compute the number of joints.
            self.num_joints = length(self.joints);
            
            % Set the screw axes.
            self = self.set_screw_axes( );
            
        end
        
        
        % Implement a function to initialize the limb joints.
        function self = initialize_joints( self, IDs, names, parent_link_IDs, child_link_IDs, ps, Rs, vs, ws, w_screws, thetas )
        
            % Define the number of joints.
            self.num_joints = length(IDs);
            
            % Preallocate an array of joints.
            self.joints = repmat( joint_class(), 1, self.num_joints );
            
            % Create each joint object.
            for k = 1:self.num_joints               % Iterate through each of the joints...
                
                % Create this joint.
                self.joints(k) = joint_class( IDs(k), names{k}, parent_link_IDs(k), child_link_IDs(k), ps(:, k), Rs(:, :, k), vs(:, k), ws(:, k), w_screws(:, k), thetas(k) );
                
            end            
        
        end
            
        
        % Implement a function to initialize the limb links.
        function self = initialize_links( self, IDs, names, parent_joint_IDs, child_joint_IDs, start_points, end_points, lens, widths, masses, pcms, vcms, wcms, mesh_types )
            
            % Define the number of links.
            self.num_links = length(IDs);
            
            % Preallocate an array of links.
            self.links = repmat( link_class(), 1, self.num_links );

            % Create each link object.
            for k = 1:self.num_links                    % Iterate through each of the links...
                
                % Create this link.
                self.links(k) = link_class( IDs(k), names{k}, parent_joint_IDs(k), child_joint_IDs(k), start_points(:, k), end_points(:, k), lens(k), widths(k), masses(k), pcms(:, k), vcms(:, k), wcms(:, k), mesh_types{k} );
                
            end
            
        end
        
        
        % Implement a function to set the screw axes.
        function self = set_screw_axes( self, Ss )
            
            % Determine how to set the screw axes.
            if nargin < 2                                       % If a screw axis matrix was not provided...
                
                % Preallocate the screw axes.
                self.Ss = zeros(6, self.num_joints);

                % Compute each screw axis.
                for k = 1:self.num_joints                       % Iterate through each of the screw axes...

                    % Store this screw axis.
                    self.Ss(:, k) = self.joints(k).S;

                end
            
            else                                                % Otherwise...
                
                % Define the screw axes.
               self.Ss = Ss;
                
            end
                
        end
      
        
        % Implement a function to compute the home and joint assignment matrices associated with the joints.
        function [ Ms_joints, Js_joints ] = compute_joint_home_assignment_matrices( ~, Ps_home_joints, Rs_home_joints )
            
            % Set the default home joint orientation to be the identity matrix.
            if nargin < 3, Rs_home_joints = eye(3); end
            
            % Retrieve the apparent number of joints.
            n_joints = size( Ps_home_joints, 2 );
            
            % Determine whether we need to agument the size of the home joint orientation matrix.
            if size( Rs_home_joints, 3 ) == 1
                
                % Repeat the home joint orientation matrix the appropriate number of times.
                Rs_home_joints = repmat( Rs_home_joints, [ 1, 1, n_joints ] );
            
            end
                
            % Preallocate the joint home matrix.
            Ms_joints = zeros(4, 4, 1, n_joints);
            
            % Preallocate the joint assignment matrix.
            Js_joints = zeros(1, n_joints);
            
            % Build the joint home and assignment matrices.
            for k = 1:n_joints                  % Iterate through each joint...
                
                % Determine the home matrix associated with this joint.
                Ms_joints(:, :, 1, k) = [ Rs_home_joints(:, :, k), Ps_home_joints(:, k); zeros(1, 3), 1 ];             % USE THE RP TO TRANS FUNCTION INSTEAD.
                
                % Determine the joint assignment entry associated with this joint.
                Js_joints(1, k) = k;
                
            end
            
        end
        
        
        % Implement a function to set the home and joint assignment matrices associated with the limbs joints.
        function self = set_joint_home_assignment_matrices( self )
            
            % Preallocate an array to store the joint home locations.
            Ps_home_joints = zeros( 3, self.num_joints );
            
            % Preallocate an array to store the joint home orientations.
            Rs_home_joints = zeros( 3, 3, self.num_joints );
            
            % Retrieve the joint home positions and orientations.
            for k = 1:self.num_joints                   % Iterate through each joint...

                % Retrieve the home position and orientation associated with this joint.
                [ Rs_home_joints(:, :, k), Ps_home_joints(:, k) ] = TransToRp( self.joints(k).M );
                
            end            
            
            % Compute the joint home and assignment matrices.
            [ self.Ms_joints, self.Js_joints ] = self.compute_joint_home_assignment_matrices( Ps_home_joints, Rs_home_joints );
            
        end
           
        
        
        % Implement a function to compute the home and joint assignment matrices associated with the centers of mass of each joint.
        function [ Ms_coms, Js_coms ] = compute_com_home_assignment_matrices( ~ )
            
            Mcms = zeros(4, 4, 1, num_bodies);
            Jcms = zeros(1, num_bodies);
            for k = 1:num_bodies
                Mcms(:, :, 1, k) = [Rhome_cms(:, :, k), Phome_cms(:, k); zeros(1, 3), 1] ;
                Jcms(1, k) = k;
            end
            
            
            
        end
        
    end
end

