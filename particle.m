% particle class
classdef particle 
    properties
        x; z; rhor;
    end
    methods
        function obj=particle(D)
            obj.x = 2000*rand();
            obj.z = -2.5*rand();
            obj.rhor = D; % Solid to Water Density Ratio
        end
    end
end