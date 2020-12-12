classdef Cell < handle
    properties
        x         % center position
        l = 2     % length
        a = 0     % auxin
        h         % handle
        ldiv = 4  % cell length for division
    end
    
    methods
        
        % -----------------------------------------------------------------
        
        function obj = Cell(x)
            obj.x = x;
            obj.a = 1 + rand*0.01;
        end
        
        % -----------------------------------------------------------------
        
        function cells = divide(obj, ovule)
            cells = ovule.cells;
            l = obj.l; x = obj.x; a = obj.a;
            if l < obj.ldiv; return; end
            obj.l = l/2;
            obj.x = x - l/4;
            obj.a = a/2;
            
            newcel = Cell(x+l/4);
            newcel.l = l/2;
            newcel.a = a/2;
            
            cells(end+1) = newcel;
            ovule.cells = cells;
        end
        
        % -----------------------------------------------------------------
        
        function obj = grow(obj)
            dt = 1e-3; r = 1e-2;
            obj.l = obj.l + r*obj.l*dt;
        end
        
        % -----------------------------------------------------------------
        
        function obj = plot(obj)
            x = obj.x + obj.l/2 * [-1, -1, 1, 1];
            if ishandle(obj.h)
                set(obj.h, 'xdata', x, 'cdata', obj.a);
            else
                y = [-1 1 1 -1];
                obj.h = patch(x, y, obj.a); 
            end
        end
        
        % -----------------------------------------------------------------
        
    end
end
