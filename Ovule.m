classdef Ovule<handle
    properties
        cells
        x
        lp
        a
        h
        tdiv = 0;
    end
    
    methods
        % -----------------------------------------------------------------
        
        function obj = Ovule(x)
            l = diff(x);
            l = [l(1), (l(1:end-1)+l(2:end))/2, l(end)];
            for i = 1:length(x)
                cells(i) = Cell(x(i));
                cells(i).l = l(i);
            end
            
            obj.cells = cells;
        end
        
        % -----------------------------------------------------------------

        function obj = auxin(obj)
            dt = 0.5e-3;
            Da = 0.5;
            p = 1;
            Ep = 1;
            Ga = 0.1;
            
            Ra = Da/Ep/p;
            mu = (1+Ra)/2;
            L = 2*pi/acos(mu);
            
            cells = obj.cells;
            a = [cells.a];
            
            x = [cells.x];
            l = [cells.l];
            
            N = length(a);
            
            i = 1:N;
            h = circshift(i,+1);  % h = i-1
            g = circshift(i,+2);  % g = i-2
            j = circshift(i,-1);  % j = i+1
            k = circshift(i,-2);  % k = i+2
            
            pij = 2*p*a(j)./(a(j)+a(h));
            pji = 2*p*a(i)./(a(i)+a(k));
            pih = 2*p*a(h)./(a(h)+a(j));
            phi = 2*p*a(i)./(a(g)+a(i));
            
            fij = Ep*( pij.*a(i) - pji.*a(j) );
            fih = Ep*( pih.*a(i) - phi.*a(h) );
            
            da = Ga*(1-a(i))-(fij+fih) + ...
                          Da./abs(x(j)-x(i)).*(a(j)./l(j)-a(i)./l(i)) + ...
                          Da./abs(x(h)-x(i)).*(a(h)./l(h)-a(i)./l(i));
            
            a = a + da*dt;
            
            for i = 1:length(cells)
               cells(i).a = a(i); 
            end
            
        end
        % -----------------------------------------------------------------
        
        function obj = grow(obj)
            cells = obj.cells;
            for i = 1:length(cells); cells(i).grow; end
            l = [obj.cells.l];
            x = cumsum([0, l]);
            for i = 1:length(cells)
                cells(i).x = x(i)+l(i)/2;
            end
        end
        
        % -----------------------------------------------------------------
        
        function obj = divde(obj)
            cells = obj.cells;
            for i = 1:length(cells)
                cells(i).divide(obj);
            end
            cells = obj.cells;
            [~, i] = sort([cells.x]);
            obj.cells = cells(i);
        end
        
        % -----------------------------------------------------------------
        
        function plot(obj,t)
            cells = obj.cells;
            for i = 1:length(cells); cells(i).plot; end
            obj.x = [cells.x];
            obj.a = [cells.a];
            
            [pks, locs] = findpeaks(obj.a);
            
            locs = locs(pks>1.5);
            
            x = obj.x(locs);
            
            
            if length(obj.a)==100 & obj.tdiv==0
                obj.tdiv = t;
                obj.lp = obj.lp*2;
            end
            
            if length(locs)>0 & obj.tdiv==0
                obj.lp = locs;
            end

            id = zeros(size(locs));
            minv = zeros(size(locs));
            mini = zeros(size(locs));
            lp = obj.lp;    
            if obj.tdiv>0 %& length(locs)>length(lp)
                for i = 1:length(locs)
                    [minv,mini] = min(abs(locs(i)-lp));
                    if minv<3 
                        lp(mini) = locs(i);
                    else
                        id(i)=1;
                    end
                end
            end
            
%             fid = fopen(['dat/', num2str(t), '.1'],'wt');
%             fprintf(fid, '%8s\t %8s\t\n', 'h', 'x');
%             fprintf(fid, '%8.4f\t %8.4f\t\n', [ones(size(obj.x(lp)))', obj.x(lp)']');
%             fclose(fid);
%             
%             fid = fopen(['dat/', num2str(t), '.2'],'wt');
%             fprintf(fid, '%8s\t %8s\t\n', 'h', 'x');
%             fprintf(fid, '%8.4f\t %8.4f\t\n', [ones(size(x(id==1)))', x(id==1)']');
%             fclose(fid);
            
            if ~isempty(locs)
                if ishandle(obj.h)
                    set(obj.h(1), 'xdata',obj.x(lp), 'ydata', 1.5*ones(size(obj.x(lp))));
                    set(obj.h(2), 'xdata',x(id==1), 'ydata',-1.5*ones(size(x(id==1))));
                else
                    hold on
                    obj.h(1) = plot(obj.x(lp), 1.5*ones(size(obj.x(lp))), 'ob');
                    obj.h(2) = plot(NaN,NaN, 'or');
                end
            end
            
            
            
            x = [cells.x]';
            a = [cells.a]';
            l = [cells.l]';
            
            dx = 0.25;
            n = ceil((max(x)-min(x))/dx);
            
            xi = linspace(min(x),max(x),n);
            
            ai = interp1(x,a,xi,'spline');
%             fid = fopen(['dat/', num2str(t), '.xa'],'wt');
%             fprintf(fid, '%8s\t %8s\t\n', 'x', 'a');
%             fprintf(fid, '%8.4f\t %8.4f\t\n', [xi' ai']');
%             fclose(fid);
        end
        
        % -----------------------------------------------------------------
                
    end
end