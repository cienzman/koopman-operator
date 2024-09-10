classdef PIRBFLifting < LiftingStrategy
    % PIRBFLIFTING Physics-Informed Thin-Plate Spline RBFs.
    % Enforces the theoretical constraint that the origin is an equilibrium
    % by shifting the basis functions such that phi(0) = 0.
    % This guarantees f(0,0)=0, yielding proportional error bounds around the origin.
    
    properties
        Nrbf
        cent
        phi_zero  % Stores the value of the basis functions at the origin
    end
    
    methods
        function obj = PIRBFLifting(n, Nrbf)
            obj.n = n;
            obj.Nrbf = Nrbf;
            obj.Nlift = n + Nrbf;
        end
        
        function fit(obj, Xdata)
            try
                % Use K-means to place RBF centers efficiently
                options = statset('MaxIter', 500);
                [~, C] = kmeans(Xdata', obj.Nrbf, 'Options', options);
                obj.cent = C';
            catch
                % Fallback if Statistics toolbox isn't installed
                disp('K-means unavailable. Randomly pick Nrbf data points from training set to be centers.');
                idx = randperm(size(Xdata, 2), obj.Nrbf);
                obj.cent = Xdata(:, idx);
            end
            
            % Compute the value of the RBFs at the origin x = 0
            obj.phi_zero = zeros(obj.Nrbf, 1);
            for i = 1:obj.Nrbf
                r0 = norm(obj.cent(:,i));
                if r0 > 0
                    obj.phi_zero(i) = r0^2 * log(r0);
                end
            end
        end
        
        function z = lift(obj, x)
            N = size(x, 2);
            phi = zeros(obj.Nrbf, N);
            for i = 1:obj.Nrbf
                r = sqrt(sum((x - obj.cent(:,i)).^2, 1));
                mask = r > 0;
                % Thin-plate spline: r^2 * log(r)
                phi(i, mask) = r(mask).^2 .* log(r(mask));
                
                % PHYSICS-INFORMED CONSTRAINT: Shift by phi_zero
                phi(i, :) = phi(i, :) - obj.phi_zero(i);
            end
            z = [x; phi]; 
        end
    end
end
