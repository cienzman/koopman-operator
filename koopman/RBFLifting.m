classdef RBFLifting < LiftingStrategy
% RBFLIFTING Thin-Plate Spline RBFs. 
    % Uses K-Means clustering on training data for optimal center placement,
    % which prevents dead/useless RBFs in unvisited regions of the state space.

    properties
        Nrbf
        cent
    end
    
    methods
        function obj = RBFLifting(n, Nrbf)
            obj.n = n;
            obj.Nrbf = Nrbf;
            obj.Nlift = n + Nrbf;
        end
        
        function fit(obj, Xdata)
            try
                % Use K-means to place RBF centers exactly where data lives
                options = statset('MaxIter', 500);
                [~, C] = kmeans(Xdata', obj.Nrbf, 'Options', options);
                obj.cent = C';
            catch
                % Fallback if Statistics toolbox isn't installed
                disp('K-means unavailable. Randomly pick Nrbf data points from the actual training set to be centers.');
                idx = randperm(size(Xdata, 2), obj.Nrbf);
                obj.cent = Xdata(:, idx)
            end
        end
        
        function z = lift(obj, x)
            N = size(x, 2);
            phi = zeros(obj.Nrbf, N);
            for i = 1:obj.Nrbf
                r = sqrt(sum((x - obj.cent(:,i)).^2, 1));
                mask = r > 0;
                phi(i, mask) = r(mask).^2 .* log(r(mask));
            end
            z = [x; phi]; 
        end
    end
end 



    