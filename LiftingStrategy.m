classdef LiftingStrategy < handle
    % LIFTINGSTRATEGY Abstract interface for Koopman observable functions.
    % By standardizing the interface, we can swap RBFs, Polynomials, or 
    % Neural Networks without changing the core EDMD regression code.
    
    properties
            Nlift   % Total dimension of the lifted state z
            n       % Original state dimension
        end
        methods (Abstract)
            fit(obj, Xdata)       
            z = lift(obj, x)      
        end
    end