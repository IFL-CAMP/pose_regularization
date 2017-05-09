classdef proximalMapRegFirstOrder < handle
    
    properties (SetAccess = private)
        logHandle
        expHandle
        prodHandle
        q = 0;
    end
    
    methods
        
        %constructor
        function pmr = proximalMapRegFirstOrder( logHandle, expHandle, prodHandle, q )
            pmr.logHandle = logHandle;
            pmr.expHandle = expHandle;
            pmr.prodHandle = prodHandle;
            pmr.q = q;
        end
        
    end
        
	methods %(Access = {?proximalMapRegFirstOrder})
            
        % evaluation function
        function [x1, x2] = eval( pmr, x1, x2, alpha, lambda )
        
        % compute logarithm
        cache1 = pmr.logHandle( x1, x2 );
        cache2 = pmr.logHandle( x2, x1 );
        
        % compute l2-norm of logarithm
        dist = sqrt( pmr.prodHandle( cache1, cache1, x1 ) );
        
        % calculate time for data term geodesic 
        t = geodesicLength( lambda, dist, pmr.q, alpha );
                
        % multiply computed logarithm with t and
        % compute proximal mapping for data term
        x1 = pmr.expHandle( x1, t*cache1 );
        x2 = pmr.expHandle( x2, t*cache2 );
        
        end % function
        
    end % methods
    
end % classdef