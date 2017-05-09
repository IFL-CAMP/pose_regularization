classdef proximalMapData < handle
    
    properties (SetAccess = private)
        logHandle
        expHandle
        prodHandle
        p = 0;
    end
    
    methods 
        
        %constructor
        function pmd = proximalMapData( logHandle, expHandle, prodHandle, p )
            pmd.logHandle = logHandle;
            pmd.expHandle = expHandle;
            pmd.prodHandle = prodHandle;
            pmd.p = p;
        end
       
    end
       
    methods %(Access = {?proximalMapData})
        
        % evaluation function
        function x = eval( pmd, x, f, lambda )
        
        % compute logarithm
        cache = pmd.logHandle( x, f );
        
        % compute l2-norm of logarithm
        dist = sqrt( pmd.prodHandle( cache, cache, x ) );
        
        % calculate time for data term geodesic 
        t = geodesicLength( lambda, dist, pmd.p );
                
        % multiply computed logarithm with t and
        % compute proximal mapping for data term
        x = pmd.expHandle( x, t*cache );
        
        end % function
        
    end % methods
    
end % classdef