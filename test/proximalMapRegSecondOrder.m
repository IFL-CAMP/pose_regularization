classdef proximalMapRegSecondOrder < handle
    
    properties (SetAccess = private)
        logHandle    % logarithm map
        expHandle    % exponential map
        prodHandle   % inner product 
        frameHandle  % moving frame of geodesic
        parHandle    % parallel transport
        adjustHandle % adjust coefficients
        r = 1;       % exponent for regularization
        steps = 1;  % number of iteration steps
    end
    
    methods
        
        %constructor
        function pmr = proximalMapRegSecondOrder( logHandle, ...
                                                  expHandle, ...
                                                  prodHandle, ...
                                                  frameHandle, ...
                                                  parHandle, ...
                                                  adjustHandle, ...
                                                  r, steps )
            pmr.logHandle = logHandle;
            pmr.expHandle = expHandle;
            pmr.prodHandle = prodHandle;
            pmr.frameHandle = frameHandle;
            pmr.parHandle = parHandle;
            pmr.adjustHandle = adjustHandle;
            pmr.r = r;
            pmr.steps = steps;
        end
        
    end
        
	methods %(Access = {?proximalMapRegSecondOrder})
            
        % evaluation function
        function [y1, y2, y3] = eval( pmr, x1, x2, x3, beta, lambda )
        
            % initialization
            y1 = x1;
            y2 = x2;
            y3 = x3;
            % initialize y1 y2 y3 as points on the line passing tru the
            % intersection of medians
            M = pmr.expHandle(y1, 0.5 * pmr.logHandle( y1, y3 ));
            y2 = pmr.expHandle(y2, 2/3 * pmr.logHandle( y2, M ));
            T3 = pmr.expHandle(y1, 2 * pmr.logHandle( y1, y2 ));
            T1 = pmr.expHandle(y3, 2 * pmr.logHandle( y3, y2 ));
            y1 = pmr.expHandle(y1, 0.5 * pmr.logHandle( y1, T1 ));
            y3 = pmr.expHandle(y3, 0.5 * pmr.logHandle( y3, T3 ));
            % adjust beta according to damping parameters
            beta = beta*lambda;
            
            % main iteration (gradient descent)
            for i = 0:pmr.steps-1
                
                %%%% common stuff %%%%
                % compute mean of y1 and y3, i.e., [y1, y3]_(d/2)
                tau = 1.0/(max(beta,1.0)*power(i+1,0.95 + 0.5*power(i+1,-0.18)));
%                 tau = min(tau,10/beta);
                cache = pmr.logHandle( y1, y3 );
                cache = cache*0.5;
                ym = pmr.expHandle( y1, cache );
                
                % compute direction at xm
                log_ym_y2 = pmr.logHandle( ym, y2 );
                %%%%%%%%%%%%%%%%%%%%%%
                
                %%%% gradient for y1 %%%%
                % compute direction of geodesic
                log_y3_y1 = pmr.logHandle( y3, y1 );
                
                % compute basis of eigenvectors at x3
                [basis, eigvals] = pmr.frameHandle( log_y3_y1, y3);
                
                % move basis to xm via parallel transport
                basis = pmr.parHandle( basis, y3, ym );
                
                % compute coefficients at ym
                coeffs = zeros(1,size(basis,2));
                for j = 1:size(basis,2)
                    coeffs(j) = pmr.prodHandle( basis(:,j), log_ym_y2, ym );
                end
                
                % move basis towards y1
                basis = pmr.parHandle( basis, ym, y1 );
                
                % adjust coefficients at y1
                geoLength = sqrt( pmr.prodHandle( log_y3_y1, log_y3_y1, y3 ) );
                coeffs = pmr.adjustHandle( coeffs, eigvals, geoLength );
                
                % compose gradient at y1
                baco = basis*coeffs';
                if pmr.r == 2
                    % nothing to do
                    L = beta;
                end
                if pmr.r == 1
                    L = sqrt(pmr.prodHandle( baco, baco, y1 ));
                    if L > 0.000001
                        L=beta/L;
                    else
                        L = 0;
                    end
                end
%                 grad_1 = tau*(L*baco + pmr.logHandle( y1, x1 ));
                grad_1 = tau*(beta*baco + pmr.logHandle( y1, x1 ));

                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%% gradient for y2 %%%%
                log_y2_ym = pmr.logHandle( y2, ym );
                if pmr.r == 2
                    % nothing to do
                end
                if pmr.r == 1
                    L = sqrt(pmr.prodHandle( log_y2_ym, log_y2_ym, y2 ));
                    if L > 0.000001
%                         log_y2_ym=log_y2_ym/L;
                        log_y2_ym=log_y2_ym;
                    else
                        log_y2_ym = 0;
                    end
                end
                grad_2 = tau*(beta*log_y2_ym + pmr.logHandle( y2, x2 ));
                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%% gradient for y3 %%%%
                % compute direction of geodesic
                log_y1_y3 = pmr.logHandle( y1, y3 );
                
                % compute basis of eigenvectors at y3
                [basis, eigvals] = pmr.frameHandle( log_y1_y3, y1 );
                
                % move basis to ym via parallel transport
                basis = pmr.parHandle( basis, y1, ym );
                
                % compute coefficients at ym
                coeffs = zeros(1,size(basis,2));
                for j = 1:size(basis,2)
                    coeffs(j) = pmr.prodHandle( basis(:,j), log_ym_y2, ym );
                end
                
                % move basis towards y3
                basis = pmr.parHandle( basis, ym, y3 );
                
                % adjust coefficients at y3
                geoLength = sqrt(pmr.prodHandle( log_y1_y3, log_y1_y3, y1 ));
                coeffs = pmr.adjustHandle( coeffs, eigvals, geoLength );
                
                % compose gradient at x3
                baco = basis*coeffs';
                if pmr.r == 2
                    % nothing to do
                end
                if pmr.r == 1
                    L = sqrt(pmr.prodHandle( baco, baco, y3 ) );
                    if L > 0.000001
%                         baco=baco/L;
                        baco=baco;
                    else
                        baco = 0;
                    end
                end
                grad_3 = tau*(beta*baco + pmr.logHandle( y3, x3 ));
                %%%%%%%%%%%%%%%%%
                
                %%%%% update %%%%
                %tau = 1.0/power(i+1,0.95 + 0.5*power(i+1,-0.18));
%                 grad_1 = grad_1/pmr.prodHandle(grad_1, grad_1, y1);
%                 grad_2 = grad_2/pmr.prodHandle(grad_2, grad_2, y2);
%                 grad_3 = grad_3/pmr.prodHandle(grad_3, grad_3, y3);
                y1 = pmr.expHandle( y1, grad_1 );
                y2 = pmr.expHandle( y2, grad_2 );
                y3 = pmr.expHandle( y3, grad_3 );
%                   y1 = grad_1;
%                   y2 = grad_2;
%                   y3 = grad_3;
                %%%%%%%%%%%%%%%%%
                
            end % main iteration
            
        end % function
        
    end % methods
    
end % classdefs