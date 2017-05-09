function tau = geodesicLength(lambda, dist, p, alpha)

% initialize tau
tau = 0.0;

%% reference implementation for dataPathLength
if ( nargin == 3 )
    
    % change tau in case the geodesic distance is positive
    if ( dist > 0 )
        
        %tau determines how far we shoot !
        switch (p)
            
            % L^2-norm
            case 2
                
                tau = lambda/(1+lambda);
                
            % L^1-norm => soft thresholding
            case 1
                
                if ( lambda < dist )
                    tau = lambda/dist;
                else
                    tau = 1.;
                end
                
            % Huber-norm
            case 0
                
                % old implementation of Laurent, should be checked in a future version
                tau_Huber = 1.0;
                omega = 1.0;
                
                l_t2_huber = lambda*tau_Huber*tau_Huber;
                
                if ( dist*sqrt(2.0)*tau_Huber < omega*(1+2*l_t2_huber) )
                    tau = (2*l_t2_huber)/(1+2*l_t2_huber);
                else
                    if ( sqrt(2.0)*lambda*omega*tau_Huber  < dist )
                        tau = sqrt(2.0)*lambda*omega*tau_Huber/dist;
                    else
                        tau = 1.0;
                    end
                    
                end % if(dist*sqrt(2.0)...
                
        end % switch (p)
    end % if ( dist > 0 )
end % if ( nargin == 3 )

%% reference implementation for regPathLength
if ( nargin == 4 )
    
    % change tau in case the geodesic distance is positive */
    if ( dist > 0.0 )
        
        switch (p)
            
            % L^2-norm
            case 2
                
                tau = lambda*alpha/(2*alpha*lambda+1);
                
            % L^1-norm => soft thresholding
            case 1
                tau = lambda*alpha/dist;
                if(tau > 0.5)
                    tau = 0.5;
                end
                
            % Huber-norm
            case 0
                
                % old implementation of Laurent, should be checked in a future version
                tau_Huber = 1.0;
                omega = 1.0;
                
                if ( dist < omega*(1+4*lambda*alpha*tau_Huber*tau_Huber)/(sqrt(2.0)*tau_Huber) )
                    
                    tau = (2*alpha*lambda*tau_Huber*tau_Huber)/(1+4*alpha*lambda*tau_Huber*tau_Huber);
                    
                else
                    
                    tau = sqrt(2.0)*alpha*lambda*omega*tau_Huber/dist;
                    if(tau > 0.5)
                        tau = 0.5;
                    end
                end % if ( dist < omega*...
                
        end % end switch (p)
    end % end if ( dist > 0.0 )
    
end % if ( nargin == 4 )

