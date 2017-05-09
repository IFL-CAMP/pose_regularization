function x = proximalPointAlgorithm( f, gridDimX, gridDimY, gridDimZ,...
                                     p, q, r, alpha, beta, gamma, steps,...
                                     logHandle, expHandle, prodHandle,...
                                     frameHandle, parHandle, adjustHandle )
                                            
% define number of inner steps for second order regularization
innerSteps = 1;
                                            
% initialize proximal mappings
proxData = proximalMapData( logHandle, expHandle, prodHandle, p );
proxReg1st = proximalMapRegFirstOrder( logHandle, expHandle, prodHandle, q );
proxReg2nd = proximalMapRegSecondOrder( logHandle, expHandle, prodHandle,...
                                        frameHandle, parHandle, adjustHandle, ...
                                        r, innerSteps );
% initialize x
x = f;

% get element length
elementLength = size(f,1);

% main iteration
for ll = 0:1:steps-1
    
    % compute lambda
    lambda = 1.0/power(ll+1,0.95 + 0.5*power(ll+1,-0.18));
    
    % proximal mapping of data term
    for k = 0:1:gridDimZ-1
        for j = 0:1:gridDimY-1
            for i = 0:1:gridDimX-1
                index = (k*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                cache = proxData.eval( x(index:index+elementLength-1),...
                    f(index:index+elementLength-1),...
                    lambda );
                x(index:index+elementLength-1) = cache;
            end
        end
    end
    
    % proximal mapping for 1st order regularizer
    if alpha > 0
        
        % proximal mapping of 1st order regularizer w.r.t. to x-direction
        for k = 0:1:gridDimZ-1
            for j = 0:1:gridDimY-2
                for i = 0:1:gridDimX-1
                    index1 = (k*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    index2 = (k*gridDimX*gridDimY + (j+1)*gridDimX + i)*elementLength + 1;
                    [cache1, cache2] = proxReg1st.eval( x(index1:index1+elementLength-1),...
                                                     x(index2:index2+elementLength-1),...
                                                     alpha, lambda );
                    x(index1:index1+elementLength-1) = cache1;
                    x(index2:index2+elementLength-1) = cache2;
                end
            end
        end % proximal mapping of 1st order regularizer w.r.t. to x-direction

        % proximal mapping of 1st order regularizer w.r.t. to y-direction
        for k = 0:1:gridDimZ-1
            for j = 0:1:gridDimY-1
                for i = 0:1:gridDimX-2
                    index1 = (k*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    index2 = (k*gridDimX*gridDimY + j*gridDimX + (i+1))*elementLength + 1;
                    [cache1, cache2] = proxReg1st.eval( x(index1:index1+elementLength-1),...
                                                     x(index2:index2+elementLength-1),...
                                                     alpha, lambda );
                    x(index1:index1+elementLength-1) = cache1;
                    x(index2:index2+elementLength-1) = cache2;
                end
            end
        end % proximal mapping of 1st order regularizer w.r.t. to y-direction

        % proximal mapping of 1st order regularizer w.r.t. to z-direction
        for k = 0:1:gridDimZ-2
            for j = 0:1:gridDimY-1
                for i = 0:1:gridDimX-1
                    index1 = (k*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    index2 = ((k+1)*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    [cache1, cache2] = proxReg1st.eval( x(index1:index1+elementLength-1),...
                                                     x(index2:index2+elementLength-1),...
                                                     alpha, lambda );
                    x(index1:index1+elementLength-1) = cache1;
                    x(index2:index2+elementLength-1) = cache2;
                end
            end
        end % proximal mapping of 1st order regularizer w.r.t. to z-direction
    end % if alpha > 0
    
    % proximal mapping of 2nd order regularizer
    if ( beta > 0)
        
      
        % proximal mapping of 2nd order regularizer w.r.t. to x-direction
        for k = 0:1:gridDimZ-1
            for j = 0:1:gridDimY-1
                for i = 1:1:gridDimX-2
                    index1 = (k*gridDimX*gridDimY + j*gridDimX + (i-1))*elementLength + 1;
                    index2 = (k*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    index3 = (k*gridDimX*gridDimY + j*gridDimX + (i+1))*elementLength + 1;
                    [cache1, cache2, cache3] = proxReg2nd.eval( x(index1:index1+elementLength-1),...
                                                                x(index2:index2+elementLength-1),...
                                                                x(index3:index3+elementLength-1),...
                                                                beta, lambda );
                    x(index1:index1+elementLength-1) = cache1;
                    x(index2:index2+elementLength-1) = cache2;
                    x(index3:index3+elementLength-1) = cache3;
                end
            end
        end % proximal mapping of 2nd order regularizer w.r.t. to x-direction
        
        % proximal mapping of 2nd order regularizer w.r.t. to y-direction
        for k = 0:1:gridDimZ-1
            for j = 1:1:gridDimY-2
                for i = 0:1:gridDimX-1
                    index1 = (k*gridDimX*gridDimY + (j-1)*gridDimX + i)*elementLength + 1;
                    index2 = (k*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    index3 = (k*gridDimX*gridDimY + (j+1)*gridDimX + i)*elementLength + 1;
                    [cache1, cache2, cache3] = proxReg2nd.eval( x(index1:index1+elementLength-1),...
                                                                x(index2:index2+elementLength-1),...
                                                                x(index3:index3+elementLength-1),...
                                                                beta, lambda );
                    x(index1:index1+elementLength-1) = cache1;
                    x(index2:index2+elementLength-1) = cache2;
                    x(index3:index3+elementLength-1) = cache3;
                end
            end
        end % proximal mapping of 2nd order regularizer w.r.t. to y-direction

        % proximal mapping of 2nd order regularizer w.r.t. to z-direction
        for k = 1:1:gridDimZ-2
            for j = 0:1:gridDimY-1
                for i = 0:1:gridDimX-1
                    index1 = ((k-1)*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    index2 = (k*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    index3 = ((k+1)*gridDimX*gridDimY + j*gridDimX + i)*elementLength + 1;
                    [cache1, cache2, cache3] = proxReg2nd.eval( x(index1:index1+elementLength-1),...
                                                                x(index2:index2+elementLength-1),...
                                                                x(index3:index3+elementLength-1),...
                                                                beta, lambda );
                    x(index1:index1+elementLength-1) = cache1;
                    x(index2:index2+elementLength-1) = cache2;
                    x(index3:index3+elementLength-1) = cache3;
                end
            end
        end % proximal mapping of 2nd order regularizer w.r.t. to z-direction
    end % if(beta > 0)
    
    % REMARK: mixed proximal mappings (for diagonal
    % directions) are not implemented for pose dnoising!
    
end % main iteration

end