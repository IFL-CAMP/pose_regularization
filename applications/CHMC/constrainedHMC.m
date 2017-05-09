function [q0, qs, stats] = constrainedHMC(q0,nllFunc,conFunc,Mfunc,N,L,h,opts)
    baseOpts = struct('printMod',1,'nllGradFunc',[],'doPrint',1,'distanceFunc',@(x0,x1) norm(x0-x1));
    if nargin < 8 || isempty(opts)
        opts = baseOpts;
    else
        opts = getOpts(baseOpts,opts);
    end
    printMod = opts.printMod;
    nllGradFunc = opts.nllGradFunc;
    doPrint = opts.doPrint;
    distanceFunc = opts.distanceFunc;

    q0 = reshape(q0,[],1);
    if nargout > 1
        qs = zeros(length(q0),N-1);
    end
    if nargout > 2
        stats.Ls = zeros(1,N);
        stats.hs = zeros(1,N);
        stats.dHs = zeros(1,N);
        stats.Us = zeros(1,N);
        stats.as = zeros(1,N);
        stats.dxs = zeros(1,N);
        stats.accepted = false(1,N);
    end
    stats.N = N;
    stats.conFuncCalls = 0;
    stats.nllFuncCalls = 0;
    stats.nllGradFuncCalls = 0;
    stats.MFuncCalls = 0;
    stats.MGradFuncCalls = 0;

    startT = tic;

    
    if ~isempty(conFunc)
        [c0,dc0] = conFunc(q0);
        stats.conFuncCalls = stats.conFuncCalls + 1;
    else
        c0 = [];
        dc0 = [];
    end
    if isempty(nllGradFunc)
        [U0,dU0] = nllFunc(q0);
        stats.nllFuncCalls = stats.nllFuncCalls + 1;
        stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
    else
        U0 = nllFunc(q0);
        dU0 = nllGradFunc(q0);
        stats.nllFuncCalls = stats.nllFuncCalls + 1;
        stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
        assert(length(dU0) == length(q0))
    end

    if isnumeric(Mfunc)
        M0 = Mfunc;
    else
        [M0,dM0] = Mfunc(q0);
        stats.MFuncCalls = stats.MFuncCalls + 1;
        stats.MGradFuncCalls = stats.MGradFuncCalls + 1;
    end
    [cholM0,Minv0] = factorizeMassMatrix(M0);
    logNormConst0 = computeLogNormConst(cholM0,Minv0,dc0);

    if isnumeric(h)
        ch = h;
    end
    
    pStepSum = 0;
    for j = 1:N
        prevq = q0;
        prevU = U0;
        if nargout > 1 && j > 1
            qs(:,j-1) = q0;
        end
        
        p0 = cholM0*randn(size(q0));
        if ~isempty(conFunc)
            D = dc0*Minv0;
            p0 = p0 - D.'*((D*D.')\(D*p0));
        end
        
        H0 = computeH(p0,Minv0,U0,logNormConst0);

        if doPrint == 1 && (mod(j,printMod)==0)
            fprintf('%d: ', j);
        end
        
        if ~isnumeric(L)
            cL = L();
            if doPrint == 1 && (mod(j,printMod)==0)
                fprintf('L = %d, ', cL);
            end
        else
            cL = L;
        end
        if ~isnumeric(h)
            ch = h();
            if doPrint == 1 && (mod(j,printMod)==0)
                fprintf('h = %g, ', abs(ch));
            end
        end
        if isnumeric(Mfunc)
            [q,p,c,dc,U,dU,stats] = constantMInt_Leapfrog(nllFunc,conFunc,cL,ch,q0,p0,Minv0,c0,dc0,U0,dU0,nllGradFunc,stats);
            cholM = cholM0;
            Minv = Minv0;
        else
            [q,p,M,dM,c,dc,U,dU,stats] = variableMInt_Leapfrog(nllFunc,conFunc,Mfunc,cL,ch,q0,p0,M0,dM0,c0,dc0,U0,dU0,nllGradFunc,stats);
            [cholM,Minv] = factorizeMassMatrix(M);
        end
        
        if norm(c) > sqrt(1e-3) || isempty(p)
            H = inf;
            deltaU = inf;
            deltaH = inf;
        else
            logNormConst = computeLogNormConst(cholM,Minv,dc);
            H = computeH(p,Minv,U,logNormConst);
            %             H = computeH(p,cholM,Minv,U,dc);
            deltaU = U-U0;
            deltaH = H-H0;
        end
        
        currA = min([1,exp(-deltaH)]);
        pStepSum = pStepSum + currA;
        if (H < inf) && (~isnan(H)) && (deltaH < 0 || erand > deltaH)
            accept = 1;
        else
            accept = 0;
        end
        
        if accept
            q0 = q;
            c0 = c;
            dc0 = dc;
            U0 = U;
            dU0 = dU;
            logNormConst0 = logNormConst;
            if ~isnumeric(Mfunc)
                % only update mass matrix if it's state dependant
                cholM0 = cholM;
                M0 = M;
                Minv0 = Minv;
                dM0 = dM;
            end
        end
        
        cdx = distanceFunc(prevq,q0);
        if (doPrint == 1) && (mod(j,printMod)==0)
            if ~accept
                if isempty(q)
                    accStr = 'convergence failure, ';
                else
                    accStr = sprintf('dPropU = %g, ||x-propx|| = %g, ',deltaU,distanceFunc(q0,q));
                end
            else
                accStr = '';
            end
            if ~isempty(conFunc)
                projdU = dU0 - dc0.'*((dc0*dc0.')\(dc0*dU0));
            else
                projdU = dU0;
            end
            fprintf('U = %g, dU = %g, dH = %g, ||dudx|| = %g, ||c|| = %g, ||dx|| = %g, %sa = %g, E[a] = %g\n',U0,U0-prevU,deltaH,norm(projdU),norm(c0),cdx,accStr,currA,pStepSum/j);
            drawnow;
        end
        
        if nargout > 2
            stats.Ls(j) = cL;
            stats.hs(j) = ch;
            stats.dHs(j) = deltaH; % Proposed dH
            stats.Us(j) = U0;
            stats.dxs(j) = cdx;
            stats.as(j) = min([1,exp(-deltaH)]);
            stats.accepted(j) = accept;
        end
    end
    if doPrint == 2
        fprintf('E[p(accept)] = %g\n',pStepSum/N);
    end
    stats.t = toc(startT);
end

function opts = getOpts(opts,newOpts)
    S = fieldnames(newOpts);
    for i = 1:length(S)
        opts.(S{i}) = newOpts.(S{i});
    end
end

function [cholM,Minv] = factorizeMassMatrix(M,isinvM)
    if nargin < 2
        isinvM = false;
    end
    if isinvM
        if isscalar(M)
            cholM = 1./sqrt(M);
            Minv = 1./M;
        elseif isvector(M)
            cholM = spdiags(reshape(1./sqrt(M),[],1),0,length(M),length(M));
            Minv = spdiags(reshape(1./M,[],1),0,length(M),length(M));
        else
            cholM = chol(M,'lower');
            Minv = (cholM.'\(cholM\speye(size(M))));
            cholM = chol(M,'lower');
        end
    else
        if isscalar(M)
            cholM = sqrt(M);
            Minv = 1./M;
        elseif isvector(M)
            cholM = spdiags(reshape(sqrt(M),[],1),0,length(M),length(M));
            Minv = spdiags(reshape(1./M,[],1),0,length(M),length(M));
        else
            cholM = chol(M,'lower');
            Minv = (cholM.'\(cholM\speye(size(M))));
        end
    end
end

function logNormConst = computeLogNormConst(cholM,Minv,dc)
    if isscalar(cholM)
        logNormConst = (size(dc,2) - size(dc,1))*log(cholM);
    else
        D = dc*Minv;
        if isvector(cholM)
            cholMhat = spdiag(cholM) - D.'*((D*D.')\(D*spdiag(cholM)));
        else
            cholMhat = cholM - D.'*((D*D.')\(D*cholM));
        end
        % FIXME: there is almost certainly a better way to do this...
        s = eigs(cholMhat*cholMhat.',size(dc,2)-size(dc,1),'lm');
        logNormConst = 0.5*sum(log(s));
    end
end

function H = computeH(p,Minv,U,logNormConst)
    H = 0.5*p.'*Minv*p + logNormConst + U;
end

function [q,p,c,dc,U,dU,stats] = constantMInt_Leapfrog(nllFunc,conFunc,L,h,q,p,Minv,c,dc,U,dU,nllGradFunc,stats)
    fThresh = 1e-5;
    dxThresh = 1e-6;
    maxIts = 100;

    if ~isempty(conFunc) && (isempty(c) || isempty(dc))
        [c,dc] = conFunc(q);
        stats.conFuncCalls = stats.conFuncCalls + 1;
    end
    if isempty(U) || isempty(dU)
        if isempty(nllGradFunc)
            [U,dU] = nllFunc(q);
        else
            U = nllFunc(q);
            dU = nllGradFunc(q);
        end
        stats.nllFuncCalls = stats.nllFuncCalls + 1;
        stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
    end


    for i = 1:L
        if ~isempty(conFunc)
            A = -(h/2)*(Minv*dc.');
            %A = -((h*h)/2)*(Minv*dc.');
            p12hat = p - (h/2)*dU;
            q12hat = q + h*(Minv*p12hat);
            
            %%%%%% START simple_fsolve %%%%%%
%             curr_q = q12hat + A*lambda;
            lambda = zeros([length(c),1]);
            curr_q = q12hat;
            [curr_c,curr_dc] = conFunc(curr_q);
            stats.conFuncCalls = stats.conFuncCalls + 1;
            it = 0;
            while max(abs(curr_c)) > fThresh
                it = it + 1;

                dlambda = (curr_dc*A)\curr_c;
                lambda = lambda - dlambda;

                curr_q = q12hat + A*lambda;
                [curr_c,curr_dc] = conFunc(curr_q);
                stats.conFuncCalls = stats.conFuncCalls + 1;
%                 curr_q
%                 curr_c

                if max(abs(dlambda)) < dxThresh && max(abs(curr_c)) > fThresh
%                     warning('Change in x is below threshold, bailing out');
%                     break;
                    q = [];
                    p = [];
                    c = [];
                    dc = [];
                    U = [];
                    dU = [];
                    return;
                end

                if it >= maxIts && max(abs(curr_c)) > fThresh
%                     warning('Maximum number of iterations hit, bailing out');
%                     break;
                    q = [];
                    p = [];
                    c = [];
                    dc = [];
                    U = [];
                    dU = [];
                    return;
                end
            end
            %%%%%% END simple_fsolve %%%%%%

            p12 = p12hat - dc.'*((1/2)*lambda);
            %p12 = p12hat - dc.'*((h/2)*lambda);
            
            if max(abs(curr_c)) > fThresh
                break;
            end
            
            q = curr_q;
            c = curr_c;
            dc = curr_dc;
        else
            p12 = p - (h/2)*dU;
            q = q + h*(Minv*p12);
            if(any(isnan(p12)))
                q = [];
                p = [];
                U = [];
                dU = [];
                return;
            end
        end

        if isempty(nllGradFunc)
            [U,dU] = nllFunc(q);
            stats.nllFuncCalls = stats.nllFuncCalls + 1;
            stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
        else
            dU = nllGradFunc(q);
            stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
        end

        if ~isempty(conFunc)
            A = (h/2)*dc*(Minv*dc.');
            b = dc*(Minv*(p12 - (h/2)*dU));
            mu = A\b;
            p = p12 - (h/2)*(dU + dc.'*mu);
        else
            p = p12 - (h/2)*dU;
        end
    end
    if ~isempty(nllGradFunc)
        U = nllFunc(q);
        stats.nllFuncCalls = stats.nllFuncCalls + 1;
    end
end

function [q,p,M,dM,c,dc,U,dU,stats] = variableMInt_Leapfrog(nllFunc,conFunc,Mfunc,L,h,q,p,M,dM,c,dc,U,dU,nllGradFunc,stats)
    fThresh = 1e-4;
    dxThresh = 1e-6;
    maxIts = 50;

    if ~isempty(conFunc)&&(isempty(c) || isempty(dc))
        [c,dc] = conFunc(q);
        stats.conFuncCalls = stats.conFuncCalls + 1;
    end
    if isempty(U) || isempty(dU)
        if isempty(nllGradFunc)
            [U,dU] = nllFunc(q);
        else
            U = nllFunc(q);
            dU = nllGradFunc(q);
        end
        stats.nllFuncCalls = stats.nllFuncCalls + 1;
        stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
    end
    if isempty(M) || isempty(dM)
        [M,dM] = Mfunc(q);
        stats.MFuncCalls = stats.MFuncCalls + 1;
        stats.MGradFuncCalls = stats.MGradFuncCalls + 1;
    end

    D = size(dc,2);

    Minv = inv(M);
    MinvdM = vec_matmult(Minv,0,dM,0);
    pMat = reshape(vec_matmult(MinvdM,0,Minv,0),D^2,D).';
    MinvdM = reshape(MinvdM,D^2,D);
    trMinvdM = sum(MinvdM(1:(D+1):D^2,:),1).';

    for i = 1:L
        p12 = p - (h/2)*(-0.5*pMat*reshape(p*p.',[],1) + 0.5*trMinvdM + dU);
        q1 = q + h*(Minv*p12);

        lambda = zeros([length(c),1]);
        M1inv = Minv; %inv(Mfunc(q1));
        [c1,dc1] = conFunc(q1);
        stats.MFuncCalls = stats.MFuncCalls + 1;
        stats.conFuncCalls = stats.conFuncCalls + 1;
        
        nIts = 0;
        while true
            %J = dc1*((h/2)*(Minv + M1inv)*(-(h/2)*dc.'));
            J = dc1*((h/4)*(Minv + M1inv)*(-dc.'));
            dlambda = J\(-c1);
            lambda = lambda + dlambda;

            M1inv = inv(Mfunc(q1));

            dp12 = (p - (h/2)*(-0.5*pMat*reshape(p12*p12.',[],1) + 0.5*trMinvdM + dU) - (1/2)*dc.'*lambda) - p12;
            p12 = p12 + dp12;
            
            dq1 = (q + (h/2)*(Minv + M1inv)*p12) - q1;
            q1 = q1 + dq1;
            
            [c1,dc1] = conFunc(q1);
            stats.MFuncCalls = stats.MFuncCalls + 1;
            stats.conFuncCalls = stats.conFuncCalls + 1;

            if max(abs(c1)) < fThresh && norm(dq1) < dxThresh && norm(dp12) < dxThresh
                break;
            end
            
            nIts = nIts + 1;
            if nIts >= maxIts
                break;
            end
        end
        if nIts >= maxIts
            q = [];
            p = [];
            c = [];
            dc = [];
            U = [];
            dU = [];
            M = [];
            dM = [];
            break;
        end
        
        q = q1;
        c = c1;
        dc = dc1;

        [M,dM] = Mfunc(q);
        stats.MFuncCalls = stats.MFuncCalls + 1;
        stats.MGradFuncCalls = stats.MGradFuncCalls + 1;

        Minv = inv(M);
        M1inv = Minv;
        if isempty(nllGradFunc)
            [U,dU] = nllFunc(q);
            stats.nllFuncCalls = stats.nllFuncCalls + 1;
            stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
        else
            dU = nllGradFunc(q);
            stats.nllGradFuncCalls = stats.nllGradFuncCalls + 1;
        end

        MinvdM = vec_matmult(Minv,0,dM,0);
        pMat = reshape(vec_matmult(MinvdM,0,Minv,0),D^2,D).';
        MinvdM = reshape(MinvdM,D^2,D);
        trMinvdM = sum(MinvdM(1:(D+1):D^2,:),1).';
        dH = -0.5*pMat*reshape(p12*p12.',[],1) + 0.5*trMinvdM + dU;
    
        A = (h/2)*dc*M1inv*dc.';
        b = dc*M1inv*(p12 - (h/2)*dH);
        mu = A\b;
        p = p12 - (h/2)*(dH + dc.'*mu);
    end
    
    if ~isempty(nllGradFunc)
        U = nllFunc(q);
        stats.nllFuncCalls = stats.nllFuncCalls + 1;
    end
end


function x = erand(sz)
    if nargin > 0
        x = -log(rand(sz));
    else
        x = -log(rand);
    end
end
