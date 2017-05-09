%%
close all
clear classes
clear all


%% Setup the target posterior
dim = 4;

V = eye(dim);

% scale the eigenvalues of A and the magnitude of c
Ascale = 1000; cscale = 10;
A = Ascale*V*diag(linspace(-1,1,dim))*V.';
c = zeros([dim,1]); c(1) = cscale;
c = V.'*c;

% randomly pick a starting point by uniformly sampling from the space of
% unit length vectors.
q0 = binghamVonMisesFisherGibbsSampler(A,c,100);

nllFunc = @(x) binghamVonMisesFisherNegLikelihood(x,A,c);
conFunc = @(x) quatConst(x);


%% Setup and run multiple samplers
s = eig(A);

scaleN = 10000;

hmcN = 2*scaleN;
langevinN = 4*scaleN;
metropolisN = 6*scaleN;
gibbsN = scaleN;

samplers = struct([]);

h = 1;
L = 4;
M = max(s) - min(s);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,struct('doPrint',2));
samplers(end).N = hmcN;
samplers(end).style = 'b-';

L = 3;
M = max(s) - min(s);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,struct('doPrint',2));
samplers(end).N = hmcN;
samplers(end).style = 'b--';

L = 2;
M = max(s) - min(s);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,struct('doPrint',2));
samplers(end).N = hmcN;
samplers(end).style = 'b-.';

h = 1;
M = (max(s) - min(s));
samplers(end+1).name = 'CLangevin';
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,1,h,struct('doPrint',2));
samplers(end).N = langevinN;
samplers(end).style = 'g';

h = 0.4;
M = (max(s) - min(s));
samplers(end+1).name = 'CMetropolis';
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,1,h,struct('nllGradFunc',@(q) zeros(size(q)),'doPrint',2));
samplers(end).N = metropolisN;
samplers(end).style = 'r-';

% samplers(end+1).name = 'Gibbs (Rejection Sampler)';
% samplers(end).func = @(q,N) binghamVonMisesFisherGibbsSampler(A,c,N,q,true);
% samplers(end).N = gibbsN;

samplers(end+1).name = 'Gibbs';
samplers(end).desc = [];
samplers(end).func = @(q,N) binghamVonMisesFisherGibbsSampler(A,c,N,q,false);
samplers(end).N = gibbsN;
samplers(end).style = 'c';

nTrials = 1;
samplerQs = cell(length(samplers),nTrials);
samplerStats = cell(length(samplers),nTrials);
samplerNLLs = cell(length(samplers),nTrials);
samplerTs = zeros(length(samplers),nTrials);
for i = 1:numel(samplers)
% for i = 6
    for j = 1:nTrials
        fprintf('\n%s\n',samplers(i).name);
        [q1,qs,samplerStats{i,j}] = samplers(i).func(q0,samplers(i).N);
        samplerQs{i,j} = [qs q1];
        samplerTs(i,j) = samplerStats{i,j}.t/samplerStats{i,j}.N;
        samplerNLLs{i,j} = nllFunc(samplerQs{i,j});
    end
end

%% Plot autocorrelation of the NLL and print stats based on the NLL
autocorrDispTime = 0.05;
essautocorrTime = 0.1;
printFigs = false;

samplerESSStats = zeros(dim,numel(samplers),nTrials);

samplerNLLAutoCorrs = cell(1,numel(samplers));
samplerNLLAutoCorrStds = cell(1,numel(samplers));
samplerNLLAutoCorrTimes = cell(1,numel(samplers));
for i = 1:numel(samplers)
    nlls = cat(2,samplerNLLs{i,:}).';

    samplerNLLAutoCorrTimes{i} = linspace(0,samplerStats{i,1}.t,samplerStats{i,1}.N);
    samplerNLLAutoCorrs{i} = mean(autocorr(nlls,2),1);
    samplerNLLAutoCorrStds{i} = std(autocorr(nlls,2),[],1);
end

plotArgs = [samplerNLLAutoCorrTimes;samplerNLLAutoCorrs;{samplers.style}];
h = figure(1);
op = get(h,'Position');
op(4) = 0.4*op(3);
set(h,'Position',op);
set(h,'PaperPositionMode','auto');
plot(plotArgs{:},'LineWidth',2);
axis([0 autocorrDispTime -0.1 1]);
legend(samplers.name);
grid on;
set(gca,'FontSize',14);
drawnow
if printFigs
    print -depsc -noui BinghamExpt-AutoCorrelation.eps
end

