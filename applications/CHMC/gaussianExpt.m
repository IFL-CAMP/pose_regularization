%%
close all
clear classes
clear all

%% Setup a constrained gaussian distribution
% R = quatrot(quataxisangle([1;1;1]./sqrt(3),pi/4));
% R = eye(3);
R = [ 0.8047   -0.3106    0.5059;
      0.5059    0.8047   -0.3106;
     -0.3106    0.5059    0.8047 ];
mu = [1;1;2];
Sigma = R*diag(100./[1, 100, 10000])*R.';

A = [ 1 1 0 ];
b = 0;

% HMC setup
nllFunc = @(x) gaussianNegLikelihood(x,mu,Sigma,false);
conFunc = @(x) linearConstraint(x,A,b);
% q0 = [0;0;0];
% q0 = mu;

% Ground truth comparison
conR = [orth(A.').';null(A).'];
[mup,Sigmap] = guassianConditional(conR*mu,conR*Sigma*conR.',((length(b)+1):length(mu)).',(1:length(b)).',(A*conR(1:length(b),:).')\b);
P = conR((length(b)+1):length(mu),:);
q0 = P.'*mup;

nllpFunc = @(x) gaussianNegLikelihood(x,mup,Sigmap,false);
sampFunc = @(sz) gaussianRand(mup,Sigmap,sz);

Nscale = [1000 10000];
nRuns = 1+10;

%% Setup a constrained gaussian distribution
dim = 4;
R = eye(dim);
mu = zeros([dim,1]);
Sigma = R*diag(1./[1, 1, 100, 100])*R.';

A = [ 1 1 1 1; 1 1 -1 -1 ];
b = [ 0; 0 ];

% HMC setup
nllFunc = @(x) gaussianNegLikelihood(x,mu,Sigma,false);
conFunc = @(x) linearConstraint(x,A,b);

% Ground truth comparison
conR = [orth(A.').';null(A).'];
[mup,Sigmap] = guassianConditional(conR*mu,conR*Sigma*conR.',((length(b)+1):length(mu)).',(1:length(b)).',(A*conR(1:length(b),:).')\b);
P = conR((length(b)+1):length(mu),:);
q0 = {};
q0{1} = P.'*[2;-20];
q0{2} = P.'*mup;
% q0 = [0;0;0;0];
% q0 = mu;

nllpFunc = @(x) gaussianNegLikelihood(x,mup,Sigmap,false);
sampFunc = @(sz) gaussianRand(mup,Sigmap,sz);

%%
% Nscale = [1000 10000];
% nRuns = 1+10;
Nscale = 5000;
runTime = 15; % Time per run in seconds
nTrials = 10;

doPrint = false;

samplers = struct([]);

h = .15;
M = 1;
L = 10;
hmcOpts = struct('intMethod',0,'doPrint',doPrint,'printMod',10);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,hmcOpts);
samplers(end).N = ceil(Nscale./L);
samplers(end).style = 'b-';

h = .15;
M = 1;
L = 7;
hmcOpts = struct('intMethod',0,'doPrint',doPrint,'printMod',10);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,hmcOpts);
samplers(end).N = ceil(Nscale./L);
samplers(end).style = 'b-.';

h = .15;
M = 1;
L = 5;
hmcOpts = struct('intMethod',0,'doPrint',doPrint,'printMod',10);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,hmcOpts);
samplers(end).N = ceil(Nscale./L);
samplers(end).style = 'b:';

h = .15;
M = 1;
L = 3;
hmcOpts = struct('intMethod',0,'doPrint',doPrint,'printMod',10);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,hmcOpts);
samplers(end).N = ceil(Nscale./L);
samplers(end).style = 'b--';

h = .15;
M = 1;
L = 1;
hmcOpts = struct('intMethod',0,'doPrint',doPrint,'printMod',10);
samplers(end+1).name = sprintf('CLangevin');
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,hmcOpts);
samplers(end).N = ceil(Nscale./L);
samplers(end).style = 'g';

h = .1;
M = 1;
L = 1;
hmcOpts = struct('intMethod',0,'doPrint',doPrint,'printMod',10,'nllGradFunc',@(x) zeros(size(q0{1})));
samplers(end+1).name = sprintf('CMetropolis');
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,hmcOpts);
samplers(end).N = 2*ceil(Nscale./L);
samplers(end).style = 'r';

samplerQs = cell(length(samplers),1+nTrials);
samplerStats = cell(length(samplers),1+nTrials);
samplerNLLs = cell(length(samplers),1+nTrials);
samplerTs = zeros(length(samplers),1+nTrials);
for i = 1:numel(samplers)
    fprintf('%s...',samplers(i).name);
    [q1,qs,samplerStats{i,1}] = samplers(i).func(q0{1},samplers(i).N);
    samplerQs{i,1} = [q0{1} qs q1];
    samplerTs(i,1) = samplerStats{i,1}.t/samplerStats{i,1}.N;
    samplerNLLs{i,1} = nllFunc(samplerQs{i,1});
    
    currN = ceil(runTime./samplerTs(i,1));
    fprintf('%d samples: ',currN);
    for j = 1:nTrials
        fprintf('%d...',j);
        [q1,qs,samplerStats{i,1+j}] = samplers(i).func(q0{2},currN);
        samplerQs{i,1+j} = [q0{2} qs q1];
        samplerTs(i,1+j) = samplerStats{i,1+j}.t/samplerStats{i,1+j}.N;
        samplerNLLs{i,1+j} = nllFunc(samplerQs{i,1+j});
    end
    fprintf('\n');
end

%%
autocorrDispTime = 0.05;
tracePlotTimeStep = 0.025;
tracePlotNSteps = 20;
essautocorrTime = 0.1;

baseFigNum = 5000;
printFigs = false;

levelContours = 0.5*[1 2 3 4 5].^2;

s = 10*sqrt(max(eig(Sigmap)));
[X,Y] = meshgrid(mup(1) + linspace(-s,s,120),mup(2) + linspace(-s,s,120));
Z = reshape(nllpFunc([X(:) Y(:)].'),size(X));

samplerQsESSStats = zeros(length(q0{1}),numel(samplers),nTrials);
samplerQsESSPerSecStats = zeros(length(q0{1}),numel(samplers),nTrials);
samplerNLLESSStats = zeros(1,numel(samplers),nTrials);
samplerNLLESSPerSecStats = zeros(1,numel(samplers),nTrials);

nllacPlotParams = {};
samplerNLLAutoCorrs = cell(1,numel(samplers));
samplerNLLAutoCorrStds = cell(1,numel(samplers));
samplerNLLAutoCorrTimes = cell(1,numel(samplers));

meanCoordPlotsParams = cell([1,dim]);
for i = 1:numel(samplers)
    ishmc = isfield(samplerStats{i,1+1},'as');
    if ishmc
        allas = zeros([nTrials,samplerStats{i,1+j}.N]);
    end

    % Trace Plots
    currt = mean(samplerTs(i,:));
    cqs = samplerQs{i,1};
    qInds = [1 round((1:tracePlotNSteps)*(tracePlotTimeStep./currt))];
    
    h = figure(baseFigNum+i);
    op = get(h,'Position');
    op(3:4) = [280 210];
    plot(P(1,:)*cqs(:,qInds),P(2,:)*cqs(:,qInds),'k:o','LineWidth',1);
    axis([-10 10 -15 5]);
    set(gca,'FontSize',20);
    hold on;
    contour(X,Y,Z,levelContours);
    hold off;
    if printFigs
        set(h,'Position',op);
        set(h,'PaperPositionMode','auto');
        print(h,'-depsc','-noui',sprintf('GaussExpt-%s.eps',samplers(i).name));
    else
        title(samplers(i).name);
    end

    % AC Plots
    nlls = cat(3,samplerNLLs{i,2:end});
    samplerNLLAutoCorrTimes{i} = linspace(0,mean(samplerTs(i,:))*samplerStats{i,2}.N,samplerStats{i,2}.N+1);
    samplerNLLAutoCorrs{i} = mean(autocorr(nlls,2),3);
    samplerNLLAutoCorrStds{i} = std(autocorr(nlls,2),[],3);
    nllacPlotParams = [ nllacPlotParams, samplerNLLAutoCorrTimes{i}, samplerNLLAutoCorrs{i}, samplers(i).style ];
    
    % Mean coordinate plots
    qs =  cat(3,samplerQs{i,1});
    runavgqs = bsxfun(@rdivide,cumsum(qs,2),reshape(1:size(qs,2),[1,size(qs,2),1]));
    for j = 1:dim
        meanCoordPlotsParams{j} = [ meanCoordPlotsParams{j}, {linspace(0,mean(samplerTs(i,:))*samplerStats{i,1}.N,samplerStats{i,1}.N+1), runavgqs(j,:,:), samplers(i).style} ];
    end
end

h = figure(baseFigNum+100);
plot(nllacPlotParams{:}, 'LineWidth',2);
legend(samplers.name);
axis([0,0.25,-0.1,1]);
ylabel('Autocorrelation');
grid on;
set(gca,'FontSize',20);
drawnow
if printFigs
    op = get(h,'Position');
    op(4) = 0.4*op(3);
    set(h,'Position',op);
    set(h,'PaperPositionMode','auto');
    print(h,'-depsc','-noui','GaussExpt-NLLAutoCorrelation.eps');
else
    title('NLL Autocorrelation');
end

for j = 1:dim
    minQ = min(qs(j,:));
    maxQ = max(qs(j,:));
    
    h = figure(baseFigNum + 200 + j);
    plot(meanCoordPlotsParams{j}{:}, 'LineWidth',2);
    legend(samplers.name);
    ax = axis; ax(1:2) = [0,1];
    axis(ax);
    ylabel(sprintf('E[q_{%d}]',j));
    set(gca,'FontSize',20);
    if printFigs
        op = get(h,'Position');
        op(4) = 0.4*op(3);
        set(h,'Position',op);
        set(h,'PaperPositionMode','auto');
        print(h,'-depsc','-noui',sprintf('GaussExpt-MeanCoord%d.eps',j));
    end
end
