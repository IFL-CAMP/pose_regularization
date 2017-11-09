function [Pgt, Pnoise] = genTestData(testCase)

% define parameters for noise distributions globally
kappa = [10 50 100 500 1000];
sigma = [1 0.5 0.25 0.1 0.05];

close all

N = 99;

alpha = [zeros(1,N+1) 0:(pi)/N:pi];

beta = [zeros(1,50) -pi/2:pi/(2*49):0 0:2*pi/N:2*pi];

gamma = [zeros(1,N+1) 0:(pi)/(2*N+1):pi];

x = [5*cos(pi:pi/49:2*pi) 5*ones(1,50) 5*cos(0:2*pi/N:2*pi)];

y = [[5*sin(pi:pi/49:2*pi)-20-20/49] -20-20/49:20/49:-20/49 5*sin(0:2*pi/N:2*pi)];

z = [zeros(1,N+1) 0:10/N:10];

Pgt = zeros(4,4,2*(N+1));
Pnoise = zeros(4,4,2*(N+1));

% generate ground truth pose
for i = 1:2*(N+1)
    
    Pgt(:,:,i) = genPose(0,0,0,x(i),y(i),z(i))*genPose(alpha(i),beta(i),gamma(i),0,0,0);
    
end

% add noise to the poses, i.e. Gaussian Noise to the translational
% component and von Mises-Fisher noise to the rotational component

% extract rotational and translational components
R = Pgt(1:3,1:3,:);
t = squeeze(Pgt(1:3,4,:));

% sample from von Mises-Fisher distribution
cd CHMC
[~,qs] = binghamVonMisesFisherGibbsSampler(zeros(4),kappa(testCase)*[1,0,0,0]',2*(N+1)+1);
cd ..

% sample from Gaussian distribution and apply noise
Sigma = sigma(testCase)*eye(3);
C = chol(Sigma);
t = (t'+ randn(size(t'))*C)';

% apply noise
for i = 1:2*(N+1)
    
    % transform rotational part to quaternion
    q = rotm2quat(R(:,:,i));
    
    % apply von Mises-Fisher noise
    q = quatmultiply(q,qs(:,i)');
    
    % compose new matruix
Pnoise(:,:,i) = [[quat2rotm(q) t(:,i)]; 0 0 0 1];

end

% plotPoseMatrix(Pgt,2.0)
% figure
% plotPoseMatrix(Pnoise,2.0)
% axis equal

end % function