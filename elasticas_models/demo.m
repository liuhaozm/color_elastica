clear

%% read image
f=double(imread('onion.jpg'))/255;

[m,n,k]=size(f);
f0=f;

%%% Gaussian noise
noiseLevel=0.15;
% f=f0+randn(m,n,k)*sqrt(noiseLevel);
f=f0+randn(m,n,k)*(noiseLevel);

%%
% choose the model you want to test in the following and comment out the
% other model

% First modified model

alpha1=5e-4;
beta=50; % weight of color elastica
eta=20; % parameter of the fidelity term
tol=5e-5; % tolorance of the full algorithm
tol1=1e-5; % tolorance of the fixed point iteration for p
tau=5e-2; % time step
tic
u=FirstModifiedModel(f,alpha1,beta,eta,tau,tol,tol1);
toc
%% second modified model
% alpha1=5e-3;
% beta=30; % weight of color elastica
% eta=2.4; % parameter of the fidelity term
% tol=1e-5; % tolorance of the full algorithm
% tol1=1e-5; % tolorance of the fixed point iteration for p
% tau=5e-2; % time step
% 
% tic
% u=SecondModifiedModel(f,alpha1,beta,eta,tau,tol,tol1);
% toc
%%

figure
subplot(1,3,1)
imshow(f,[0,1])
title('noisy image')
subplot(1,3,2)
imshow(f0,[0,1])
title('original image')
subplot(1,3,3)
imshow(u,[0,1])
title('denoised image')