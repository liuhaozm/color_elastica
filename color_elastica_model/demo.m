clear

f=double(imread('chipCrop.jpg'))/255;

[m,n,k]=size(f);
f0=f;

%%% Gaussian noise
noiseLevel=0.06;
% f=f0+randn(m,n,k)*sqrt(noiseLevel);
f=f0+randn(m,n,k)*(noiseLevel);

% %%%%%%%% Poisson noise
% fac=1*10^10;
% f=double(imnoise(f0/fac,'poisson')*fac);


f=min(f,1);
f=max(f,0);

alpha1=1e-2; 
beta=0.005;
eta=1;
tol=1e-2;
tol1=1e-5;
tau=5e-2;

u=CE(f,alpha1,beta,eta,tau,tol,tol1);
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
title('our result')