% MATLAB code for "Elastica Models for Color Image Regularization" by Liu, Hao 
% and Tai, Xue-Cheng and Kimmel, Ron and Glowinski, Roland .
% https://arxiv.org/abs/2203.09995
%
% If this code is useful, pleae cite our paper.
%

% Copyright (c) 2022 Hao Liu (haoliu AT hkbu.edu.hk)
% Department of Mathematics,
% Hong Kong Baptist University
% Kowloon, Hong Kong
% https://www.math.hkbu.edu.hk/~haoliu/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

function u=FirstModifiedModel(f,alpha1,beta,eta,tau,tol,tol1)
% The first modified model in our paper

% parameters:
% alpha1: alpha in the model 
% beta: beta in the model
% eta: eta in the model
% tau: time step.
% tol: stopping critetion for the algorithm
% tol1: stopping criterion for the fixed point iteration for p

[m,n,k]=size(f);
gamma1=1;
gamma2=3;



h=1;
u=f;
[p1,p2]=gradu(u);


[G11,G12,G22]=Mp(p1,p2,alpha1);



        g=G11.*G22-G12.^2;



[lambda1,lambda2]=initialLambda(g,G11,G12,G22,p1,p2);

xx=(1:m)-1;
yy=(1:n)-1;
[y1,x1]=meshgrid(yy,xx);
zx=2*pi*x1/m;
zy=2*pi*y1/n;

w0=(1-exp(-sqrt(-1)*zx)).*(exp(sqrt(-1)*zx)-1)+...
    (1-exp(-sqrt(-1)*zy)).*(exp(sqrt(-1)*zy)-1);
w=eta*((1-exp(-sqrt(-1)*zx)).*(exp(sqrt(-1)*zx)-1)+...
    (1-exp(-sqrt(-1)*zy)).*(exp(sqrt(-1)*zy)-1))-tau;


ss=10;
kk=0;
err=[];
energy=[];

while ss>tol
    kk=kk+1;
    %%%%%% step 1
    dLambda=divLambda(lambda1,lambda2);
    s1=sum(dLambda.^2,3);
    s=zeros(m,n);
    [p1_1,p2_1,Iterp]=updatepFirst(p1,p2,s1,alpha1,beta,tau,tol1);
    [G11temp,G12temp,G22temp]=Mp(p1_1,p2_1,alpha1);
    G11_1=exp(-gamma2*tau)*G11+(1-exp(-gamma2*tau))*G11temp;
    G12_1=exp(-gamma2*tau)*G12+(1-exp(-gamma2*tau))*G12temp;
    G22_1=exp(-gamma2*tau)*G22+(1-exp(-gamma2*tau))*G22temp;

            g_1=G11_1.*G22_1-G12_1.^2;

    
    tau_lambda_temp=2*tau*beta.*sqrt(g_1);
    tau_lambda=max(abs(tau_lambda_temp(:)));

    gg=G11temp.*G22temp-G12temp.^2;
    
    rhs=(tau_lambda-2*tau*beta.*sqrt(gg)).*dLambda;
    rhsex=expandf3(rhs);
    rhs_1=gamma1*lambda1-(rhsex(3:end,2:end-1,:)-rhsex(2:end-1,2:end-1,:));
    rhs_2=gamma1*lambda2-(rhsex(2:end-1,3:end,:)-rhsex(2:end-1,2:end-1,:));
    
    rhs_1F=zeros(m,n,k);
    rhs_2F=zeros(m,n,k);
    for ii=1:3
        rhs_1F(:,:,ii)=fft2(rhs_1(:,:,ii));
        rhs_2F(:,:,ii)=fft2(rhs_2(:,:,ii));
    end
    
    a11=gamma1-2*tau_lambda*(cos(zx)-1);
    a12=tau_lambda*(cos(zx)-1+sqrt(-1)*sin(zx)).*(cos(zy)-1-sqrt(-1)*sin(zy));
    a21=tau_lambda*(cos(zy)-1+sqrt(-1)*sin(zy)).*(cos(zx)-1-sqrt(-1)*sin(zx));
    a22=gamma1-2*tau_lambda*(cos(zy)-1);
    D=gamma1^2+2*gamma1*tau_lambda*(2-cos(zx)-cos(zy));
    
    lambda1_1=zeros(m,n,k);
    lambda2_1=zeros(m,n,k);
    for ii=1:3
        lambda1_1(:,:,ii)=real(ifft2((a22.*rhs_1F(:,:,ii)-a12.*rhs_2F(:,:,ii))./D));
        lambda2_1(:,:,ii)=real(ifft2((-a21.*rhs_1F(:,:,ii)+a11.*rhs_2F(:,:,ii))./D));
    end
    

    %%%%%%%%%%% step 2
    A11=2*G11_1.^2./g_1+2*G12_1.^2./g_1+2*gamma1;
    A12=2*G11_1.*G12_1./g_1+2*G12_1.*G22_1./g_1;
    A22=2*G12_1.^2./g_1+2*G22_1.^2./g_1+2*gamma1;
    
    b1=2*G11_1./sqrt(g_1).*p1_1+2*G12_1./sqrt(g_1).*p2_1+2*gamma1*lambda1_1;
    b2=2*G12_1./sqrt(g_1).*p1_1+2*G22_1./sqrt(g_1).*p2_1+2*gamma1*lambda2_1;
    
    detA=A11.*A22-A12.^2;
    
    lambda1_2=(A22.*b1-A12.*b2)./detA;
    lambda2_2=(-A12.*b1+A11.*b2)./detA;
    
    p1_2=(G11_1.*lambda1_2+G12_1.*lambda2_2)./sqrt(g_1);
    p2_2=(G12_1.*lambda1_2+G22_1.*lambda2_2)./sqrt(g_1);
    
    [G11temp,G12temp,G22temp]=Mp(p1_2,p2_2,alpha1);
    G11_2=exp(-gamma2*tau)*G11_1+(1-exp(-gamma2*tau))*G11temp;
    G12_2=exp(-gamma2*tau)*G12_1+(1-exp(-gamma2*tau))*G12temp;
    G22_2=exp(-gamma2*tau)*G22_1+(1-exp(-gamma2*tau))*G22temp;

            g_2=G11_2.*G22_2-G12_2.^2;

    dp=divLambda(p1_2,p2_2);

    %%%%%%%%%% step 3
    dp=divLambda(p1_2,p2_2);
    
    rhs=eta*dp-tau*f;
    
    rhs_F=zeros(m,n,k);
    for ii=1:3
        rhs_F(:,:,ii)=fft2(rhs(:,:,ii))./w;
    end
    
    u_3=zeros(m,n,k);
    for ii=1:3
        u_3(:,:,ii)=real(ifft2(rhs_F(:,:,ii)));
    end

    
    [p1_3,p2_3]=gradu(u_3);
    [G11temp,G12temp,G22temp]=Mp(p1_3,p2_3,alpha1);
    G11_3=exp(-gamma2*tau)*G11_2+(1-exp(-gamma2*tau))*G11temp;
    G12_3=exp(-gamma2*tau)*G12_2+(1-exp(-gamma2*tau))*G12temp;
    G22_3=exp(-gamma2*tau)*G22_2+(1-exp(-gamma2*tau))*G22temp;

            g_3=G11_3.*G22_3-G12_3.^2;

    %%%%%%%%%%%%%% energy
    gtemp=G11temp.*G22temp-G12temp.^2;
    mu1=(G22temp.*p1_3-G12temp.*p2_3)./sqrt(gtemp);
    mu2=(-G12temp.*p1_3+G11temp.*p2_3)./sqrt(gtemp);
    dmu=divLambda(mu1,mu2);
    ener1=sqrt(gtemp)+beta*sum(dmu.^2,3).*sqrt(gtemp);
    ener2=sum((u_3-f).^2,3)/(2*eta);
    ener3=sum(abs(p1_3)+abs(p2_3),3);
    ener=sum(ener1(:))+sum(ener2(:));
    ss=sqrt(sum(abs(u(:)-u_3(:)).^2))/sqrt(sum(abs(u(:)).^2));
    err=[err;ss];
    energy=[energy;ener];
    sprintf('%i %i %e %e',kk,Iterp,ener,ss) % print errors
    %%%%%%%%% update
    u=u_3;
    p1=p1_3;
    p2=p2_3;
    lambda1=lambda1_2;
    lambda2=lambda2_2;
    G11=G11_3;
    G12=G12_3;
    G22=G22_3;
    g=g_3;
end

