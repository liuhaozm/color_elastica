function [lm1,lm2]=initialLambda(g,G11,G12,G22,p1,p2)

% [m,n]=size(g);

% lm1=zeros(m,n,3);
% lm2=zeros(m,n,3);
% 
% lm11=zeros(m,n,3);
% lm12=zeros(m,n,3);
% lm21=zeros(m,n,3);
% lm22=zeros(m,n,3);
% lm31=zeros(m,n,3);
% lm32=zeros(m,n,3);
sqg=sqrt(g);
fac=sqg./(g+1e-10);
lm1=fac.*(G22.*p1-G12.*p2);
lm2=fac.*(-G12.*p1+G11.*p2);
% lm11=fac.*(G22.*p1(:,:,1)-G12.*p2(:,:,1));
% lm12=fac.*(-G12.*p11+G11.*p12);
% lm21=fac.*(G22.*p21-G12.*p22);
% lm22=fac.*(-G12.*p21+G11.*p22);
% lm31=fac.*(G22.*p31-G12.*p32);
% lm32=fac.*(-G12.*p31+G11.*p32);