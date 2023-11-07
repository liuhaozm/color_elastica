function [q1,q2,kk]=updatepSecond(p1,p2,s,alpha1,beta,tau,tol)
% fix point method to update p in the second modified model

IterNum=200;



q1=p1;
q2=p2;

qq1=q1;
qq2=q2;

err=10;
qq=p1(:,:,1);
kk=0;
while err>tol && kk<IterNum
    kk=kk+1;
    %     q1temp(:,:,1)=qq;
    [G11,G12,G22]=Mp(q1,q2,alpha1);
    m=G11.*G22-G12.^2-alpha1^2;
    
    ww=(1+beta*s)./(sqrt(m)+1e-3);
    
    nomi=p1(:,:,1)/tau+G12.*q2(:,:,1).*ww;
    denomi=1/tau+G22.*ww;
    qq1(:,:,1)=nomi./denomi;
    
    nomi=p1(:,:,2)/tau+G12.*q2(:,:,2).*ww;
    denomi=1/tau+G22.*ww;
    qq1(:,:,2)=nomi./denomi;
    
    nomi=p1(:,:,3)/tau+G12.*q2(:,:,3).*ww;
    denomi=1/tau+G22.*ww;
    qq1(:,:,3)=nomi./denomi;
    
    nomi=p2(:,:,1)/tau+G12.*q1(:,:,1).*ww;
    denomi=1/tau+G11.*ww;
    qq2(:,:,1)=nomi./denomi;
    
    nomi=p2(:,:,2)/tau+G12.*q1(:,:,2).*ww;
    denomi=1/tau+G11.*ww;
    qq2(:,:,2)=nomi./denomi;
    
    nomi=p2(:,:,3)/tau+G12.*q1(:,:,3).*ww;
    denomi=1/tau+G11.*ww;
    qq2(:,:,3)=nomi./denomi;
    
    err1=max(abs(qq1(:)-q1(:)));
    err2=max(abs(qq2(:)-q2(:)));
    err=max(err1,err2);
    q1=qq1;
    q2=qq2;
end







