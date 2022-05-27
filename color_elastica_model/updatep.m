function [q1,q2]=updatep(p1,p2,s,s1,alpha1,beta,tau,tol)
IterNum=50;
scale1=1;
q1=zeros(size(p1));
q2=zeros(size(p2));

q1temp=p1;
q2temp=p2;
err=10;
qq=p1(:,:,1);
kk=0;
while err>tol && kk<IterNum
    kk=kk+1;
    q1temp(:,:,1)=qq;
    [G11,G12,G22]=Mp(q1temp,q2temp,alpha1);
    m=G11.*G22-G12.^2;
    dq1=1/tau*(qq-p1(:,:,1))+1/2*((1+s1)./sqrt(m)-beta.*s./(m.^(3/2))).*(2*G22.*qq-2*G12.*p2(:,:,1));
    dq2=1/tau+1/2*((1+s1).*m.^(-1/2)-beta*s.*m.^(-3/2)).*(2*G22-2*p2(:,:,1).^2)+...
        1/2*((1+s1).*(-1/2).*m.^(-3/2)+3/2*beta*s.*m.^(-5/2)).*(2*G22.*qq-2*G12.*p2(:,:,1)).^2;
    qq1=qq-scale1*dq1./dq2;
    err=max(abs(qq1(:)-qq(:)));
    qq=qq1;
end
q1(:,:,1)=qq1;

q1temp=p1;
q2temp=p2;
err=10;
qq=p1(:,:,2);
kk=0;
while err>tol && kk<IterNum
    kk=kk+1;
    q1temp(:,:,2)=qq;
    [G11,G12,G22]=Mp(q1temp,q2temp,alpha1);
    m=G11.*G22-G12.^2;
    dq1=1/tau*(qq-p1(:,:,2))+1/2*((1+s1)./sqrt(m)-beta.*s./(m.^(3/2))).*(2*G22.*qq-2*G12.*p2(:,:,2));
    dq2=1/tau+1/2*((1+s1).*m.^(-1/2)-beta*s.*m.^(-3/2)).*(2*G22-2*p2(:,:,2).^2)+...
        1/2*((1+s1).*(-1/2).*m.^(-3/2)+3/2*beta*s.*m.^(-5/2)).*(2*G22.*qq-2*G12.*p2(:,:,2)).^2;
    qq1=qq-scale1*dq1./dq2;
    err=max(abs(qq1(:)-qq(:)));
    qq=qq1;
end
q1(:,:,2)=qq1;

q1temp=p1;
q2temp=p2;
err=10;
qq=p1(:,:,3);
kk=0;
while err>tol && kk<IterNum
    kk=kk+1;
    q1temp(:,:,3)=qq;
    [G11,G12,G22]=Mp(q1temp,q2temp,alpha1);
    m=G11.*G22-G12.^2;
    dq1=1/tau*(qq-p1(:,:,3))+...
        1/2*((1+s1)./sqrt(m)-beta.*s./(m.^(3/2))).*(2*G22.*qq-2*G12.*p2(:,:,3));
    dq2=1/tau+1/2*((1+s1).*m.^(-1/2)-beta*s.*m.^(-3/2)).*(2*G22-2*p2(:,:,3).^2)+...
        1/2*((1+s1).*(-1/2).*m.^(-3/2)+3/2*beta*s.*m.^(-5/2)).*(2*G22.*qq-2*G12.*p2(:,:,3)).^2;
    qq1=qq-scale1*dq1./dq2;
    err=max(abs(qq1(:)-qq(:)));
    qq=qq1;
end
q1(:,:,3)=qq1;

q1temp=p1;
q2temp=p2;
err=10;
qq=p2(:,:,1);
kk=0;
while err>tol && kk<IterNum
    kk=kk+1;
    q2temp(:,:,1)=qq;
    [G11,G12,G22]=Mp(q1temp,q2temp,alpha1);
    m=G11.*G22-G12.^2;
    dq1=1/tau*(qq-p2(:,:,1))+...
        1/2*((1+s1)./sqrt(m)-beta.*s./(m.^(3/2))).*(2*G11.*qq-2*G12.*p1(:,:,1));
    dq2=1/tau+1/2*((1+s1).*m.^(-1/2)-beta*s.*m.^(-3/2)).*(2*G11-2*p1(:,:,1).^2)+...
        1/2*((1+s1).*(-1/2).*m.^(-3/2)+3/2*beta*s.*m.^(-5/2)).*(2*G11.*qq-2*G12.*p1(:,:,1)).^2;
    qq1=qq-scale1*dq1./dq2;
    err=max(abs(qq1(:)-qq(:)));
    qq=qq1;
end
q2(:,:,1)=qq1;

q1temp=p1;
q2temp=p2;
err=10;
qq=p2(:,:,2);
kk=0;
while err>tol && kk<IterNum
    kk=kk+1;
    q2temp(:,:,2)=qq;
    [G11,G12,G22]=Mp(q1temp,q2temp,alpha1);
    m=G11.*G22-G12.^2;
    dq1=1/tau*(qq-p2(:,:,2))+...
        1/2*((1+s1)./sqrt(m)-beta.*s./(m.^(3/2))).*(2*G11.*qq-2*G12.*p1(:,:,2));
    dq2=1/tau+1/2*((1+s1).*m.^(-1/2)-beta*s.*m.^(-3/2)).*(2*G11-2*p1(:,:,2).^2)+...
        1/2*((1+s1).*(-1/2).*m.^(-3/2)+3/2*beta*s.*m.^(-5/2)).*(2*G11.*qq-2*G12.*p1(:,:,2)).^2;
    qq1=qq-scale1*dq1./dq2;
    err=max(abs(qq1(:)-qq(:)));
    qq=qq1;
end
q2(:,:,2)=qq1;

q1temp=p1;
q2temp=p2;
err=10;
qq=p2(:,:,3);
kk=0;
while err>tol && kk<IterNum
    kk=kk+1;
    q2temp(:,:,3)=qq;
    [G11,G12,G22]=Mp(q1temp,q2temp,alpha1);
    m=G11.*G22-G12.^2;
    dq1=1/tau*(qq-p2(:,:,3))+...
        1/2*((1+s1)./sqrt(m)-beta.*s./(m.^(3/2))).*(2*G11.*qq-2*G12.*p1(:,:,3));
    dq2=1/tau+1/2*((1+s1).*m.^(-1/2)-beta*s.*m.^(-3/2)).*(2*G11-2*p1(:,:,3).^2)+...
        1/2*((1+s1).*(-1/2).*m.^(-3/2)+3/2*beta*s.*m.^(-5/2)).*(2*G11.*qq-2*G12.*p1(:,:,3)).^2;
    qq1=qq-scale1*dq1./dq2;
    err=max(abs(qq1(:)-qq(:)));
    qq=qq1;
end
q2(:,:,3)=qq1;