function dLambda=divLambda(lambda1,lambda2)

lambda1ex=expandf3(lambda1);
lambda2ex=expandf3(lambda2);

dLambda=(lambda1ex(2:end-1,2:end-1,:)-lambda1ex(1:end-2,2:end-1,:)+...
    lambda2ex(2:end-1,2:end-1,:)-lambda2ex(2:end-1,1:end-2,:));