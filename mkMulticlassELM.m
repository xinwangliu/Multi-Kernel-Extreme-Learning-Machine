function [gamma,alpha,obj] = mkMulticlassELM(K,C,Y,qnorm)

numKer = size(K,3);
numTrn = size(K,1);
classIndx = unique(Y);
numClass = length(classIndx); 
%%% YF coding class information.
YF = zeros(numTrn,numClass);
for i =1:numTrn
    for ic =1:numClass
        if Y(i)==classIndx(ic)
            YF(i,ic) = 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma0 = ones(numKer,1)*numKer^(-1/qnorm);
flag =1;
iter =1;
objold = inf;
maxIter = 500;
while flag
    KC = sumKbeta(K,gamma0);
    %%%%%% alpha step %%%%%%%%%%%
    alpha = (KC+eye(numTrn)/C)\YF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp0 = 0;
    for ic =1:numClass
         tmp0 = tmp0 + 0.5*alpha(:,ic)'*(KC+eye(numTrn)/C)*alpha(:,ic);
    end
    obj(iter) = tmp0;
    tmp = zeros(numKer,numClass);
    for ik =1:numKer
        for ic =1:numClass
            tmp(ik,ic) = gamma0(ik)^2*(alpha(:,ic)'*K(:,:,ik)*alpha(:,ic));
        end
    end
    tmp1 = sum(tmp,2);
    gamma = (tmp1.^(1/(1+qnorm)))/(sum(tmp1.^(qnorm/(1+qnorm))))^(1/qnorm);
    if (max(abs(gamma-gamma0))<1e-4)|| (abs(obj(iter)-objold)/obj(iter))<1e-4 || (iter>maxIter)
        flag =0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma0 = gamma;
    objold = obj(iter);
    iter = iter +1;
end