function [alpha] = AveragemulticlassMKELM(K,C,Y)

classIndx = unique(Y);
numClass = length(classIndx); 
numTrn = size(K,1);
%%% YF coding class information.
YF = zeros(numTrn,numClass);
for i =1:numTrn
    for c =1:numClass
        if Y(i)==classIndx(c)
            YF(i,c) = 1;
        end
    end
end
numKer = size(K,3);
gamma0 = ones(numKer,1)/numKer;
KC = sumKbeta(K,gamma0);
%%%% alpha should be numTrn*numClass
alpha = (KC+eye(numTrn)/C)\YF;