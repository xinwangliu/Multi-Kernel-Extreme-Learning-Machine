function [alpha] = multiclassELM(K,C,Y)

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
%%%% alpha should be numTrn*numClass
alpha = (K+eye(numTrn)/C)\YF;