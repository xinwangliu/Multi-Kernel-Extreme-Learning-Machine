function [alpha] = AverageMKELM(K,C,Y)
numKer = size(K,3);
numTrn = size(K,1);
gamma0 = ones(numKer,1)/numKer;
KC = sumKbeta(K,gamma0);

alpha = (KC+eye(numTrn)/C)\Y;