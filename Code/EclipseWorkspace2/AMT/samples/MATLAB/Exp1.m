%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sean Goldberg
% 2/13/2012
% CrowdPillar Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 1:
%
% Test the correlation between accuracy and entropy for each token
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Exp = csvread('Exp1_acc_vs_ent.csv');

numtokens = size(Exp,1);

Exp_0 = Exp(:,2);
Exp_0 = Exp_0(Exp(:,1)==0);
Exp_1 = Exp(:,2);
Exp_1 = Exp_1(Exp(:,1)==1);

[n0, out0] = hist(Exp_0,100);
[n1, out1] = hist(Exp_1,100);

figure; bar(out0,n0,'r'); grid;
figure; bar(out1,n1,'b'); grid;