function [dydt] = enz_kin_inh(t,y,k1,k_1,k2,k3,k_3)
E = y(1);
S = y(2);
C = y(3);
P = y(4);
CI = y(5);

% dEdt = -k1*E*S + k_1*C + k2*C;
% dSdt = -k1*E*S + k_1*C;
% dCdt = k1*E*S - k_1*C - k2*C + k_3*CI - k3*C*S;
% dPdt = k2*C;
% dCIdt = k3*C*S - k_3*CI;


dEdt = -k1*E*S + k_1*C + k2*C;
dSdt = -k1*E*S + k_1*C - k3*C*S + k_3*CI;
dCdt = k1*E*S - k_1*C - k2*C + k_3*CI - k3*C*S;
dPdt = k2*C;
dCIdt = k3*C*S - k_3*CI;


dydt = [dEdt;dSdt;dCdt;dPdt;dCIdt];
end
