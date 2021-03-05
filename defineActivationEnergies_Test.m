function [Qg,Qc,Qm] = defineActivationEnergies_Test(T,Qgm,Qgp,tc)

Temp = T-273;

% Activation Energy for Creep
tp = 0;
tm = -20;
Qcp = 100; %kJ/mol
Qcm = 60; %kJ/mol
c1 = (Qcp-Qcm)./(tanh(tp-tc)-tanh(tm-tc));
c2 = Qcp-c1*tanh(tp-tc);

% Activation Energy for Grain Growth
% define Arrhenius relation for grain growth
tp = 0;
tm = -20;
g1 = (Qgp-Qgm)./(tanh(tp-tc)-tanh(tm-tc));
g2 = Qgp-g1*tanh(tp-tc);

% Activation Energy for Grain Boundary Mobility
% define Arrhenius relation for grain boundary mobility
tp = 0;
tm = -20;
Qmp = Qgm; %kJ/mol
Qmm = Qgp; %kJ/mol
m1 = (Qmp-Qmm)./(atan(tp-tc)-atan(tm-tc));
m2 = Qmp-m1*atan(tp-tc);

Qg = zeros(size(T));
Qg = g1*tanh(Temp-tc)+g2;
Qg = Qg.*1e3;

Qc = zeros(size(T));
Qc = c1*tanh(Temp-tc)+c2;
Qc = Qc.*1e3;

% define Arrhenius relation for grain boundary mobility
Qm = zeros(size(T));
Qm = m1*atan(Temp-tc)+m2;
Qm = Qm.*1e3;

end

