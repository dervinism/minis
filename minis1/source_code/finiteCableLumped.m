function [V, riseTime, t10, t90, t100] = finiteCableLumped(nmax, L, lambda, D, C_m, t, tau_m, X, Q, tau_sy1, tau_sy2)
% finiteCableLumped simulates a single shape of a mini of 1mV amplitude
% (standardised). The function uses a finite cable and implements lumped
% terms.
%
% The input variables are:
%
% nmax - a number of exponentials;
% L - electrotonic length of a cylinder;
% lambda - the dendritic length constant;
% D - the diameter of a dendritic cylinder, cm;
% C_m - the membrane capacitance;
% t - the duration of a simulated minis (10 taus);
% tau_m - the dendritic time constant, ms;
% X - the location of a point charge;
% Q - the input charge;
% tau_sy1 - the synaptic rise constant;
% tau_sy2 - the synaptic decay constant.
%
% The output variables are:
%
% V - the membrane potential of the soma across the simulation time;
% riseTime - the 10-90% rise time of a simulated minis

% Fourier coefficients:
c_m = C_m*pi*D; % capacitance per unit length
g_inf = c_m*lambda/tau_m; % keeps units consistent at any rate! ie ms vs s
B_0 = 1000*Q/(L*lambda*c_m); % mV
B = 2*B_0;

% Membrane potential for given L and X:
if tau_sy1==tau_m
    disp('tau_sy1 == tau_m')
    tau_sy1=0.95*tau_sy1;
end
n = 0;
A_0 = B_0*cos((n*pi*X)/L);
V = A_0*exp(-t./tau_m)*tau_m^2/( (tau_m-tau_sy1)*(tau_m-tau_sy2) );
for n = 1:nmax
    tau_n = tau_m/(1+((n*pi)/L)^2);
    if tau_n < t(1)/4
        break;
    end
    A_n=B*cos((n*pi*X)/L);
    V_n = A_n*exp(-t./tau_n)*tau_n^2/( (tau_n-tau_sy1)*(tau_n-tau_sy2) );
    V = V + V_n;
end

A_lump1 = Alump(tau_m,tau_sy1,L,X,g_inf);
A_lump1= -1000*Q*A_lump1/(tau_sy1-tau_sy2); % since convention of biophys J paper reversed here: tau_sy1 is the slowest

A_lump2 = Alump(tau_m,tau_sy2,L,X,g_inf);
A_lump2= -1000*Q*A_lump2/(tau_sy2-tau_sy1);

V_lump = A_lump1*exp(-t./tau_sy1) + A_lump2*exp(-t./tau_sy2);
V = [0 V + V_lump];

[Vmax, maxindex] = max(V);
V = V/Vmax;

% Calculating 10-90% amplitude rise time:
V10 = V - 0.1;
V90 = V - 0.9;
[~, arrayPosition10] = min(abs(V10(1:maxindex)));
[~, arrayPosition90] = min(abs(V90(1:maxindex)));
t = [0 t];
riseTime = t(arrayPosition90) - t(arrayPosition10);
t10 = t(arrayPosition10);
t90 = t(arrayPosition90);
t100 = t(maxindex);
end

function A_lump = Alump(tau_m,tau_sy,L,X,g_inf)
w=sqrt( (tau_m/tau_sy) - 1);
A_lump = cos(w*(L-X))/(w*g_inf*sin(w*L));
end