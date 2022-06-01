ktt=vpa(3.6e8);
kt=vpa(2e2);
ksens=vpa(1.63e9);
kfl=vpa(1.836e8);
kNR=vpa(5.037e5);
kex=vpa(5);
S0=1.4e-4;
A0=0.1;
FL=vpa((ktt*kfl)/((kfl+kNR))); % Coefficient that precedes the steady state annihilator triplet concentration term
epi=vpa(1.4);
kic=vpa(2e8);
kts=2e3;

% Calculating steady state fluorescence rate for the given rate constants
i0=-2:0.001:8; % Defines the range of irradiances to be substituted in logarithmic units
syms I % Defines a symbolic variable representing irradiance. All other rate constants are predefined.

% The three coefficients of the quadratic equation for the steady state concentration of annihilator triplets are given as:
a=1.25*ktt+(0.25*ktt*ksens*kex*I.*S0)./((kfl+kNR).*(kex*I+ksens.*A0+kts))+(0.75*ktt*ksens*kex*I.*S0)./((kic.*(kex*I+ksens.*A0+kts)));
b=kt+(ksens*kex.*I*S0)./(kex*I+ksens.*A0+kts);
c=-(ksens*kex*I.*S0.*A0)./(kex*I+ksens.*A0+kts);   
Fss=(0.25*FL)*((-b+sqrt((b.^2)-4*a.*c))./(2*a)).^2;

% The local slope is calculated as irradiance divided by the steady-state rate of fluorescence at that irradiance, multiplied by the slope of the
% linear fluorescence vs. irradiance curve

N=1 %Value of n for which we want to find irradiance. Example: We would like to know at what irradiance does the slope reach 1 

DFss=diff(Fss,I);n=(I./Fss).*DFss; %Symbollicaly performes differentiation of Fss w.r.t to I and calculates n
P=(subs(n,I,10.^i0)); % substitues a range of values for I into the expression computed above
plot(i0,P) % Plots a graph of Log(I) vs. n
eqn=n==N; % Establishes an equation for n, where N is an value for n between 2 and 0
vpasolve(eqn,I); % Solves the established equation for I 
