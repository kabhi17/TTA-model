% Defining the rate constants involved in TTA-UC 
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

syms I % Expressions for the TTA-UC quantum yield are calculated in terms of I

a=1.25*ktt+(0.25*ktt*ksens*kex*I.*S0)./((kfl+kNR).*(kex*I+ksens.*A0+kts))+(0.75*ktt*ksens*kex*I.*S0)./((kic.*(kex*I+ksens.*A0+kts)));
b=kt+(ksens*kex.*I*S0)./(kex*I+ksens.*A0+kts);
c=-(ksens*kex*I.*S0.*A0)./(kex*I+ksens.*A0+kts);   
Fss=(0.25*FL)*((-b+sqrt((b.^2)-4*a.*c))./(2*a)).^2;

Neweffic=Fss./(((ksens+kts)*AC*kex*I*S)./(kex*I*S+ksens*AC+kts)); %Calculates an expression for TTA-UC QY
nEff=(I./Neweffic).*(diff(Neweffic,I)); %Calculates an expression for the slope of the QY curve
DFss=diff(Fss,I);n=(I./Fss).*DFss;eqn=n==1; %Defines an equation to find I for when n = 1.
MaxNeweffic=subs(Neweffic,I,vpasolve(eqn,I)); %Calculates maximum TTA-UC QY 

EffArr=subs(Neweffic,I,10.^(-2:0.01:8)); %Generates an array of QY values for the system between irradiances of 0.01 and 10^8 