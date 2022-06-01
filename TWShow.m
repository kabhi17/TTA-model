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

% The three coefficients of the quadratic equation for the steady state concentration of annihilator triplets are given as:
a=1.25*ktt+(0.25*ktt*ksens*kex*I.*S0)./((kfl+kNR).*(kex*I+ksens.*A0+kts))+(0.75*ktt*ksens*kex*I.*S0)./((kic.*(kex*I+ksens.*A0+kts)));
b=kt+(ksens*kex.*I*S0)./(kex*I+ksens.*A0+kts);
c=-(ksens*kex*I.*S0.*A0)./(kex*I+ksens.*A0+kts);   
Fss=(0.25*FL)*((-b+sqrt((b.^2)-4*a.*c))./(2*a)).^2;

DFss=diff(Fss,I);n=(I./Fss).*DFss;eqn=n==1.1;eqn1=n==0.9;%Defines 2 equations for the slope when n = 1.1 and 0.9
TW=log10(vpasolve(eqn1(i),I))-log10(vpasolve(eqn(i),I)); %Calculates the transition width 

% For a range of different system parameters, say for kT values between
% 0.002 s^-1 and 200 s^-1, the transition width is calculated as below

TW=[] % Defines an empty vector called 'Transition Width' (TW)
for i=1:length(eqn) % For loop solves each iteration of kT value
    TW=[TW,log10(vpasolve(eqn1(i),I))-log10(vpasolve(eqn(i),I))];
end