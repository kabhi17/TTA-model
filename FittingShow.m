AD=readmatrix('excel/notepad/csv/ file');
xaxis=AD(:,1); %Chooses column of irradiance values
yaxis=AD(:,2)/(max(AD(:,2))); %Chooses column of fluorescence values and normalizes them
Weights=1./yaxis; %Default weights are defined as 1/normalized fluorescence emission values
syms kfl ktt kt kNR ksens kex I %Here all system parameters save for the ground state concentrations are presumed unknown. If any are known, they should be removed from this line
A0=; %Ground state annihilator concentration
S0=; %Ground state sensitizer concentration
kic=vpa(2e8); %Internal conversion and sensitizer triplet decay constants can also be set as unknown and fit for
kts=2e3; %Sensitizer triplet quenching term. If known can be provided here or else defined as a symbol earlier to solve for during fitting

a=1.25*ktt+(0.25*ktt*ksens*kex*I.*S0)./((kfl+kNR).*(kex*I+ksens.*A0+kts))+(0.75*ktt*ksens*kex*I.*S0)./((kic.*(kex*I+ksens.*A0+kts)));
b=kt+(ksens*kex.*I*S0)./(kex*I+ksens.*A0+kts);
c=-(ksens*kex*I.*S0.*A0)./(kex*I+ksens.*A0+kts);
Fss=(((0.25*ktt*kfl)/(kfl+kNR))*(((-b+sqrt((b.^2)-4*a.*c))./(2*a)).^2)); %Generates an expression for Fss which will be used in fitting