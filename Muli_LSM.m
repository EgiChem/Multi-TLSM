function [D12calc] = Muli_LSM(T, Density, x2, M, dLJ, eLJ, AD)
% FOR TERNARY SYSTEMS ONLY
% Calculates the tracer diffusivity of a given solute in a mixture of 2 solvents.
% T - Absolute temperature in K in column format.
% Density - Density of the solvent in g/cucm  in column format.
% x2 - Molar fraction of the componente 2 (solvent 2)  in column format.
% M - Molecular mass of the compounds by the following order solute, solvent 1, solvent 2 in vector or column format.
% dLJ - Molecular LJ diameter of the compounds in angstrom by the following order solute, solvent 1, solvent 2 in vector or column format.
% eLJ - Molecular LJ energy of the compounds in K by the following order solute, solvent 1, solvent 2 in vector or column format.
% AD - (optional) adjustable parameter of the model, if unknown it is set to 1.

% Check if AD exists
if nargin==6
    AD=1;
elseif nargin<6
    error('Not enough input arguments')
end

% Check input consistency
Tsize = size(T);
Dsensitysize = size(Density);
x2size = size(x2);

if Tsize(2)~=1 || Dsensitysize(2)~=1 || x2size(2)~=1
    error('T, density and x2 should be column')
elseif Tsize(1)~=Dsensitysize(1) || Tsize(1)~=x2size(1) || Dsensitysize(1)~=x2size(1)
    error('T, density and x2 column size not consistent')
end

numM = numel(M);
numdLJ = numel(dLJ);
numeLJ = numel(eLJ);
numAD = numel(AD);

if numM~=numdLJ || numM~=numeLJ || numeLJ~=numM
    error('Inconsistent size of M, LJ diamenter and LJ energy')
end

if numAD~=1
    error('AD should only have 1 value')
end
    
% constants
Rg=8.3144; 
Na=6.022e23; %1/mol

    fi1 = 1.2588;
    fi2 = 0.27862;
    fi3 = 0.75; 
    fi4 = 1/1.3229;


%% input organization
data.T = T; %K
data.Density = Density;
data.x2 = x2;

solute.M = M(1); 
solute.dlj = dLJ(1)*10^-8;
solute.elj = eLJ(1);

solvent_1.M = M(2);
solvent_1.dlj = dLJ(2)*10^-8;
solvent_1.elj = eLJ(2);

solvent_2.M = M(3);
solvent_2.dlj = dLJ(3)*10^-8;
solvent_2.elj = eLJ(3);



%% calculation
dlj_1s = (solute.dlj + solvent_1.dlj)/2;
dlj_2s = (solute.dlj + solvent_2.dlj)/2;
dlj_12 = (solvent_1.dlj + solvent_2.dlj)/2;

elj_1s = sqrt((solvent_1.elj*solvent_1.dlj^3 * solute.elj*solute.dlj^3))/(dlj_1s^3);
elj_2s = sqrt((solvent_2.elj*solvent_2.dlj^3 * solute.elj*solute.dlj^3))/(dlj_2s^3);
elj_12 = sqrt((solvent_1.elj*solvent_1.dlj^3 * solvent_2.elj*solvent_2.dlj^3))/(dlj_12^3);

Teff_1s = data.T./elj_1s;
Teff_2s = data.T./elj_2s;
Teff_12 = data.T./elj_12; 

Teff_1 = data.T./solvent_1.elj;
Teff_2 = data.T./solvent_2.elj;

T_mix = (1-data.x2).*Teff_1s + data.x2.*Teff_2s;

dlj_1s_eff = 2^(1/6)*dlj_1s.*(1+(Teff_1s/fi4).^0.5).^(-1/6);
dlj_2s_eff = 2^(1/6)*dlj_2s.*(1+(Teff_2s/fi4).^0.5).^(-1/6);
dlj_12_eff = 2^(1/6)*dlj_12.*(1+(Teff_12/fi4).^0.5).^(-1/6);

dlj_1_eff = 2^(1/6)*solvent_1.dlj.*(1+(Teff_1/fi4).^0.5).^(-1/6);
dlj_2_eff = 2^(1/6)*solvent_2.dlj.*(1+(Teff_2/fi4).^0.5).^(-1/6);

dlj_mix_a = (1-data.x2).^2.*dlj_1_eff.^3 + 2*data.x2.*(1-data.x2).*dlj_12_eff.^3 + data.x2.^2.*dlj_2_eff.^3;
dlj_mix_2 = (1-data.x2).*dlj_1s_eff.^2 + data.x2.*dlj_2s_eff.^2;

M_solvent_med = data.x2.*solvent_2.M + (1-data.x2).*solvent_1.M;
M_eff = 2*(M_solvent_med * solute.M)./(M_solvent_med + solute.M);

density_d = data.Density.*Na./M_solvent_med;
density_r = density_d.*dlj_mix_a;

D12calc = AD.*21.16./(density_d.*dlj_mix_2) .*(1000*Rg*data.T./M_eff).^0.5 .*exp(-fi3.*density_r./(fi1-density_r)-fi2./T_mix);
