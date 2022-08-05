function [dlj, elj] = Lennard_Jones(Tc, Pc, Vc)
% Calculates the molecular LJ diameter (dlj) in A and energy (eLJ) in K
% TC - critical temperature in K
% PC - critical pressure in bar
% VC - critical volume in cucm/mol

if Tc/Pc <= 100
        dlj = (0.17791 + 11.779*Tc/Pc - 0.049029*(Tc/Pc)^2).^(1/3); %angstrom
elseif Tc/Pc > 100 && sum(isnan(Vc)) == 0
        dlj = 0.809*Vc^(1/3); % angstrom
end

elj = 0.774*Tc; %K