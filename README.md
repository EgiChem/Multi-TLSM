# Multi Tracer Liu-Silva-Macedo to calculate diffusivities in ternary systems

Calculates the tracer diffusivity of a given solute in a mixture of 2 solvents. If used please cite:"xxx"

Requiered data to use:

* `T` - Absolute temperature in K in column format.
* `Density` - Density of the solvent in g/cucm  in column format.
* `x2` - Molar fraction of the componente 2 (solvent 2)  in column format.
* `M` - Molecular mass of the compounds by the following order solute, solvent 1, solvent 2 in vector or column format.
* `dLJ` - Molecular LJ diameter of the compounds in angstrom by the following order solute, solvent 1, solvent 2 in vector or column format.
* `eLJ` - Molecular LJ energy of the compounds in K by the following order solute, solvent 1, solvent 2 in vector or column format.
* `AD` - (optional) adjustable parameter of the model, if unknown it is set to 1.


## Examples
Copy paste on Matlab to run, tested in Matlab 2021b

### Example 1 benzoic acid in CO2/methanol using Multi-TLSM

```matlab
T = 308; %K
Density = 0.83949; % g/cum
x2 = 0.05; % mol/mol fraction of solvent 2

% sequence is solute; solvent_1; solvent_2

M = [122.124; 44.01; 32.042]; % g/mol
dLJ = [5.65763; 3.26192; 3.79957];  % angstrom
eLJ = [582.05; 500.71; 685.96]; % K;

[D12calc_Multi_LSM] = Muli_LSM(T, Density, x2, M, dLJ, eLJ) % sqrtcm/s
```


### Example 2 benzoic acid in CO2/methanol using Multi-TLSM_AD

```matlab
T = [308; 318; 328]; %K

Density = [0.83949; 0.78356; 0.72005]; % g/cum

x2 = [0.05; 0.05; 0.05]; % mol/mol fraction of solvent 2

M = [122.12; 44.01; 32.042]; % g/mol

%Function used when some of the LJ constants are unknown

[dlj_solute, elj_solute] = Lennard_Jones(752.00, 45.60, 341); %Tc in K;  Pc in bar; Vc in cm3 mol-1

% sequence is solute; solvent_1; solvent_2

dLJ = [dlj_solute; 3.26192; 3.79957];  % angstrom

eLJ = [elj_solute; 500.71; 685.96]; % K

AD = 0.8674; 

[D12calc_Multi_LSM_AD] = Muli_LSM(T, Density, x2, M, dLJ, eLJ, AD) % sqrtcm/s
```
