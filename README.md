# HyperSym: An Educational MATLAB Code for Hyperelasticity

This repository contains **HyperSym**, a **MATLAB**-based educational tool that uses symbolic differentiation to derive stress and tangent stiffness tensors for hyperelastic materials. It is designed to support students and researchers in learning and implementing nonlinear finite element methods, and integrates directly with the NLFEA framework.

We also presente **HyperFit**, another educational code that adjusts hyperlastic material parameters to experimental test data using symbolic derivation.

This file provides a quick start, so for further information on `HyperSym`, please refer to our article in the following reference.

> 📄 *Vinicius O. Fontes, André X. Leitão, Anderson Pereira*  
> **HyperSym: An Educational MATLAB Code for Hyperelasticity**  
> *Computer Applications in Engineering Education*, 2025  
> [DOI: 10.1002/cae.70037](https://doi.org/10.1002/cae.70037)

## 🚀 Features

- Symbolic derivation of second Piola–Kirchhoff stress and material tangent tensors
- Support for compressible and nearly incompressible hyperelastic models
- Automatic generation of **MATLAB** functions for integration with FE frameworks
- Seamless integration with the educational **NLFEA** code

## 📦 Prerequisites

- **MATLAB** R2022a or later  
- **Symbolic Math Toolbox**
- A supported **C/C++ compiler** (e.g., Microsoft Visual Studio or MinGW-w64) — required *only* for generating **MEX** files

## 🛠 Usage: HyperSym

To use **HyperSym**, you must first download the source code, unpack the main folder `HyperSym-main`, and open it in **MATLAB**.

The example below show these steps for the nine hyperelastic material models presented in our article:

1. Download the source code from this reppository and extract the folder `HyperSym-main` to your working directory.
2. Open the `HyperSymScript.m` file in **MATLAB** in the `HyperSym + NLFEA/HyperSym (modified)` folder.
3. Specify the symbolic variables: deformation metrics and material parameters (`ModHyperSymScript.m` lines 27-35).

```matlab
%% Symbolic variables
syms I1 I2 J J1 J2 lambda mu A10 A20 A30 A01 K D1 D2 D3 real
% Cell array to hold the strain energy density functions
MID = {'SVK'; 'mSVK1'; 'mSVK2'; 'mSVK3'; 'nH1'; 'nH2'; 'nH3'; 'MR'; 'Yeoh'};
W = cell(length(MID),1);
PROP = cell(length(MID),1);       % Array of properties for each model
[PROP{1:7}] = deal([lambda mu]);  % SVK and nH-based models
PROP{8} = [A10 A01 K];            % MR
PROP{9} = [A10 A20 A30 D1 D2 D3]; % Yeoh 3 parameters
```

4. Define your strain energy function `W` using previously defined symbolic variables (`ModHyperSymScript.m` lines 36-38).
   - For example, line 38 shows the strain energy function for the SVK model in terms of the principal invariants of the Cauchy-Green tensor.

```matlab
%% Compressible models, i.e. W = W(I1,I2,J)
% Materials 1-4: St. Venant-Kirchhoff
W{1} = 1/8*lambda*(I1 - 3)^2 + mu/4*(I1^2 - 2*I2 - 2*I1 + 3);
```

5. Run the script to perform symbolic differentiation and generate (`ModHyperSymScript.m` lines 65-72):
   - `Sv` – second Piola–Kirchhoff stress tensor (Voigt vector)
   - `Dv` – material tangent stiffness matrix (Voigt matrix)
  
```matlab
[Sv,Dv,C] = HyperSym(W{m});        
% Store tensors as column vector (for MEX compatibility)
Dv = Dv(:); 
% Generate Matlab functions
Controls.File = [Dir MID{m}];
matlabFunction(Sv,Dv,'Vars',{C,PROP{m}},Controls);
% Generate MEX function
if generateMEX,Generate_mex_models(MID{m},PROP{m},false); end
```
 
6. The generated functions will be saved in the `HyperSym (modified)/Functions/` directory and can be directly used in FE simulations.
   - If your call the function `Generate_mex_models`, both **MATLAB** (`.m`) and **MEX** (`.mexw64`) functions will be generated in the directory.

![image](https://github.com/user-attachments/assets/a9b3975d-4533-479e-98a8-aa8f80ef5e27)

   - If you defined a nearly incompressible hyperelastic model, isochoric (`_iso`) and volumetric (`_vol`) functions will also be created.

![image](https://github.com/user-attachments/assets/19d7c7ea-1e75-4e3b-ae35-a405bd06e6aa)

## 🔗 HyperSym Integration with NLFEA

The scripts in the `HyperSym-main` folder are integrated with the educational finite element code, **NLFEA**. The original code for **NLFEA** is available at https://web.mae.ufl.edu/nkim/INFEM/.

When running one of the example simulations described in the next section, **NLFEA** invokes functions generated by **HyperSym**. These functions are located in the `HyperSym (modified)/Functions` directory. The preferred format for these functions is **MEX** (`_.mexw64`), provided they are available.

As a result, an `Output` folder is automatically created, where the results from the example simulations are saved.

## 📊 Examples: HyperSym

Example problems included in this repository:

1. Cube under uniaxial traction: Compare the response of several compressible hyperelastic models in tension.

<p align="center">
  <img src="https://github.com/user-attachments/assets/2a205fd7-412b-47a7-aefc-6f6a0077756d" height="250px" />
  <img src="https://github.com/user-attachments/assets/a51462e7-6eb1-4323-a99c-325bea211c3a" height="250px" />  
</p>


2. Thin square strip: Demonstrates large deformation and sensitivity to compressibility. There are two scripts.

- `Strip.m`: Visualize the deformation using the MR model for different mesh sizes, and compare it to solutions obtained with **ANSYS**.

<p align="center">
  <img src="https://github.com/user-attachments/assets/15328a70-ee54-4d84-b497-9845af23d809" height="250px" />
  <img src="https://github.com/user-attachments/assets/eb0c90fa-cbea-49f3-a583-81451c0eb7c6" height="250px" />
</p>

- `Strip_models.m`: Compare the response of several compressible hyperelastic models for the same mesh discretization.

<p align="center">
  <img src="https://github.com/user-attachments/assets/9cfc9b40-6013-4728-ad1a-202ebe1e7e23" height="250px" />
  <img src="https://github.com/user-attachments/assets/c61c1baa-83f1-47ac-a554-f621e4df8380" height="250px" />
</p>

3. Thick-walled cylinder: Compare the deformation of a pressure vessel simulation to the results from ANSYS.

<p align="center">
  <img src="https://github.com/user-attachments/assets/19026140-a8c9-4fb0-bc5e-203449894c9d" height="250px" />
  <img src="https://github.com/user-attachments/assets/eb2b8b6a-6f30-4cac-b839-685d54baa83e" height="250px" />
</p>

## 🛠 Usage: HyperFit

To use **HyperFit**, you must first download the source code, unpack the main folder `HyperSym-main`, and open it in **MATLAB**.

The example below show these steps for fitting Treloar's data to Ogden hyperelastic model:

1. Download the source code from this reppository and extract the folder `HyperSym-main` to your working directory.
2. Open the `HyperFitScript_Ogden.m` file in **MATLAB** in the `HyperFit` folder.
3. For each experimental test, add the test data to a two-column array and store it in the ``Test`` struct (`HyperFitScript_Ogden.m` lines 23-28).
   - For example, the uniaxial stress data is input as follows.
  
```matlab
% Digitalized data from [1], with stress converted to MPa.
% ua: Simple Elongation (uniaxial stress)
ua_exp = readmatrix('Treloar Data/Treloar_uniaxial.csv');
ua_exp(:,2) = ua_exp(:,2)*kg2MPa; % Converting units
Test(1) = struct(...
    'name','uniaxial',...
    'data',ua_exp);
```

4. Define the hyperelastic material model as a function of symbolic variables (`HyperFitScript_Ogden.m` lines 41-51).
   - The strain energy density must be function of the principal stretches `lambda1`, `lambda2` and `lambda3`.
   - Material parameters must be added to the array `PROP`.

```matlab
%% Define the material model: Ogden
syms lambda1 lambda2 lambda3
M = 4; % Number of pairs alpha and mu
W = sym(0);                           
alpha = sym('alpha',[M 1],'real');    mu = sym('mu',[M 1],'real');
for i = 1:M
    W = W + (mu(i)/alpha(i))*(...
        lambda1^alpha(i) + lambda2^alpha(i) + lambda3^alpha(i) - 3);
end
% Array of properties
PROP = reshape([alpha,mu]',1,2*M); % PROP = [alpha1, mu1, alpha2, mu2, ...]
```

5. Define optimization parameters for the curve-fitting algorithm and call **HyperFit** (`HyperFitScript_Ogden.m` lines 52-61).

```matlab
%% Call HyperFit
opt = struct(...
    'N',5,...     % Number of optimization runs
    'pL',-10,...  % Lower boundary of p
    'pU',+10,...  % Upper boundary of p
    'p0',[],...   % First initial guess
    'cc',true ... % Check for consistency conditions
    );
% HyperFit returns fitted parameters (p), and error metrics (Error, res)
[p,Error,res] = HyperFit(W,PROP,Test,opt); 
```

**Important:** A function must be created for each type of experimental test used in the curve fitting algorithm. For instance, for tests with the uniaxial stress condition, we use the following (`HyperFit` lines 103-109):

```matlab
%% Uniaxial test
function t1 = uniaxial(W)
% Kinematic constraints for each uniaxial test.
syms lambda1 lambda2
W = subs(W,lambda2,1/sqrt(lambda1)); % Apply kinematic constraint
t1 = simplify(diff(W,lambda1));       % t = dW/dlambda1
end
```

Equibiaxial stress and pure shear functions are already included. The name of the test/condition must match that in the `Test(i).name` string (where `i` denotes the `i`-th test).

Upon running `HyperFitScript_Ogden.m`, the curve fitting algorithm is run 5 times, and the optimal parameters are output to the array `p_opt`.

<div align="center">
  <img src="https://github.com/user-attachments/assets/22d5a76b-a846-49eb-a012-87c0cc1a1f36" height="250px" />
</div>

Additionally, graphs comparing experimental and fitted data are plotted/

<p align="center">  
  <img src="https://github.com/user-attachments/assets/5506d6c4-7cb3-4294-95cf-a615a79c5421" height="250px" />
  <img src="https://github.com/user-attachments/assets/6962fc38-7e41-4f99-9fb9-a601dabac920" height="250px" />
</p>

Additional examples using MR and nH models are provided in the scripts `HyperFitScript_MR.m` and `HyperFitScript_nH.m`, respectively.

## 📖 Citation
If you use this software in academic work, please cite:

```bibtex
@article{fontes2025hypersym,
  title   = {HyperSym: An Educational MATLAB Code for Hyperelasticity},
  author  = {Fontes, Vinicius O. and Leitão, André X. and Pereira, Anderson},
  journal = {Computer Applications in Engineering Education},
  year    = {2025},
  doi     = {10.1002/cae.70037}
}
```

## 🏛️ Funding

Funding: This study was financed by Carlos Chagas Filho Foundation for Research Support of the State of Rio de Janeiro (FAPERJ), Process SEI 260003/008107/2022 and E‐26/202.573/2022 (277082), the Coordination for the Improvement of Higher Education Personnel – Brazil (CAPES) – Finance Code 001, and National Council for Scientific and Technological Development – Brazil (CNPq) (grant number 317319/2021‐3).

The Article Processing Charge for the publication of this research was funded by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior ‐ Brasil (CAPES) (ROR identifier: 00x0ma614).

## 🙏 Acknowledgemets

The authors kindly thank Professor Nam‐Ho Kim for authorizing the use and modification of his finite element source code **MATLAB** code, NLFEA, and Professor Raymond Ogden for sharing the experimental data used for the curve fitting algorithm. 

## 📜 Version History

- v1.0 (April 24, 2025): README added.
- v1.1 (April 25, 2025): HyperFit files added to repository and README. DOI added to the header of **MATLAB** files.
