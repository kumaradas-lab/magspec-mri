function J_total = get_TotalMagneticPolarization(HW, mol_atoms, T_C, B0, gamma, I)
%% Calculate total magnetic polarization of completely relaxed spin system
%
%   J_total = get_TotalMagneticPolarization(HW, mol_atoms, T_C, B0, gamma, I)
%
% Calculate the total magnetic polarization of a spin system following basic
% thermodynamic equations.
%
%
% INPUT:
%
%   If the function is called with fewer input arguments or they are empty,
%   default values might be used.
%
%   HW
%       HW object or structure
%
%   mol_atoms
%       Total number of atoms in mol.
%       (Default: 1)
%
%   T_C
%       Thermodynamic equilibrium temperature in degrees Celsius.
%       (Default: 30.1)
%
%   B0
%       External magnetic field in Tesla.
%       (Default: HW.B0)
%
%   gamma
%       Gyromagnetic ratio of the atoms in consideration in rad/s/T.
%       (Default: HW.GammaDef)
%
%   I
%       Nuclear spin quantum number of the atoms in consideration.
%       (Default: 1/2)
%
%
% OUTPUT:
%
%   J_total
%       Expected total magnetic polarization in T*m^3 at the conditions
%       according to the input parameters.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input parameters
if (nargin <2) || isempty(mol_atoms)
  mol_atoms = 1;
end
if (nargin <3) || isempty(T_C)
  T_C = 30.1;
end
if (nargin <4) || isempty(B0)
  B0 = HW.B0;
end
if (nargin <5) || isempty(gamma)
  gamma = HW.GammaDef;
end
if (nargin <6) || isempty(I)
  I = 1/2;
end


%% magnetic polarization at expected population inversion
num_atoms = mol_atoms .* HW.Constant.Avogadro;
x = gamma .* HW.Constant.h_Planck .* B0 ./ HW.Constant.Boltzmann ./ (T_C + HW.Constant.TZero);
if I == 1/2
  % Taylor series expansion (up to 6-th order) for spin 1/2
  % NOTE: The contributions of the terms of 5-th or higher order are below the
  %       limit of double precision floating point numbers for reasonable
  %       equilibrium temperatures.
  a = 0;
  for n = 1:2:6
    a = a + x.^n./(2.^n .* factorial(n));
  end
  b = 2;
  for n = 2:2:6
    b = b + x.^n./(2.^(n-1) .* factorial(n));
  end
else
  % FIXME: The following loop effectively calculates very small differences
  %        between close floating point numbers. That could potentially lead to
  %        large numeric deviations. Consider using an approximation (e.g.,
  %        Taylor series expansion). See, e.g., the case for spin 1/2 above.
  a = 0;
  b = 0;
  for Im = -I:1:I
    a = a + Im .* exp(Im .* x);
    b = b + exp(Im .* x);
  end
end
J_total = HW.Constant.My0 .* num_atoms .* gamma .* HW.Constant.h_Planck .* a ./ b;

end

% AtomesPerVolume=Thermometer.mol_1H(:).*Avogadro./WaterVolume;
%
% I=1/2;
% % % M0=AtomesPerVolume*HW.GammaDef^2*(HW.Constant.Planck/2/pi)^2*I*(I+1)/3/HW.Constant.Boltzmann/(T+HW.Constant.TZero)*HW.B0
% % M0=AtomesPerVolume .* HW.B0 ./ (T+HW.Constant.TZero) .* HW.GammaDef.^2 .* (I.^2+I)./3 .* HW.Constant.h_Planck.^2 ./ HW.Constant.Boltzmann ; % in A/m
% % J=M0*HW.Constant.My0; % in T
%
% x=HW.GammaDef.*HW.Constant.h_Planck.* HW.B0./HW.Constant.Boltzmann./(T+HW.Constant.TZero);
% a=0;
% b=0;
% for Im=-I:1:I
%   a=a+Im.*exp(Im.*x);
%   b=b+exp(Im.*x);
% end
% M0=AtomesPerVolume.*HW.GammaDef.*HW.Constant.h_Planck.*a./b;
% % M0=AtomesPerVolume*HW.GammaDef*HW.Constant.h_Planck.*(-0.5*exp(-0.5*x)+0.5*exp(0.5*x))/(exp(-0.5*x)+exp(0.5*x)) % I=0.5
%
% MM=M0.*WaterVolume; % in A*m^2
%
% J=M0.*HW.Constant.My0; % in T
%
% J_Total=J.*WaterVolume; % in T*m^3
