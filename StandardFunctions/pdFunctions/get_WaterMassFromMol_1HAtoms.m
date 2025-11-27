function waterMass_g = get_WaterMassFromMol_1HAtoms(mol_1HAtoms)
%% Get water mass in gram from the number of 1H atoms in mol
%
%   waterMass_g = get_WaterMassFromMol_1HAtoms(mol_1HAtoms)
%
% The number of 1H atoms can be derived from the measurement signal amplitude
% (see get_CalibrationSampleTotalMagneticPolarization). Assuming that the
% measurement sample consists of water, this function can be used to calculate
% the mass of that water sample from the number of 1H atoms in it.
%
%
% INPUT:
%
%   mol_1Hatoms
%       Number of 1H atoms in mol.
%
%
% OUTPUT:
%
%   waterMass_g
%       Mass of water for the given number of 1H atoms in gram.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


mmH2O = 18.01528;  % molar mass of water in g/mol
% mm1H = 1.00784;  % molar mass of hydrogen atoms in g/mol
HtoD = 6400;  % 1H to deuterium 6.400

% molH2O = WaterMass_g / mmH2O;  % moles of 1H atoms
% mol_1HAtoms=molH2O*2*(1-1/HtoD);  % moles of 1H atoms
molH2O = mol_1HAtoms / (2 * (1 - 1/HtoD));  % moles of 1H atoms
waterMass_g = molH2O * mmH2O;  % moles of 1H atoms

end
