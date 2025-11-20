function mol_1HAtoms = get_mol_1HAtomsFromWaterMass(waterMass_g)
%% Get number of 1H atoms in mol from the water mass in gram
%
%   mol_1HAtoms = get_mol_1HAtomsFromWaterMass(waterMass_g)
%
% The measurement signal amplitude can be predicted by the number of 1H atoms in
% the sample (see get_CalibrationSampleTotalMagneticPolarization). Assuming that
% the measurement sample consists of water, this function can be used to
% calculate the number of 1H atoms from the mass of that water sample.
%
%
% INPUT:
%
%   waterMass_g
%       Mass of water in gram.
%
%
% OUTPUT:
%
%   mol_1Hatoms
%       Number of 1H atoms in mol for the given mass.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


mmH2O = 18.01528;  % molar mass of water in g/mol
% mm1H = 1.00784;  % molar mass of hydrogen atoms in g/mol
HtoD = 6400;  % 1H to deuterium 6.400

molH2O = waterMass_g / mmH2O;  % moles of 1H atoms
mol_1HAtoms = molH2O * 2 * (1 - 1/HtoD);  % moles of 1H atoms

end
