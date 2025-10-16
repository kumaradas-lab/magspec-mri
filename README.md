# Magspec-MRI

This respository contains MATLAB code to interface with the Osensa temperature sensor, the MagSpec MRI system and process the acquired data. Before you run the main script, `mrtherm_averaged_region.m`, you will have to:

- [configure your Matlab system to use Python](https://www.mathworks.com/help/matlab/matlab-and-python.html).
- Discover the COM port of the Osensa device and set it in `mrtherm_averaged_region.m` file.
- One can, as a first setp, test Osensa alone by runing `test_osensa.m` file.
