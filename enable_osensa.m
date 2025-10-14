
function [transmitter] = enable_osensa(port)
    ospy=py.importlib.import_module('osensaMatlab');
    transmitter = ospy.Transmitter(port, uint16(247));
    disp("Osensa Transmitter ON")
end
