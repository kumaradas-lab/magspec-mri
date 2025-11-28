%% Demo Sequence "DigitalIO"
% This demo sequence demonstrates how to use DigitalIOs

%% DigitalIO
    LoadSystem                                      % Load System parameters

    Seq.tRep = 0.5;                                 % repetition time
    Seq.plotSeq = 0;                                % display pulse program

    Seq.DigitalIO.SetTime = [  0e-6; ...            % times for IO events in s
                              10e-6; ...
                              20e-3; ...
                              30e-3; ...
                             100e-3; ...
                             200e-3];

    Seq.DigitalIO.SetValue = [ 0; ...               % set all DigitalIOs to 0
                               1; ...               % set only Out 1 to 1
                               8; ...               % set only Out 4 to 1
                               2 + 32; ...          % set Out 2 to 1 and Out 6 to 1
                               63; ...              % set all DigitalIOs to 1
                               0];                  % set all DigitalIOs to 0


    % DigitalIOs
    % The digital outputs correspond to the following bits in Seq.DigitalIO.SetValue
    % Out 1: bit 0 (high: 2^0 = 1)
    % Out 2: bit 1 (high: 2^1 = 2)
    % Out 3: bit 2 (high: 2^2 = 4)
    % Out 4: bit 3 (high: 2^3 = 8)
    % Out 5: bit 4 (high: 2^4 = 16)
    % Out 6: bit 5 (high: 2^5 = 32)


    % Start measurement
    [ Raw, SeqOut, data, data_1D ] = set_sequence(HW, Seq, AQ, TX, Grad);



%% -----------------------------------------------------------------------------
% (C) Copyright 2017-2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
