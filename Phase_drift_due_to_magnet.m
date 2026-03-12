%% Phase Drift Mapping using repeated FLASH images (2D)
LoadSystem;
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0, [], 1);

% --- Imaging parameters (edit to match your stable setup) ---
Seq.AQSlice(1).nRead       = 64;
Seq.AQSlice(1).nPhase(2)   = 64;
Seq.AQSlice(1).sizeRead    = 0.010;
Seq.AQSlice(1).sizePhase(2)= 0.010;
Seq.AQSlice(1).thickness   = 0.005;

Seq.tEcho = 5e-3;   % TE (s)  <-- set to your TE
Seq.tRep  = 50e-3;  % TR (s)  <-- set to your TR (must be > TE)
Seq.average = 1;

% Loop settings
Nloops     = 120;       % e.g. 120 frames
LoopBreak  = 2.0;       % seconds between frames
refLoop    = 1;         % reference frame index

% Storage
StartTime  = nan(Nloops,1);
PhiStack   = [];        % will become (Nx, Ny, Nloops)
MagStack   = [];

hf = figure; colormap gray;

for k = 1:Nloops
    % schedule time like Pure Devices style (optional)
    if k == 1
        Seq.StartSequenceTime = now*24*3600 + LoopBreak + 1;
    else
        Seq.StartSequenceTime = SeqOut.StartSequenceTime + LoopBreak;
        HW.tRepInit = 0.05; % reduce re-init overhead for short breaks
    end

    % --- Run FLASH acquisition ---
    % Choose ONE that you already know works:
    % [SeqOut, mySave] = Example_Flash(HW, Seq, AQ, TX, Grad, mySave);
    % OR:
    [SeqOut, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

    % Record time
    StartTime(k) = SeqOut.StartSequenceTime;

    % --- Extract complex image ---
    % You may need to adjust this line depending on the struct fields.
    % Common possibilities:
    % img = SeqOut.data.Image;      % might be complex
    % img = SeqOut.data.ImageZ;     % alternative reconstructed image
    img = SeqOut.data.Image;  % <-- edit if needed

    % Build stacks
    if k == 1
        [Nx, Ny] = size(img(:,:,1)); % adjust if image has singleton dims
        PhiStack = nan(Nx,Ny,Nloops);
        MagStack = nan(Nx,Ny,Nloops);
    end

    img2 = squeeze(img(:,:,1));  % ensure 2D
    PhiStack(:,:,k) = angle(img2);
    MagStack(:,:,k) = abs(img2);

    % --- Quick live plot of Δphi vs first frame ---
    if k == refLoop
        phiRef = PhiStack(:,:,k);
    end
    dphi = angle(exp(1i*(PhiStack(:,:,k) - phiRef))); % wrapped difference

    subplot(1,2,1);
    imagesc(MagStack(:,:,k)); axis image off;
    title(sprintf('Magnitude (loop %d)', k));

    subplot(1,2,2);
    imagesc(dphi); axis image off;
    title('\Delta\phi vs reference (wrappedad)');
    drawnow;
end

% Convert time axis to seconds from start
t = StartTime - StartTime(1);

%% ROI drift plots (center vs edge)
cx = round (Nx/2);
cy = round(Ny/2);
half=1;
xIdx=max(1,cx-half): min(Nx, cx+half);
yIdx= max(1,cy-half): min(Ny, cy+half);
roiCenter={xIdx, yIdx};
xE=1:min(Nx,3);
yE=1:min (Ny,3);
roiEdge ={xE,yE};

phiRef = PhiStack(:,:,refLoop);
dphiStack = angle(exp(1i*(PhiStack - phiRef)));

centerMean = squeeze(mean(dphiStack(roiCenter{1}, roiCenter{2}, :), [1 2], 'omitnan'));
edgeMean   = squeeze(mean(dphiStack(roiEdge{1}, roiEdge{2}, :), [1 2], 'omitnan'));

figure;
plot(t, centerMean, '-o'); hold on;
plot(t, edgeMean, '-o');
grid on;
xlabel('Time (s)');
ylabel('Mean \Delta\phi (rad)');
legend('Center ROI','Edge ROI');
title('Phase drift over time (relative to reference frame)');

