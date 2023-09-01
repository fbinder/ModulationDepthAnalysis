% MODULATIONDEPTHANALYSIS script to analyse a reconstructed paper stack
% volume.
% author: Felix Binder
% e-mail: felix.binder@fmt.fau.de
% organisation: Chair of Manufacturing Metrology,
%               Friedrich-Alexander University,
%               91052 Erlangen, Germany
% date (dd.mm.yyyy): 01.09.2023
% version: 1.0
% description: This is the accompanied script for the modulation depth
% analysis in the submitted publication: "Preserving Fragile History:
% Assessing the Feasibility of Segmenting Digitized Historical Documents
% with Modulation Depth Analysis" to Heritage. Please note that the code
% needs a reconstructed uint16 volume provided to the variable "rawVolume".
% Measured data from the publication can be requested via e-mail.
% dependency: LineProfile.m, FctHelper.m, PaperSegmentation.m

%% PREPERATION
clearvars;
close all;
clc;
tic;

%% USER INPUT
% File import and variable definition
fileName = 'zpp_R254_Acrylpresse_vMuster3_v2_170kV250uA2s8x2BM025Cu2050Proj280X 2023-1-27 20-18-29.vgi';
filePath = 'C:\Users\bnd.FMT\Desktop\PaperSegmentation_v10\Measurement\Wdh\';
expectedPeaks   = 38;     % [-] - Number of expected peaks
paperThickness  = 130;    % [um] - Minimal paper thickness
p2pMin          = paperThickness;    % [um] - Minimal peak-to-peak distance
p2pMax          = 396;               % [um] - Constructed max peak-to-peak distance
voxelSize       = 42.1205204786074;  % [um] - Reconstructed voxel size
angularIterations = 360;  % [-] - Angular iterations steps in between 0°-360°

% PLEASE ASSIGN A RECONSTRUCTED UINT16 VOLUME TO THE VARIABLE RAWVOLUME
% If you need to switch dimension, use "permute()".
rawVolume = zeros(2048,2048,2048, 'uint16');

%% Main
minion = bnd.FctHelper();

% Rescale volume from manual area (manualMIN - manualMAX [uint16]) to max uint16
rescaledVolume = minion.RescaleVolumeUINT16(rawVolume, 0, 65535);
clear rawVolume recoVolume;

% Equalize the x and y dimensions (added zeros to the smaller dimension)
equalVolume = minion.CreateXYEqualSizeVolume(rescaledVolume);
clear rescaledVolume;

% Rotate volume (if necessary)
equalVolume = imrotate3(equalVolume, 90, [0, 1, 0], 'nearest', 'crop', 'FillValues', 0); % Rotate along z

% Display centre slice and define cylindrical ROI
[vSizeX, vSizeY, vSizeZ] = size(equalVolume);
imageSlice = minion.GetSlice(equalVolume, floor(vSizeZ/2), 'xy');
figureNr = minion.ViewSlice(imageSlice, floor(vSizeZ/2), 'xy', 1024, 768, 16);
figureHandle = figure(figureNr);
circleROI = minion.GetCurrentFigureROI('circle', 'Please select the inner paper stack ROI and press any key to continue.');
circleCenter = floor(circleROI.Center);
circleRadius = floor(circleROI.Radius);
close(figureHandle);

% Translate ROI centre to the image centre for rotation
circleCenter = fliplr(circleCenter); % because image() has order (y,x,z) instead of (x,y,z);
equalVolume = imtranslate(equalVolume, [floor(vSizeX/2)-circleCenter(1), floor(vSizeY/2)-circleCenter(2), 0], 'linear');

% Define a slice ROI to avoid the acrylic fixture
imageSlice = minion.GetSlice(equalVolume, floor(vSizeY/2), 'xz');
figureNr = minion.ViewSlice(imageSlice, floor(vSizeY/2), 'xz', 1024, 768, 16);
figureHandle = figure(figureNr);
partialSliceROI = minion.GetCurrentFigureROI('rectangle', 'Please select the paper stack area without the fixture. Press any key to continue.');
partialSliceROIPosition = partialSliceROI.Position; % [xmin, ymin, width, height]
close(figureHandle);

partialSliceZStart = floor(partialSliceROIPosition(2));
partialSliceZEnd = floor(partialSliceROIPosition(2) + partialSliceROIPosition(4));
partialSliceXStart = floor(partialSliceROIPosition(1));
partialSliceXEnd = floor(partialSliceROIPosition(1) + partialSliceROIPosition(3));

% Get paper reference ROI
imageSlice = minion.GetSlice(equalVolume, floor(vSizeY/2), 'xz');
figureNr = minion.ViewSlice(imageSlice, floor(vSizeY/2), 'xz', 1024, 768, 16);
figureHandle = figure(figureNr);
referenceROI = minion.GetCurrentFigureROI('rectangle', 'Please select the paper reference ROI and press any key to continue.');
referenceROIPosition = referenceROI.Position; % [xmin, ymin, width, height]
close(figureHandle);

referenceROIPosition = [referenceROIPosition(2), referenceROIPosition(1), referenceROIPosition(4), referenceROIPosition(3)]; % rearrange to keep x,y order
referencePaperROI = imcrop(imageSlice, referenceROIPosition);
[meanReferencePaperROI, ~] = minion.CalcMeanConfidenceImage(referencePaperROI, 0.95);

[iSizeX, iSizeZ] = size(imageSlice);
lineProfiles = zeros(iSizeX-2, iSizeZ-2, angularIterations);
allP2P = [];
allDIPs = [];

rotationAngles = linspace(0, 360, angularIterations+1);
rotationAngles = rotationAngles(1:end-1); % Drop double angles

% Interate angular for line segmentation XZ plane: x(x) and y(z)
for idx = 1:1:angularIterations

    rotatedVolume = imrotate3(equalVolume, rotationAngles(idx), [0, 0, 1], 'nearest', 'crop', 'FillValues', 0); % Rotate along z

    midSliceIdx = floor(vSizeY/2);
    frontSlice = minion.GetSlice(rotatedVolume, midSliceIdx-1, 'xz');
    midSlice = minion.GetSlice(rotatedVolume, midSliceIdx, 'xz');
    backSlice = minion.GetSlice(rotatedVolume, midSliceIdx+1, 'xz');

    % Calculate the 3x3 mean for every coordinate
    [iSizeX, iSizeZ] = size(midSlice);
    partialVolume = zeros(3,3,3);
    pxOffset = 1;

    % kdx is the running index along x
    for kdx = partialSliceXStart+pxOffset:1:partialSliceXEnd-pxOffset
        % pdx is the running index along z
        for pdx = partialSliceZStart+pxOffset:1:partialSliceZEnd-pxOffset
            % Assign partial volume
            partialVolume(:,:,3) = frontSlice(kdx-1:kdx+1, pdx-1:pdx+1);
            partialVolume(:,:,2) = midSlice(kdx-1:kdx+1, pdx-1:pdx+1);
            partialVolume(:,:,1) = backSlice(kdx-1:kdx+1, pdx-1:pdx+1);

            % Calculate mean and add to lineprofiles
            partialMean = mean(partialVolume, 'all');
            lineProfiles(kdx-pxOffset, pdx-pxOffset, idx) = partialMean;
        end
    end

    % Reduce matrix data to nonzero elements
    profilesROI = lineProfiles(partialSliceXStart:partialSliceXEnd-2*pxOffset, :, :);

    paper = bnd.PaperSegmentation(profilesROI(:,:,idx), voxelSize, expectedPeaks, paperThickness, meanReferencePaperROI);
    minion.ExtDisp(['Rotation angle: ', num2str(rotationAngles(idx)), '°']);
    paper.AnalyzeSingleImageSegements();

    allP2P = [allP2P, paper.allPeakDistances]; %#ok<AGROW>
    allDIPs = [allDIPs, paper.allDIPs]; %#ok<AGROW>
end

histStart = 140-7;  % in um
histEnd = 378+7;    % in um
binWidth = 14;      % in um
histBins = floor((histEnd-histStart)/binWidth);
if mod(histEnd-histStart, binWidth) ~= 0
    disp('Be careful. The bin widths are uneven distributed.');
end

% Apply a custom histogram and calculate the mean and the std for each bin.
sizeBeforeReduction = length(allP2P);
validIdx = find(allP2P >= histStart & allP2P <= histEnd);
wdhP2P = allP2P(validIdx);
wdhDIPs = allDIPs(validIdx);

usedPercentage = length(validIdx)/sizeBeforeReduction*100;
clear allDIPs allP2P;

% Calculate mean and std
meanDIP = zeros(histBins,1);
stdDIP = zeros(histBins,1);
binP2P = histStart+binWidth/2:binWidth:histEnd-binWidth/2;
binP2P = binP2P';

% Iterate over every bin available
for idx = 1:1:histBins
    lowBinBorder = binP2P(idx)-binWidth/2;
    uppBinBorder = binP2P(idx)+binWidth/2;

    validIdx = find(wdhP2P > lowBinBorder & wdhP2P <= uppBinBorder);
    binDIP = wdhDIPs(validIdx); %#ok<FNDSB> 

    meanDIP(idx) = mean(binDIP);
    stdDIP(idx) = std(binDIP);
end

% Remove possible NaN values
validIdx = find(~isnan(meanDIP));
meanDIP = meanDIP(validIdx);
binP2P = binP2P(validIdx);
stdDIP = stdDIP(validIdx);

warning('off','curvefit:fit:noStartPoint')
fitFunction = fittype('a*x^2+b*x+c');
[fittedCurve, goodness, ~] = fit(binP2P, meanDIP, fitFunction, 'MaxIter', 2000, 'MaxFunEvals', 5000, 'Lower', [-Inf, -Inf, -Inf], 'Upper', [Inf, Inf, Inf], 'Start', [0, 0, 0]);
aCoefficient = fittedCurve.a;
bCoefficient = fittedCurve.b;
cCoefficient = fittedCurve.c;

% left figure
figure();
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');

nexttile;
hold on;
grid on;
plot(wdhP2P, wdhDIPs, 'xblack', 'MarkerSize', 10, 'LineWidth', 1.5);
xline(130, '--black', 'LineWidth', 1.5);
xline(392, '--black', 'HandleVisibility', 'off', 'LineWidth', 1.5);
hold off;

xlim([120, 400]);
xticks(100:40:800);
ylim([0,100]);
yticks(0:20:100);
xlabel('Peak-to-peak distance in µm');
ylabel('Modulation depth in %');

% right figure
nexttile;
hold on;
grid on;
errorbar(binP2P, meanDIP, stdDIP, 's', 'MarkerSize', 10, 'MarkerEdgeColor', [0/255, 0/255, 0/255], 'MarkerFaceColor', [0/255, 0/255, 0/255], 'LineWidth', 1.5, 'Color', [0/255, 0/255, 0/255]);
fitY = binP2P.^2.*aCoefficient + binP2P.*bCoefficient + cCoefficient;
plot(binP2P, fitY, '-.red', 'LineWidth', 2);
hold off;

xline(130, '--black', 'LineWidth', 1.5);
xline(392, '--black', 'HandleVisibility', 'off', 'LineWidth', 1.5);
xlim([120, 400]);
xticks(100:40:800);
ylim([0,100]);
yticks(0:20:100);
xlabel('Peak-to-peak distance in µm');
ylabel('Modulation depth in %');

% Resize figure
cf = gcf();
set(findall(gcf,'-property','FontSize'), 'FontSize', 20);
set(findall(gcf,'-property','FontName'), 'FontName', 'Cambria');
minion.ResizeTiledLayout(1400, 480, 18);

nexttile(2);
legend('Binned mean data (deviation: \pm 1\sigma)', ['Polynomial fit (R^2: ', num2str(goodness.rsquare, '%.2f'), ')'], 'Expected distances', 'Location', 'northwest', 'FontSize', 16);

nexttile(1);
legend('Raw data (unbinned)', 'Expected distances', 'Location', 'northwest', 'FontSize', 16);
toc;