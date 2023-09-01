classdef PaperSegmentation < handle
    %#ok<*AGROW>
    % PAPERSEGMENTATION Class to separate a image into lines and
    % analyze the extracted profiles.
    % author: Felix Binder
    % e-mail: felix.binder@fmt.fau.de
    % organisation: Chair of Manufacturing Metrology,
    %               Friedrich-Alexander University,
    %               91052 Erlangen, Germany
    % date (dd.mm.yyyy): 01.09.2023
    % version: 1.0
    % description: This class provides a basic separation of an uint16
    % image into line profiles and marks peaks and dips in order to
    % calculate the modulation depth of each peak-pair. Additionaly the
    % correlation of the found peak pairs with the modulation depth and a
    % polinomial fit function is provided.
    % dependency: LineProfile.m

    %% Public properties
    properties (Access = public)
        sourceImage = 0;
        voxelSize = 0;
        expectedPeaks = 0;
        paperThickness = 0;
        profiles = bnd.LineProfile(0);
        corrCoeff = 0;
        pValues = 0;
        nrOfProfiles = 0;
        paperThreshold = 0;
    end

    %% Private properties
    properties (Access = public)
        allPeakDistances = 0;
        allDIPs = 0;
    end

    %% Public methods
    methods (Access = public)
        % Contructor
        function obj = PaperSegmentation(sourceImage, voxelSize, expectedPeaks, paperThickness, paperThreshold)
            obj.sourceImage = sourceImage;
            obj.voxelSize = voxelSize;
            obj.expectedPeaks = expectedPeaks;
            obj.paperThickness = paperThickness;
            obj.paperThreshold = paperThreshold;
        end

        % Extract data for one image
        function AnalyzeSingleImageSegements(obj)
            obj.CreateLineProfiles();
            obj.ExtractProfilePeaks();
            obj.ExctractProfileValleys();
            obj.CalcPeakToPeakDistances();
            obj.CalcDIPs();
            obj.CalcCorrelationCoeff();
        end

        % Plot found peaks for debugging
        function varargout = PlotFoundValleys(obj, profileNumber)
            % Figure properties
            screenSize      = get(0,'ScreenSize');
            widthPixel      = screenSize(3);
            heightPixel     = screenSize(4);
            figureWidth     = 1280; % [px]
            figureHeight    = 480;  % [px]

            figure('Position',[widthPixel/2 - figureWidth/2, heightPixel/2 - figureHeight/2, figureWidth, figureHeight]);
            hold on;
            lineProfileX = 1:1:length(obj.profiles(profileNumber).lineProfile);
            plot(lineProfileX, obj.profiles(profileNumber).lineProfile, 'black', 'LineWidth', 2);
            plot(obj.profiles(profileNumber).lineValleyIdx, obj.profiles(profileNumber).lineValleys, 'oblue', 'MarkerSize', 10, 'LineWidth', 1.5);
            text(obj.profiles(profileNumber).lineValleyIdx+1, obj.profiles(profileNumber).lineValleys, num2str((1:numel(obj.profiles(profileNumber).lineValleys))'));
            hold off;
            xlabel('Line profile index');
            ylabel('Greyvalue in uint16');
            title(['Line profile (Found valleys: ', num2str(obj.profiles(profileNumber).foundPeaks-1), ', Profile index: ', num2str(profileNumber),')']);
            legend('Line profile', 'Found valleys');
            set(gca,'FontName','Arial');
            set(gca(),'FontSize',14);
            figureHandle = gcf();
            varargout = {figureHandle.Number};
        end

        % Plot found peaks for debugging
        function varargout = PlotFoundPeaks(obj, profileNumber)
            % Figure properties
            screenSize      = get(0,'ScreenSize');
            widthPixel      = screenSize(3);
            heightPixel     = screenSize(4);
            figureWidth     = 1280; % [px]
            figureHeight    = 480;  % [px]

            figure('Position',[widthPixel/2 - figureWidth/2, heightPixel/2 - figureHeight/2, figureWidth, figureHeight]);
            hold on;
            lineProfileX = 1:1:length(obj.profiles(profileNumber).lineProfile);
            plot(lineProfileX, obj.profiles(profileNumber).lineProfile, 'black', 'LineWidth', 2);
            plot(obj.profiles(profileNumber).linePeakIdx, obj.profiles(profileNumber).linePeaks, 'xred', 'MarkerSize', 10, 'LineWidth', 1.5);
            text(obj.profiles(profileNumber).linePeakIdx+1, obj.profiles(profileNumber).linePeaks, num2str((1:numel(obj.profiles(profileNumber).linePeaks))'));
            hold off;
            xlabel('Line profile index');
            ylabel('Greyvalue in uint16');
            title(['Line profile (Found peaks: ', num2str(obj.profiles(profileNumber).foundPeaks), ', Profile index: ', num2str(profileNumber),')']);
            legend('Line profile', 'Found peaks');
            set(gca,'FontName','Arial');
            set(gca(),'FontSize',14);
            figureHandle = gcf();
            varargout = {figureHandle.Number};
        end

        % Plot DIP and P2P distance with a polyfit
        function varargout = PlotDIPFit(obj, allPnPs, allDIPs, figureWidth, figureHeight, fontSize)
            try
                % Figure properties
                screenSize      = get(0,'ScreenSize');
                widthPixel      = screenSize(3);
                heightPixel     = screenSize(4);

                xValue = allPnPs';
                yValue = allDIPs';

                warning('off','curvefit:fit:noStartPoint')
                fitFunction = fittype('a*x^2+b*x+c');
                [fittedCurve, goodness, ~] = fit(xValue, yValue, fitFunction, 'MaxIter', 2000, 'MaxFunEvals', 5000, 'Lower', [-Inf, -Inf, -Inf], 'Upper', [Inf, Inf, Inf], 'Start', [0, 0, 0]);
                aCoefficient = fittedCurve.a;
                bCoefficient = fittedCurve.b;
                cCoefficient = fittedCurve.c;
                disp('Fit function: a*x^2 + b*x + c');
                disp(['a: ', num2str(aCoefficient, '%.6f'), ', b: ', num2str(bCoefficient, '%.6f'), ', c: ', num2str(cCoefficient, '%.6f')]);
                disp(['R-squared: ', num2str(goodness.rsquare, '%.6f'), ' RMSE: ', num2str(goodness.rmse, '%.6f')]);

                figure('Position',[widthPixel/2 - figureWidth/2, heightPixel/2 - figureHeight/2, figureWidth, figureHeight]);
                hold on;
                plot(xValue, yValue, 'xblack', 'MarkerSize', 8, 'LineWidth', 1.1);
                fittedPlot = plot(fittedCurve, '--red');
                hold off;
                set(fittedPlot, 'LineWidth', 2);
                xlabel('Page to page distance in µm');
                ylabel('Modulation depth in pct');
                title(['Line profile sets (artificial voxel size: ', num2str(obj.voxelSize, '%.2f'), ' µm)']);
                legend('Line profile peak pairs', 'Polynomial fit (2nd order)', 'Location', 'southeast');
                grid on;
                xlim([100, 800]);
                xticks(100:100:800);
                ylim([0, 100]);
                yticks(0:20:100);
                set(gca,'FontName','Arial');
                set(gca,'FontSize', fontSize);
                disp(['Used ', num2str(length(xValue)), ' datapoints.']);

                figureHandle = gcf();
                varargout = {figureHandle.Number};
            catch
                disp('Fit function cannot be displayed.');
            end
        end
    end

    %% Private methods
    methods (Access = private)

        % Fct: Extracts the line profiles in a given rectangular area
        function CreateLineProfiles(obj)
            % Add spline to increase resolution
            resolutionIncrease = 3;
            obj.voxelSize = obj.voxelSize/resolutionIncrease;
            [iSizeX, iSizeZ] = size(obj.sourceImage);
            zDetailed = 1:1/resolutionIncrease:iSizeZ;
            for idx = 1:1:iSizeX
                lineProfile = spline(1:iSizeZ, obj.sourceImage(iSizeX, :), zDetailed);
                lineProfile = lineProfile - obj.paperThreshold;
                lineProfile(lineProfile < 0) = 0;
                obj.profiles(idx,1) = bnd.LineProfile(lineProfile);
            end
        end

        % Fct: Extracts profile peaks of each line profile
        function ExtractProfilePeaks(obj)
            nrOfLineProfiles = length(obj.profiles);
            for idx = 1:1:nrOfLineProfiles
                % Get profile data and mean
                lineProfile = obj.profiles(idx).lineProfile;

                % Set elements below the mean to zero
                profileMean = mean(lineProfile);
                elementIdx = find(lineProfile > profileMean);
                filteredLineProfile = zeros(length(lineProfile),1);
                filteredLineProfile(elementIdx) = lineProfile(elementIdx);

                % Determine local peaks in the filtered profile
                profileZ = 1:1:length(lineProfile);
                [peaks, peakIdx, ~, ~] = findpeaks(filteredLineProfile, profileZ, "MinPeakDistance", floor(obj.paperThickness/obj.voxelSize));
                [~, sortedPeakIdx] = sort(peaks, 'descend');

                peaks = peaks(sortedPeakIdx(1:obj.expectedPeaks));
                peakIdx = peakIdx(sortedPeakIdx(1:obj.expectedPeaks));
                [peakIdx, resortedPeakIdx] = sort(peakIdx, 'ascend');
                peaks = peaks(resortedPeakIdx);
                peakIdx = peakIdx';

                % Find local peaks which are closer than a page
                % distance and set them to zero (indicator flag).
                maxPageIdxDistance = floor(obj.paperThickness/obj.voxelSize); %[um]
                endIdx = length(peakIdx)-1;
                for kdx = 1:1:endIdx
                    firstVal = peakIdx(kdx);
                    secVal = peakIdx(kdx+1);
                    valDistance = secVal-firstVal;
                    if valDistance <= maxPageIdxDistance
                        if peaks(kdx) >= peaks(kdx+1)
                            peaks(kdx+1) = 0;
                            peakIdx(kdx+1) = 0;
                        else
                            peaks(kdx) = 0;
                            peakIdx(kdx) = 0;
                        end
                    end
                end

                % Remove data with zero and assign peak values
                peakIdx = peakIdx(peakIdx~=0);
                obj.profiles(idx).linePeakIdx = peakIdx;
                peaks = peaks(peaks~=0);
                obj.profiles(idx).linePeaks = peaks;
                obj.profiles(idx).foundPeaks = length(peaks);
            end
        end

        % Fct: Extracts profile peaks of each line profile
        function ExctractProfileValleys(obj)
            nrOfLineProfiles = length(obj.profiles);
            for idx = 1:1:nrOfLineProfiles
                nrOfFoundPeaks = obj.profiles(idx).foundPeaks;
                for kdx = 1:1:nrOfFoundPeaks-1
                    startIdx = obj.profiles(idx).linePeakIdx(kdx);
                    endIdx = obj.profiles(idx).linePeakIdx(kdx+1);
                    peakToPeakData = obj.profiles(idx).lineProfile(startIdx:endIdx);
                    [minData, minIdx] = min(peakToPeakData);
                    obj.profiles(idx).lineValleys(kdx) = minData;
                    obj.profiles(idx).lineValleyIdx(kdx) = minIdx+startIdx-1;
                end
            end
        end

        % Fct: Calculates the peak to peak distance in [um] of each profile
        function CalcPeakToPeakDistances(obj)
            nrOfLineProfiles = length(obj.profiles);
            for idx = 1:1:nrOfLineProfiles
                peakToPeakDistance = diff(obj.profiles(idx).linePeakIdx)*obj.voxelSize;
                obj.profiles(idx).peakDistances = peakToPeakDistance';
            end
        end

        % Fct: Calculates the modulation depth (DIP, ASTM 2597E)
        function CalcDIPs(obj)
            % DIP = 100*(A+B-2C)/(A+B)
            nrOfLineProfiles = length(obj.profiles);
            for idx = 1:1:nrOfLineProfiles
                nrOfFoundPeaks = obj.profiles(idx).foundPeaks;
                for kdx = 1:1:nrOfFoundPeaks-1
                    firstPeak = obj.profiles(idx).linePeaks(kdx);
                    secPeak = obj.profiles(idx).linePeaks(kdx+1);
                    valley = obj.profiles(idx).lineValleys(kdx);
                    modulationDepth = 100*(firstPeak+secPeak-2*valley)/(firstPeak+secPeak);
                    obj.profiles(idx).peakDIPs(kdx) = modulationDepth;
                end
            end
        end

        % Fct: Calculates correlation coefficients (mismatches excluded)
        function CalcCorrelationCoeff(obj)
            dispLinearCorrelationCoefficient = false; % debugging
            try
                cumIdx = 1;
                nrOfLineProfiles = length(obj.profiles);
                for idx = 1:1:nrOfLineProfiles
                    nrOfFoundPeaks = obj.profiles(idx).foundPeaks;

                    for kdx = 1:1:nrOfFoundPeaks-1
                        obj.allPeakDistances(cumIdx) = obj.profiles(idx).peakDistances(kdx);
                        obj.allDIPs(cumIdx) = obj.profiles(idx).peakDIPs(kdx);
                        cumIdx = cumIdx +1;
                    end

                end
                [rCoeff,pValue] = corrcoef(obj.allPeakDistances, obj.allDIPs);
                obj.corrCoeff = rCoeff;
                obj.pValues = pValue;
                obj.nrOfProfiles = nrOfLineProfiles;
                if dispLinearCorrelationCoefficient == true
                    disp(['Correlation: R = ', num2str(rCoeff(2), '%.2f'), ', p = ', num2str(pValue(2), '%.2f')]);
                end
            catch
                disp('Correlation coefficient cannot be determined. (Class: PaperSegmentation, Function: CalcCorrelationCoeff');
            end
        end
    end
end