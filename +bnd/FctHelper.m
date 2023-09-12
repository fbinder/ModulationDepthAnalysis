classdef FctHelper < handle
    % FCTHELPER is a helper class handling uint16 data.
    % author: Felix Binder
    % e-mail: felix.binder@fmt.fau.de
    % organisation: Chair of Manufacturing Metrology,
    %               Friedrich-Alexander University,
    %               91052 Erlangen, Germany
    % date (dd.mm.yyyy): 01.09.2023
    % version: 1.0
    % description: This class is a set of useful functions handling uint16 
    % data.

    % Public properties
    properties
        colors;
        brightGreyColormap;
    end

    % Public methods
    methods
        % Constructor
        function obj = FctHelper()
            % Load colormaps
            fileLocation = fileparts(mfilename('fullpath'));
            addpath(obj.FullFilePath('libs\', fileLocation));
            obj.brightGreyColormap = load('brightGreyColormap.mat');
            obj.brightGreyColormap = obj.brightGreyColormap.BrightGreyColormap;

            % Define colors
            obj.colors.cBlack = [0.0000, 0.0000, 0.0000];
            obj.colors.cWhite = [1.0000, 1.0000, 1.0000];
            obj.colors.cYellow = [1.0000, 1.0000, 0.0000];
            obj.colors.cfauBlue = [0.0706, 0.1922, 0.3961];
            obj.colors.cfauGrey = [0.6000, 0.6000, 0.6000];
            obj.colors.ctmLightBlue = [0.5333, 0.8000, 0.9333];
            obj.colors.ctmBlueGreenish = [0.2667, 0.6667, 0.6000];
            obj.colors.ctmDarkGreen = [0.0667, 0.4667, 0.2000];
            obj.colors.ctmDarkBlue = [0.2000, 0.1333, 0.5333];
            obj.colors.ctmOcher = [0.8667, 0.8000, 0.4667];
            obj.colors.ctmLime = [0.6000, 0.6000, 0.2000];
            obj.colors.ctmLightPink = [0.8000, 0.4000, 0.4667];
            obj.colors.ctmDeepPurple = [0.5333, 0.1333, 0.3333];
            obj.colors.ctmViolet = [0.6667, 0.2667, 0.6000];
            obj.colors.ctmLightGrey = [0.8667, 0.8667, 0.8667];
        end

        % Fct: Normalizes a uint16 image to max contrast (0-65535)
        function normUint16Image = NormalizeUint16Img(~, uint16Image)
            offsetLow = double(min(min(uint16Image)));
            normUint16Image = double(uint16Image)-offsetLow;
            stretchFactor = 65535./max(max(normUint16Image));
            normUint16Image = normUint16Image.*stretchFactor;
            normUint16Image = uint16(normUint16Image);
        end

        % Fct: Normalizes a uint16 image to max contrast (0-65535)
        function rescaledVolume = RescaleVolumeUINT16(~, volume, newMin, newMax)
            volume = double(volume);
            rescaledVolume = volume-newMin;
            rescaledVolume(rescaledVolume < 0) = 0;
            stretchFactor = 65535/(newMax-newMin);
            rescaledVolume = rescaledVolume.*stretchFactor;
            rescaledVolume(rescaledVolume > 65535) = 65535;
            rescaledVolume = uint16(rescaledVolume);
        end

        % Fct: Get current screen size
        function [width, height] = GetScreenSize(~)
            screenSize = get(0,'ScreenSize');
            width = screenSize(3);
            height = screenSize(4);
        end

        % Fct: Combine full file path
        function fullFilePath = FullFilePath(~, filename, filepath)
            if(filepath(end) ~= '\')
                filepath = [filepath, '\'];
            end
            fullFilePath = [filepath, filename];
        end

        % Fct: Resize tiled layout
        function ResizeTiledLayout(obj, widthPx, heightPx, fontSize)
            [screenWidth, screenHeight] = obj.GetScreenSize();
            cf = gcf();
            set(findall(gcf,'-property','FontSize'), 'FontSize', fontSize);
            cf.Position = [screenWidth/2 - widthPx/2, screenHeight/2 - heightPx/2, widthPx, heightPx];
        end

        % Fct: Display volume with histogram
        function varargout = VolumeHistogramUINT16(obj, volume, widthPx, heightPx, visibilityDivider)
            [screenWidth, screenHeight] = obj.GetScreenSize();

            % Backgroundcolor not availbe in 2022b+
            %figHandle1 = figure('Position', [screenWidth/2 - widthPx/2, screenHeight/2 - heightPx/2, widthPx/2, heightPx]);
            %volHandle = volshow(volume, 'BackgroundColor', obj.colors.ctmLightGrey);
            volHandle = volshow(volume);
            volHandle.Colormap = flipud(obj.brightGreyColormap);
            customAlphaMap = ones(256,1);
            customAlphaMap(1:floor(256*visibilityDivider)) = 0;
            volHandle.Alphamap = customAlphaMap;

            figHandle2 = figure('Position', [screenWidth/2, screenHeight/2 - heightPx/2, widthPx/2, heightPx]);
            sizeRaw = size(volume);
            nrBins = floor(nthroot(sizeRaw(1)*sizeRaw(1)*sizeRaw(1), 3));
            minVol = min(min(min(volume)));
            maxVol = max(max(max(volume)));

            hold on;
            [counts, bins] = imhist(volume, nrBins);
            counts(1) = 0;
            counts(end) = 0;
            barHandle = bar(bins, counts);
            grid on;
            xline(floor(65535*visibilityDivider), 'LineWidth', 1.5, 'Color', obj.colors.cfauGrey);
            barHandle.BarWidth = 1;
            barHandle.FaceColor = obj.colors.cfauBlue;
            xlim([0,65535]);
            xticks(0:5000:65535);
            xlabel('Bins in uint16');
            ylabel('Voxel counts (no extremes)');
            title(['Volume histogram (min: ', num2str(minVol), ' , max: ', num2str(maxVol), ')']);
            maxCounts = max(counts);
            text(floor(65535*(visibilityDivider*1.05)), maxCounts*0.95, 'Visible');
            legend(['Histogram (N = ', num2str(nrBins), ')'], 'Visbibility line');
            %varargout{1} = num2cell(figHandle1.Number);
            varargout{1} = num2cell(figHandle2.Number);
        end

        % Fct: Create volume with equal sizes in x and y
        function evenVolume = CreateXYEqualSizeVolume(~, volume)
            [xSize, ySize, zSize] = size(volume);

            % Reduce dimension by -1 to avoid odd offsets
            if mod(xSize, 2) ~= 1
                xSize = xSize -1;
            end
            if mod(ySize, 2) ~= 1
                ySize = ySize -1;
            end

            % Find the bigger dimension to copy data into
            if xSize > ySize
                targetSize = xSize;
            else
                targetSize = ySize;
            end

            % 0 = no signal (white), 65535 max signal (dark grey)
            evenVolume = zeros(targetSize, targetSize, zSize, 'uint16');
            xSizeOffset = floor((targetSize-xSize)/2);
            ySizeOffset = floor((targetSize-ySize)/2);

            evenVolume(1+xSizeOffset:targetSize-xSizeOffset, 1+ySizeOffset:targetSize-ySizeOffset, :) = volume(1:xSize, 1:ySize, :);
        end

        % Fct: Get volume slices with order
        function slice = GetSlice(~, volume, fixed, direction)
            % According to order of [x,y,z]
            switch direction
                case 'xy'
                    slice = squeeze(volume(:,:,fixed));
                case 'yx'
                    slice = squeeze(volume(:,:,fixed))';
                case 'xz'
                    slice = squeeze(volume(:,fixed, :));
                case 'zx'
                    slice = squeeze(volume(:,fixed,:))';
                case 'yz'
                    slice = squeeze(volume(fixed,:,:));
                case 'zy'
                    slice = squeeze(volume(fixed,:,:))';
            end
        end

        % Fct: Calculate mean and confidence interval of image
        function [meanImg, confImg] = CalcMeanConfidenceImage(~, img, level)
            [imgSizeX, imgSizeY] = size(img);
            degreesOfFreedom = imgSizeX*imgSizeY-1;

            % Two sided t-statistics --> alpha/2
            alphaTwoSided = 1-(1-level)/2;
            tCritical = tinv(alphaTwoSided, degreesOfFreedom);

            % Confidence intervals
            meanImg = mean2(img);
            stdSample = std2(img);
            confImg = stdSample*tCritical;
        end

        % Fct: Display a volume slice
        function varargout = ViewSlice(obj, slice, fixed, direction, widthPx, heightPx, fontSize)
            figHandle = figure();
            tiledlayout(1,1,'Padding', 'none', 'TileSpacing', 'compact');
            nexttile;
            image(slice'); % Transform to keep [x,y,z] order
            colormap(flipud(gray(65535)));
            sliceSize = size(slice);
            maxSize = max(sliceSize);
            pbaspect([sliceSize(1)/maxSize, sliceSize(2)/maxSize, 1]);
            colorbar();
            xlabel([direction(1), '-axis in px']);
            ylabel([direction(2), '-axis in px']);
            title(['Slice: ', num2str(fixed), ' (', direction, '-plane)']);
            obj.ResizeTiledLayout(widthPx, heightPx, fontSize);
            varargout = num2cell(figHandle.Number);
        end

        % Fct: Get manual ROI from the current open figure
        function regionOfInterest = GetCurrentFigureROI(obj, type, msg)
            if ~exist('msg', 'var')
                dispMSG = 'Please define the ROI and PRESS any key to continue!';
            else
                dispMSG = msg;
            end

            switch type
                case 'rectangle'
                    obj.ExtDisp(dispMSG);
                    regionOfInterest = drawrectangle("Color", [1, 1, 0], "FaceAlpha", 0);
                    pause();
                case 'circle'
                    obj.ExtDisp(dispMSG);
                    regionOfInterest = drawcircle("Color", [1, 1, 0], "FaceAlpha", 0);
                    pause();
                otherwise
                    assert(false, 'Currently supported types: rectangle, circle');
            end
        end

        % Fct: Get a timestamp for display messages
        function ExtDisp(~, msg)
            utc_time = datetime('now', 'timezone', 'utc');
            utc_time.Format = 'yyyy-MM-dd HH:mm:ss';
            disp(['>>', convertStringsToChars(string(utc_time)), ' (UTC+0): ', msg]);
        end
    end
end

