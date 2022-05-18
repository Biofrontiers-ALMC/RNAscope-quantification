%RNAscope quantification
clearvars
clc

%Parameters
file = 'D:\Projects\ALMC Tickets\T362-McNulty\data\1c_PLVGATD2D1.nd2';

nuclThreshold = 150;  %Greyscale to threshold nuclei in GFP channel
maxNuclSize = 50;  %Largest valid nucleus in pixels

%Threshold to count puncta by maximum brightness in the cell
ch570Threshold = 800;
ch650Threshold = 3000;

%% Begin processing

%Read in nd2 file
nd2 = BioformatsImage(file);

%Compute the maximum intensity projection of each channel
for iZ = 1:nd2.sizeZ

    currImage = zeros(nd2.height, nd2.width, nd2.sizeC, 'uint16');

    for iC = 1:nd2.sizeC
        currImage(:, :, iC) = getPlane(nd2, iZ, iC, 1);
    end

    if iZ == 1
        mip = currImage;
    else 
        mip = max(mip, currImage);
    end

end

%Make a mask using the GFP channel
nuclMask = mip(:, :, 2) > nuclThreshold;
nuclMask = imopen(nuclMask, strel('disk', 1));

dd = -bwdist(~nuclMask);
dd(~nuclMask) = -Inf;
dd = imhmin(dd, 0.5);

L = watershed(dd);

nuclMask(L == 0) = false;
nuclMask = bwareaopen(nuclMask, maxNuclSize);

%Grow mask
%nuclMask = imdilate(nuclMask, [1 1; 1 1]);

%Display an image to valide the mask
figure(1)
showoverlay(mip(:, :, 2), bwperim(nuclMask))
title('GFP channel')

%% Count number of cells
for iC = 1:3

    %Calculate intensity in each channel
    chCellData{iC} = regionprops(nuclMask, mip(:, :, iC + 1), 'meanIntensity', 'maxIntensity', 'centroid');

    if iC == 1
        %Use mean intensity for FITC
        validCells{iC} = true(1, numel(chCellData{iC}));

    elseif iC == 2
        validCells{iC} = [chCellData{iC}.MaxIntensity] > ch570Threshold;

    elseif iC == 3
        validCells{iC} = [chCellData{iC}.MaxIntensity] > ch650Threshold;

    end

end

%% Count number of cells

%Number of GFP labeled cells
numVGaT = nnz(validCells{1})

%Number of GFP labeled cells with puncta ONLY in 520 channel
numVGaTplusDrd2 = nnz(validCells{1} & validCells{2} & ~validCells{3})

%Number of GFP labeled cells with puncta ONLY in 650 channel
numVGaTplusDrd1 = nnz(validCells{1} & validCells{3} & ~validCells{2})

%Number of GFP labeled cells with puncta in BOTH in 520 and 650 channels
numVGaTplusDrd1plusDrd2 = nnz(validCells{1} & validCells{2} & validCells{3})

%% Plot images for sanity checks
%
% These images show a merged image of the 520 and 650 channels. There is a
% green outline for each GFP labeled cell. Cells that meet the threshold
% criteria for puncta are marked with a red circle (puncta in 520 nm
% channel) or a white cross (puncta in 650 nm channel).

%Normalize the 520 and 650 data for better visualization
ch520Norm = double(mip(:, :, 3));
ch520Norm = (ch520Norm - min(ch520Norm(:)))/(max(ch520Norm(:)) - min(ch520Norm(:)));

ch650Norm = double(mip(:, :, 4));
ch650Norm = (ch650Norm - min(ch650Norm(:)))/(max(ch650Norm(:)) - min(ch650Norm(:)));

%Create an RGB image of the 520 (yellow) and 650 (magenta) channels
alphaBlend = 0.5;

red = alphaBlend * ch520Norm + (1 - alphaBlend) * ch650Norm;
green = alphaBlend * ch520Norm;
blue = (1 - alphaBlend) * ch650Norm;

rgbOut = cat(3, red, green, blue);

%Overlay markers to identify cells with puncta
figure(2)
showoverlay(rgbOut, bwperim(nuclMask), 'Opacity', 40)
hold on
for iC = [2, 3]
    for ii = 1:numel(chCellData{iC})
        if validCells{iC}(ii)
            if iC == 2
            plot(chCellData{iC}(ii).Centroid(1), chCellData{iC}(ii).Centroid(2), 'ro')
            else
                plot(chCellData{iC}(ii).Centroid(1), chCellData{iC}(ii).Centroid(2), 'wx')
            end
        end
    end
end
hold off
