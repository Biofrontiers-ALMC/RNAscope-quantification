%Counts number of cells

clearvars
clc

file = 'D:\Projects\ALMC Tickets\T362-McNulty\data\1c_PLVGATD2D1.nd2';

nd2 = ND2reader(file);

%Compute the maximum intensity projection of each channel
mip = getImage(nd2, 1, 1);
for iZ = 2:nd2.sizeZ

    currImage = getImage(nd2, iZ, 1);
    mip = max(mip, currImage);

end

%%
%Make a mask
% nuclMask = mip(:, :, 1) > 150;
nuclMask = mip(:, :, 2) > 100;
nuclMask = imopen(nuclMask, strel('disk', 2));

dd = -bwdist(~nuclMask);
dd(~nuclMask) = -Inf;
dd = imhmin(dd, 1);

L = watershed(dd);

nuclMask(L == 0) = false;
nuclMask = bwareaopen(nuclMask, 50);

%Grow mask
%nuclMask = imdilate(nuclMask, [1 1; 1 1]);

showoverlay(mip(:, :, 2), bwperim(nuclMask))

%% Count number of cells

for iC = 1:3

    %Calculate intensity in each channel
    chCellData{iC} = regionprops(nuclMask, mip(:, :, iC + 1), 'meanIntensity', 'maxIntensity', 'centroid');

    if iC == 1
        %Use mean intensity for FITC
        validCells{iC} = [chCellData{iC}.MeanIntensity] > 100;

    elseif iC == 2
        validCells{iC} = [chCellData{iC}.MaxIntensity] > 2000;

    elseif iC == 3
        validCells{iC} = [chCellData{iC}.MaxIntensity] > 3000;

    end

end

%% Analyze resulting data

%Count number of cells
numVGaT = nnz(validCells{1});

numVGaTplusDrd2 = nnz(validCells{1} & validCells{2});
numVGaTplusDrd1 = nnz(validCells{1} & validCells{3});
numVGaTplusDrd1plusDrd2 = nnz(validCells{1} & validCells{2} & validCells{3});

%% Plotting

iC = 2;

showoverlay(mip(:, :, iC + 1), bwperim(nuclMask), 'Opacity', 40)
hold on
for ii = 1:numel(chCellData{iC})
    if validCells{iC}(ii)
        plot(chCellData{iC}(ii).Centroid(1), chCellData{iC}(ii).Centroid(2), 'ro')
    end
end
hold off
