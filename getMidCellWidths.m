function [binnedLengthSeptumWidthEstTime] = getMidCellWidths(cellList,birthLength,divisionLength)
%This function takes the output of function getExtraDataLoop (Oufti
%cellList structure with 'extra data' appended) and:
% - for cells with 1 or 2 nucleoids between the estimated length at birth 
%(birthLength) and the estimated length at division (divisionLength)
% - calculates the minimum width in the middle 3 segments ('midcell width')
% - estimated the time since birth of the cell using the formula:
%   - length at time x = (length at birth) * 2^x/tau
% - averages the mid-cell width for cells within length bins.

%output matrix (binnedLengthSeptumWidthEstTime):
%column 1: average length within each length bin
%column 2: average mid-cell width within each length bin
%column 3: estimated time since birth of each length bin
%column 4: number of measured cells within each length bin
%column 5: the standard deviation of midcell width for each length bin
%column 6: the standard error of the mean of the midcell width for each
%length bin

%Martin White 2018,

umperpixel = 0.1;           %pixel to micrometer conversion factor
tau = 17.16;                %minutes
binEdges=3:0.1:10;          %binning for cells in micrometers

% Step 1: for each cell meeting the following requirements: contains either
% 1 or 2 nucleoids ('objects'), longer than estimated average length at
% birth and shorter than estimated average length at division
cellCount=0;
for frame = 1:length(cellList.meshData)
    for cellNum = 1:length(cellList.meshData{frame})
        if isempty(cellList.meshData{frame}{cellNum}) ...
                || ~isfield(cellList.meshData{frame}{cellNum},'mesh') ...
                ||length(cellList.meshData{frame}{cellNum}.mesh)<4 ...
                || length(cellList.meshData{frame}{cellNum}.objects.outlines)<1 ...
                || length(cellList.meshData{frame}{cellNum}.objects.outlines)>2 ...
                || cellList.meshData{frame}{cellNum}.length>divisionLength ...
                || cellList.meshData{frame}{cellNum}.length<birthLength
            continue
        end
        %Stetp 1.2: calculate cell width, in microns, along cell length
        cellList.meshData{frame}{cellNum}.width = sqrt((cellList.meshData{frame}{cellNum}.mesh(:,4)-cellList.meshData{frame}{cellNum}.mesh(:,2)).^2+(cellList.meshData{frame}{cellNum}.mesh(:,3)-cellList.meshData{frame}{cellNum}.mesh(:,1)).^2);
        cellList.meshData{frame}{cellNum}.width = umperpixel.*cellList.meshData{frame}{cellNum}.width;
        
        %Stetp 1.3: calculate cell width at midcell (as minimum width in
        %middle three width measurements
        cellList.meshData{frame}{cellNum}.septumWidth = min(cellList.meshData{frame}{cellNum}.width(ceil(end/2)-1:ceil(end/2)+1));
        
        %Step 1.4: put results in the matrix lengthAndSeptum Width
        cellCount=cellCount+1;
        lengthAndSeptumWidthEstTime(cellCount,1) = umperpixel*cellList.meshData{frame}{cellNum}.length;
        lengthAndSeptumWidthEstTime(cellCount,2) = cellList.meshData{frame}{cellNum}.septumWidth;
        
        %Step 1.5: estimated time from cell birth using the following
        %formula and put the results in column 3:
        %length at time x = (length at birth) * 2^x/tau
        %2^x/tau = (length at time x)/(length at birth)
        %x/tau = log2((length at time x)/(length at birth))
        %x = log2((length at time x)/(length at birth)) * tau
        lengthAndSeptumWidthEstTime(cellCount,3) = log2((lengthAndSeptumWidthEstTime(cellCount,1)/(birthLength*umperpixel)))*tau;

        
        %Step 2.5: sort cells by cell length
        lengthAndSeptumWidthEstTime=sortrows(lengthAndSeptumWidthEstTime);
    end
end

%Step 4: bin the data
[~,~,binIndex]=histcounts(lengthAndSeptumWidthEstTime(:,1),binEdges);
binnedLengthSeptumWidthEstTime(:,1)=accumarray(binIndex,lengthAndSeptumWidthEstTime(:,1),[],@mean);
binnedLengthSeptumWidthEstTime(:,2)=accumarray(binIndex,lengthAndSeptumWidthEstTime(:,2),[],@mean);
binnedLengthSeptumWidthEstTime(:,3)=accumarray(binIndex,lengthAndSeptumWidthEstTime(:,3),[],@mean);

%Step 5: calculate the standard deviation and standard error of the binned
%mid cell width measurements
binnedLengthSeptumWidthEstTime(:,4)=accumarray(binIndex,lengthAndSeptumWidthEstTime(:,2),[],@length); 
binnedLengthSeptumWidthEstTime(:,5)=accumarray(binIndex,lengthAndSeptumWidthEstTime(:,2),[],@std);
binnedLengthSeptumWidthEstTime(:,6)=binnedLengthSeptumWidthEstTime(:,5)./sqrt(binnedLengthSeptumWidthEstTime(:,4));

%Step 6: remove empty rows from matrix
binnedLengthSeptumWidthEstTime=binnedLengthSeptumWidthEstTime(any(binnedLengthSeptumWidthEstTime,2),:);
end



