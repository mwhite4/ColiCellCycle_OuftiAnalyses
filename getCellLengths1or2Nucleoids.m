function [cellLengthList] = getCellLengths1or2Nucleoids(cellList,umperpixel)


%% This function takes the structure cellList, and puts the lengths of all 
% cells (converted to micrometers) with either 1 or 2 nucleoids (defined as 
%Oufti segmented 'Objects') in the output vector cellLengthList.

%Martin White 2019

%Input
%cellList: cellList output of microbeTracker, or Oufti after running
%function getExtraDataLoop (the structure cellList must carry the field
%'length'

%umperpixel: conversion factor for pixel to micrometers.




%% 

% Step 1: make a matrix to contain cell lengths, frame number is columns, cells are
% rows
numberOfFrames = length(cellList.meshData);
for frameNumber = 1:length(cellList.meshData)
    cellsPerFrame(frameNumber)=length(cellList.meshData{frameNumber});
end
maxNumberofCellsperFrame = max(cellsPerFrame);
cellLengthMatrix(maxNumberofCellsperFrame,numberOfFrames)=0;

%Step 2: populate this matrix (cellLengthMatrix) with measured cell lengths

for frameNumber = 1:length(cellList.meshData)
    for cellNumber = 1:length(cellList.meshData{frameNumber})
        %MW: check if the cell is a valid cell
        if isfield(cellList.meshData{frameNumber}{cellNumber},'mesh') && length(cellList.meshData{frameNumber}{cellNumber}.mesh)>=4
            %if cell has either 1 or 2 nucleoids, enter their
            %length into the matrix
            if isfield(cellList.meshData{frameNumber}{cellNumber},'objects')
                if ~isempty(cellList.meshData{frameNumber}{cellNumber}.objects.outlines) && length(cellList.meshData{frameNumber}{cellNumber}.objects.outlines)< 3
                    cellLengthMatrix(cellNumber,frameNumber)=cellList.meshData{frameNumber}{cellNumber}.length*umperpixel;
                else
                    cellLengthMatrix(cellNumber,frameNumber)=0;
                end
            else
                cellLengthMatrix(cellNumber,frameNumber)=0;
            end
        else
            cellLengthMatrix(cellNumber,frameNumber)=0;
        end
    end
end

%Step 3: convert matrix into a single vector and remove 0s
cellLengthList=cellLengthMatrix(:);
cellLengthList=cellLengthList(cellLengthList~=0);

end
