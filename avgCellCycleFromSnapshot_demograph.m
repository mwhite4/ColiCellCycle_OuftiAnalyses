function avgCellCycleFromSnapshot_demograph(varargin) 
%-----------------------------------------------------------------------------------------------------
%% MW: This is an edit to the OUFTI function 'demograph'
%Martin White 2018


%% 
%-----------------------------------------------------------------------------------------------------
%function demograph(varargin)
%
%@author:  Jason Hocking
%@date:    October 27, 2011
%@modified:  Ahmad Paintdakhi -- August 1, 2013
%@copyright 2011-2013 Yale University
%====================================================================================================
%**********output********:
%No output arguements required, the function oly plots information.
%**********Input********:
%cellList:	cellList structure
%maxCellNum:	maximum number of cells to be included for final demograph
%maxCellLength:	maximum length of cell
%numPixelsMovingAverage:	number of pixels to be used for the moving average.  This
%							routine finds the segment where max intensity of the signal
%							is located.
%signal:	an array for signal information.  For example, to use only signal 1 the
%			array should be [1,0], for signal 2 --> [0,1] and both signal 1 and
%			signal 2 ---> [1,1].
%frameNum:  frame # to be used for analysis or [] vector (use all frames in
%           a dataset).
%
%descriptor:	The descriptor value is a key for the type of demograph to be drawn.
%				The different keys are 'randomN','randomNOriented', 'constriction_noNormalization',
%				'sort_by_constriction','constriction','normByPopulation', and 'normByPopulationOriented'.
%conversionFactor:  Pixel to micron conversion factor.
%Purpose:  script was designed to provide a colormap of relative segment intensities for every
% 		   cell in an asynchorous cellList, sorted by cell length in ascending order.
%====================================================================================================

if length(varargin) < 6 || length(varargin) > 6
    disp('A total of 6 arguments are accepted');
    return;
end
if ~isstruct(varargin{1}) && ~iscell(varargin{1})
    disp('cellList must be a struct or cell array')
    return;
end

if ~isscalar(varargin{6})
    disp('Conversion factor needs to be a scalar value such as 0.064')
    return;
end
cellList = varargin{1};
maxCellNum = varargin{2};
minCellLength=varargin{3};
maxCellLength = varargin{4};
frameNum = varargin{5};
conversionFactor = varargin{6};
warning('off','MATLAB:colon:nonIntegerIndex');

Lblock=60;

%---------------------------------------------------------------------------------
signalInfo = 'signal1';

if isempty(frameNum)
    frameList = 1:length(cellList.meshData);
else
    frameList = frameNum;
end

try
%     replacement=false;

    %MW: calculates the total number of cells shorter than the input max
    %length (=n)
    %%finds the maximum number of stepareas inside of a cell from the cellList
    maxsizelarray=[];
    n=0;
    for frame = frameList
        for cellNum = 1:length(cellList.meshData{frame})
            if isempty(cellList.meshData{frame}{cellNum}) ...
                    || ~isfield(cellList.meshData{frame}{cellNum},'mesh') ...
                    ||length(cellList.meshData{frame}{cellNum}.mesh)<4 ...
                    ||~isfield(cellList.meshData{frame}{cellNum},signalInfo) ...
                    || eval('isempty(cellList.meshData{frame}{cellNum}.(signalInfo))') ...
                    || cellList.meshData{frame}{cellNum}.length>maxCellLength ...
                    || cellList.meshData{frame}{cellNum}.length<minCellLength ...
                    || length(cellList.meshData{frame}{cellNum}.objects.outlines)<1 ...
                    || length(cellList.meshData{frame}{cellNum}.objects.outlines)>2
                continue
            end
            n=n+1;
        end
    end
    if n<=maxCellNum
        maxCellNum=n;
    end
    
    %rand is a list of length maxCellNum, containing numbers from 1:the
    %total number of cells shorter than maxCellLength.  If the number of
    %cells shorter than maxCellLength is less than maxCellNum then all such
    %cells are used
    rand=randsample(n,maxCellNum,false);
    
    n=0;
    for frame = frameList
        for cellNum = 1:length(cellList.meshData{frame})
            if isempty(cellList.meshData{frame}{cellNum}) ...
                    || length(cellList.meshData{frame}{cellNum}.mesh)<4 ...
                    || ~isfield(cellList.meshData{frame}{cellNum},signalInfo) ...
                    || eval('isempty(cellList.meshData{frame}{cellNum}.(signalInfo))') ...
                    || cellList.meshData{frame}{cellNum}.length>maxCellLength ...
                    || cellList.meshData{frame}{cellNum}.length<minCellLength ...
                    || length(cellList.meshData{frame}{cellNum}.objects.outlines)<1 ...
                    || length(cellList.meshData{frame}{cellNum}.objects.outlines)>2
                continue
            end
            n = n+1;
            b=rand==n;
            
            if sum(b)~=1
                continue
            end
            maxsizelarray=[maxsizelarray length(cellList.meshData{frame}{cellNum}.lengthvector)];%#ok<AGROW>
        end
    end
    if isempty(maxsizelarray)
        warndlg(['No field ' signalInfo ' recorded for this cell:  Use Reuse meshes toggle button to compute ' signalInfo]);
        return;
    end
    
    %using the maxima from above, a matrix consiting of zeros is created to be
    %filled in by mesh intensities
    %MW: create a matrix (# of rows = max number of steps, # columns = max number of cells) and fill it with zeros
    relintarray1=zeros(max(maxsizelarray),maxCellNum+10);
    maxsizel=max(maxsizelarray);
    if maxCellLength > maxsizel
        maxCellLength = maxsizel;
    end
    
    maxsizel2 = ceil(maxsizel); if mod(maxsizel2,2)==0, maxsizel2=maxsizel2+1; end
    %MW: estimated(?) center of longest allowed cell
    maxsizel2a = maxsizel2/2+0.5;
    n=0;
    passed=0;
    cellLength=[];
    
    %zeroarray is replaced with relative segment intensity data from the cell
    for frame = frameList
        for cellNum = 1:length(cellList.meshData{frame})

            if isempty(cellList.meshData{frame}{cellNum}) || length(cellList.meshData{frame}{cellNum}.mesh)<4 ...
                    ||~isfield(cellList.meshData{frame}{cellNum},signalInfo) ...
                    || eval('isempty(cellList.meshData{frame}{cellNum}.(signalInfo))') || cellList.meshData{frame}{cellNum}.length>maxCellLength...
                    || cellList.meshData{frame}{cellNum}.length<minCellLength ...
                    || length(cellList.meshData{frame}{cellNum}.objects.outlines)<1 ...
                    || length(cellList.meshData{frame}{cellNum}.objects.outlines)>2
                continue
            end
            n = n+1;
            b=rand==n;
            if sum(b)~=1
                continue
            end
            passed=passed+1;
            
            %MW: due to background subtraction, some segments end up with a
            %signal value of 0 (I don't know if negative values are
            %possible), since 0s get converted to NaNs I give them an
            %arbitrary small value
            for i=1:length(cellList.meshData{frame}{cellNum}.signal1)
                if cellList.meshData{frame}{cellNum}.signal1(i)<=0
                    cellList.meshData{frame}{cellNum}.signal1(i)=0.000001;
                end
            end
            
%% relint 1 are the intensity values that will be plotted, this is where you put normalization, if any.
           cellList.meshData{frame}{cellNum}.relint1 =cellList.meshData{frame}{cellNum}.signal1;
           
%Percent of total cellular signal
%             cellList.meshData{frame}{cellNum}.relint1 =cellList.meshData{frame}{cellNum}.signal1/sum(cellList.meshData{frame}{cellNum}.signal1)*100;
%normalize to maximum value (this is what Oufti uses
            %             cellList.meshData{frame}{cellNum}.relint1 = (cellList.meshData{frame}{cellNum}.signal1./max(cellList.meshData{frame}{cellNum}.signal1));
%this is a weird fraction of total signal multiplied by cell length
            %             cellList.meshData{frame}{cellNum}.relint1 =cellList.meshData{frame}{cellNum}.signal0/sum(cellList.meshData{frame}{cellNum}.signal0)*length(cellList.meshData{frame}{cellNum}.signal0);
%%           
            k = floor(cellList.meshData{frame}{cellNum}.length/2);
            %MW: temp = distance of step from midcell
            temp = cellList.meshData{frame}{cellNum}.lengthvector-cellList.meshData{frame}{cellNum}.length/2;
            %%MW: 1D linear interpolation, I think this is to take care of
            %%the issue of the steps not being equal lengths
            interpint1 = interp1(temp(1:length(cellList.meshData{frame}{cellNum}.relint1)),cellList.meshData{frame}{cellNum}.relint1,-k:k,'linear','extrap');
            
            %MW: maxsizel2a-k:maxsizel2a+k, centers the cell over the
            %center of the longest cell, 'passed' must be current cell
            relintarray1(maxsizel2a-k:maxsizel2a+k,passed)=interpint1;
            cellLength=[cellLength cellList.meshData{frame}{cellNum}.length]; %#ok<AGROW>
        end
    end
    
    %MW: I want to spike in 5 extra fake cells of fixed length (length at initiation?)and no
    %fluorescence to act as a reference size (should get a white line)
    cellLength=[cellLength [Lblock Lblock Lblock Lblock Lblock Lblock Lblock Lblock Lblock Lblock]];

    % % cells length array is concatonated with the fluorescence matrix. This matrix is then sorted by length in ascending order
    numlist=[1:1:maxCellNum+10]; %#ok
    lvint0=horzcat(numlist',cellLength');
    lvint1=horzcat(lvint0,relintarray1');
    lnumsort1=sortrows(lvint1,[2]); %#ok
%     assignin('base','testWeightedorig',lnumsort1);
    %% MW I want to subsample with weighted bias (inverse of population age structure)
    a=(0:1/(maxCellNum-1):1);
    agestructure=2*log(2)*exp(-a*log(2));
    weights=1./agestructure;
    rand2=randsample(maxCellNum,800,true,weights);
    rand2=sort(rand2); 
%     assignin('base','rand2',rand2);
    lnumsort2=zeros(length(rand2),size(lnumsort1,2));
    for j=1:length(rand2)
        lnumsort2(j,:)=lnumsort1(rand2(j),:);
    end
%     assignin('base','testWeighted',lnumsort2);
    
    %     I think that I would have to spike in the reference cells here
    %      cellLength=[cellLength [refLength refLength refLength refLength refLength refLength refLength refLength refLength refLength]];
    
    
    
    %%
%     x=[-conversionFactor*maxCellLength./2 conversionFactor*maxCellLength./2];
%     x = repmat(x(1):x(2)*2/(size(lnumsort1,2)-3):x(2),size(lnumsort1,1),1);
%     y = repmat((1:size(lnumsort1,1)),size(lnumsort1,2)-2,1)';
    x=[-conversionFactor*maxCellLength./2 conversionFactor*maxCellLength./2];
    x = repmat(x(1):x(2)*2/(size(lnumsort2,2)-3):x(2),size(lnumsort2,1),1);
    y = repmat((1:size(lnumsort2,1)),size(lnumsort2,2)-2,1)';
    
    
    %MW: dataToPlot = remove first two columns, leaving only signal
    %information
%     dataToPlot = lnumsort1(1:end,3:end);
    dataToPlot = lnumsort2(1:end,3:end);
    dataToPlot(dataToPlot==0) = NaN;
    figure
    pcolor(x,y,flipud(dataToPlot)); colormap plasma;shading flat;caxis([0 0.4])%10 %MW: add ;colorbar if you want to plot the colorbar
    xlim([-5 5])
    xticks(-5:1:5);
    set(gca,'TickDir','both')
    set(gca,'xticklabel',[])
    yticks([]);
    set(gca,'fontsize',8)
    set(gca,'LineWidth',1)

    
catch err
    if strcmpi(err.identifier,'MATLAB:catenate:dimensionMismatch')
        warndlg('Choose a smaller number for max cell number parameter');
        return;
    end
end


end







