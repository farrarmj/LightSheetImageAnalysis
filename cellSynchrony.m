function [peakSize,Responding,simIdx,idx,w,p,lcs] = cellSynchrony(traces,varargin)
%This function takes normalized calcium imaging traces and assesses the
%synchrony of calcium peaks. 
%The first part of the code automatically selects potential peaks for analysis. 
%The second part of the code prompts the user to verify that invidiual
%peaks show appropriate kinetics for calcium imaging, excluding imaging
%artifacts.
%The final part of the code calculates the similarity index of all cells to
%the template cell for each valid event. Cells with a similarlity index
%greater than similarityCutoff are considered to be "Responding". 
%For each cell and each event, a peak size (peakSize), fraction of cells
%responding (Responding), similarity index (simIdx), template cell index
%(idx), event width (w), event prominence (p), and location in the trace
%(lcs) is recorded.

%Heuristically determined threshold for a response to be synchronous
similarityCutoff = 0.7;

%Only consider peaks of at least 3 frames width
minPeakWidth = 3;

%Peak prominence
if nargin<2
    peakProm = 0.2; %Default value of 0.2
else
    peakProm = varargin{1};
end


%Compute "maximum intensity" trace from all individual cell traces
maxTrace = max(traces,[],2);

%Automated prescreen  of calcium peaks based on peak prominence and peak
%width
[vals,lcs,w,p] = findpeaks(maxTrace,'MinPeakProminence',peakProm,'MinPeakWidth',minPeakWidth);

%%Block for manual verification based on kinetics

winSize = 40; %cropping window centered onpeak
n = 0; %initialize counter for valid number of true peaks

%Consider all possible peaks
for j=1:length(vals) 
    ok = 'No';
    while strcmp(ok,'No')
        h = figure;  %create display figure
        
        %Compute peak prominences of each individual cell to find the
        %template cell
        for k=1:size(traces,2)
            [~,~,~,temp] = findpeaks(traces(max(1,lcs(j)-round(w(j))):min(lcs(j)+round(w(j)),length(traces)),k)); %Calculate Peak Prominences
        if ~isempty(temp)
             proms(k) = max(temp);
        else
            proms(k) = 0;
        end
     end
    [~,idx(j)] = max(proms); %Index of cell with higest peak prominence
    templateCell = traces(:,idx(j)); %Template cell trace
    
    %Display template cell trace
    plot(templateCell(max(1,lcs(j)-winSize/2):min(lcs(j)+winSize/2,length(traces)))),ylim([-0.1,vals(j)+0.05])
    
    %Request user input to exclude spurious events
    ok = questdlg('Valid Event?','Event ok?','Yes','No','Yes');
        if strcmp(ok,'Yes') %If the event is valid
            %User crops event to eliminate neighbouring events
            R=getrect;
            X = [R(1),R(1)+R(3),R(1)+R(3),R(1),R(1)];
            Y = [R(2),R(2),R(2)+R(4),R(2)+R(4),R(2)];
            line(X,Y,'Color','r')
            %Double check that user selection is appropriate
            ok = questdlg('Event selection ok?','Event','Yes','No','Yes');
        else %
            vals(j) = NaN;
            lcs(j) = NaN;
            w(j) = NaN;
            p(j) = NaN;
            break
        end
    end
    close(h) %close figure 
    
    %If peak is valid, crop window including peak from original data set
    if ~isnan(vals(j))
        n = n+1; %update number of valid peaks
        win(n) = {max(1,lcs(j)-winSize/2)+floor(R(1)):max(1,lcs(j)-winSize/2)+ceil(R(1)+R(3))};
    end
end
%Keep all valid peaks 
 lcs = lcs(~isnan(lcs));      
 peakSize = vals(~isnan(vals));
 w = w(~isnan(w));
 p = p(~isnan(p));
 
 
 %%Block to compute synchronicity of valid peaks
 for j=1:length(peakSize)
     for k=1:size(traces,2)
        [~,~,~,temp] = findpeaks(traces(max(1,lcs(j)-round(w(j))):min(lcs(j)+round(w(j)),length(traces)),k)); %Calculate Peak Prominences
        if ~isempty(temp)
             proms(k) = max(temp);
        else
            proms(k) = 0;
        end
     end
    [~,idx(j)] = max(proms); %ID of cell with higest peak prominence
    templateCell = traces(:,idx(j)); %template cell with highest peak value
    for k=1:size(traces,2)
        comparisonCell = traces(win{j},k); %Cell trace to compare against
        temp = xcorr(templateCell(win{j})-min(templateCell(win{j})),comparisonCell-min(comparisonCell),'Coeff'); %Normalized Cross Correlation
        similarity(k) = temp((length(temp)+1)/2); %similarity index at shift of 0 

   end
    simIdx(j) = {similarity}; %list of similarity values
    Responding(j) = sum(similarity>similarityCutoff)/size(traces,2); %Fraction of responding cells
 end
 pause(0.5);