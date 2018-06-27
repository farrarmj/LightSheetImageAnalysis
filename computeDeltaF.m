function R = computeDeltaF(trace,varargin)
%This function computes the ratio, R = delta(F)/F = (F - F0)/F0 
%for calcium signals. The raw intensity images are asssumed to have been
%processed to give the cell intensity vs time in the variable "trace"

%Parameters for measuring baseline fluorescence, F0
smoothSpan = 10; %smooth the trace over smoothSpan data points
minSpan = 100; %consider a moving window of 100 frames
cutoffPrctile = 30; %noise floor will be 30th percentile of window

%Option to use local or global baseline, F0
%Global option uses a single value of F0, and does not allow for
%photobleach or image drift

if nargin>1
    method = varargin{1};
    if isempty(method) || strcmp(method,'local')
        method = 'local';
    elseif strcmp(method,'global')
        method = 'global';
    end
    else
    method = 'local';
end

%%Block to compute F0
    Fbar = smooth(trace,smoothSpan);    %Create smoothed average
    if strcmp(method,'local')
        for k=1:length(Fbar)
                %Create window that fits in vector
                lowBound = k-minSpan/2;
                upperBound = k+minSpan/2;
                if lowBound <=1
                    lowBound = 2;
                    upperBound = minSpan+1;
                elseif upperBound>length(Fbar)
                    upperBound = length(Fbar);
                    lowBound = length(Fbar)-minSpan;
                end
                %Create F0, the local background
                F0(k) = prctile(Fbar(lowBound:upperBound),cutoffPrctile); 
        end
    else
        %Global option
        F0(1:length(Fbar)) = prctile(Fbar,cutoffPrctile);
    end
     
    
    R = (trace-F0')./F0'; %Compute ratio, R
    
end