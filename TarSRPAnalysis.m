function [] = TarSRPAnalysis(AnimalName,Date)
% TarSRPAnalysis.m
%
%  Will take data recorded in response to the stimulus created by
%   TargetSRP.m and analyze it. There, a series of increasingly larger Gabor
%   patches are flashed at the mouse. The center of the patch is the center
%   of mass of the retinotopic field calculated by Retinotopy.m and
%   MapRetinotopy.m .
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%
%OUTPUT: saved file, figures
%
% Created: 2016/08/15, 24 Cummington Mall, Boston, MA
%  Byron Price
% Updated: 2016/08/15
%  By: Byron Price

cd ('~/CloudStation/ByronExp/SRP');
set(0,'DefaultFigureWindowStyle','docked');

% read in the .plx file
EphysFileName = sprintf('TargetSRPData%d_%d',Date,AnimalName); % no file identifier
                    % because MyReadall does that for us

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    MyReadall(EphysFileName);
end

StimulusFileName = sprintf('TargetSRP%d_%d.mat',Date,AnimalName);
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName)
load(StimulusFileName)

sampleFreq = adfreq;


% tsevs are the strobed times of stimulus onset, then offset
%  Onset at tsevs{1,33}(2), offset at tsevs{1,33}(3), onset at
%  tsevs{1,33}(4), offset at 5, etc.
% allad contains the continuous data from each channel, which appear to be
%  recorded at 1000 Hz rather than 40,000

%totalAD = size(allad,2);
%totalSEVS = size(tsevs,2);

Chans = find(~cellfun(@isempty,allad));numChans = length(Chans);
strobeStart = 33;

% lowpass filter the data
dataLength = length(allad{1,Chans(1)});

ChanData = zeros(dataLength,numChans);
preAmpGain = 1;
for ii=1:numChans
    voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
    temp = smooth(voltage,0.013*sampleFreq);
    n = 30;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    ChanData(:,ii) = filter(blo,1,temp);
end

timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeTimes = tsevs{1,strobeStart};
stimLen = round((0.2)*sampleFreq); % 200 milliseconds
minWin = round(0.04*sampleFreq):1:round(0.1*sampleFreq);
maxWin = round(.1*sampleFreq):1:round(0.2*sampleFreq);

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
numPhases = round((2*pi)/phaseShift);
Response = zeros(numChans,numStimuli,numRadii,reps,stimLen);
for ii=1:numChans
    for jj=1:numStimuli
        for kk=1:numRadii
            mystr = sprintf('%d%d',jj,kk);
            stimNum = str2double(mystr);
            stimStrobes = strobeTimes(svStrobed == stimNum);
            for ll=1:reps
                stimOnset = stimStrobes(ll);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index:index+stimLen-1,ii);
                Response(ii,jj,kk,ll,:) = temp;
                clear temp;
            end
        end
    end
end

% STATISTICS OF INTEREST are T1 = min(meanResponse), T2 =
% max(meanResponse), T3 = max-min (meanResponse)
% in the interval from 0 to ~ 0.2 seconds after an image is flashed on the 
% screen, this is a measure of the size of a VEP
meanResponse = squeeze(mean(Response,4));
numStats = 3;
dataStats = struct;
for ii=1:numStats
    dataStats(ii).stdError = zeros(numChans,numStimuli,numRadii);
end
dataStats(1).mean = max(meanResponse(:,:,:,maxWin),[],4);
dataStats(2).mean = -min(meanResponse(:,:,:,minWin),[],4);
dataStats(3).mean = max(meanResponse(:,:,:,maxWin),[],4)-min(meanResponse(:,:,:,minWin),[],4);

% BOOTSTRAP FOR STANDARD ERROR OF STATISTICS IN PRESENCE OF VISUAL STIMULI
N = 1000;
for ii=1:numChans
    for jj=1:numStimuli
        for kk=1:numRadii
            Tboot = zeros(N,numStats);
            for mm=1:N
                indeces = random('Discrete Uniform',reps,[reps,1]);
                group = squeeze(Response(ii,jj,kk,indeces,:));
                meanGroup = mean(group,1);
                Tboot(mm,1) = max(meanGroup(maxWin));
                Tboot(mm,2) = -min(meanGroup(minWin));
                Tboot(mm,3) = max(meanGroup(maxWin))-min(meanGroup(minWin));
            end
            for nn=1:numStats
                dataStats(nn).stdError(ii,jj,kk) = std(Tboot(:,nn));
            end
        end
    end
end

% BOOTSTRAP FOR 95% CONFIDENCE INTERVALS OF STATISTIC IN ABSENCE OF VISUAL STIMULI
%  PLUS MEAN AND STANDARD ERRORS
%  interspersed stimulus repetitions with holdTime seconds of a blank grey
%  screen
N = 1000; % number of bootstrap samples
noStimLen = holdTime*sampleFreq-stimLen*2;

baselineStats = struct;
for ii=1:numStats
    baselineStats(ii).prctile = zeros(numChans,1);
    baselineStats(ii).mean = zeros(numChans,1);
    baselineStats(ii).stdError = zeros(numChans,1);
end

for ii=1:numChans
    Tboot = zeros(N,numStats);
    pauseOnset = strobeTimes(svStrobed == 0);nums = length(pauseOnset);
    for jj=1:N
        indeces = random('Discrete Uniform',noStimLen,[reps,1]);
        temp = zeros(reps,stimLen);
        num = random('Discrete Uniform',nums);
        [~,index] = min(abs(timeStamps-pauseOnset(num)));
        for kk=1:reps
            temp(kk,:) = ChanData(index+indeces(kk):index+indeces(kk)+stimLen-1,ii);
        end
        meanTrace = mean(temp,1);
        Tboot(jj,1) = max(meanTrace(maxWin));
        Tboot(jj,2) = -min(meanTrace(minWin));
        Tboot(jj,3) = max(meanTrace(maxWin))-min(meanTrace(minWin));
    end
    %figure();histogram(Tboot);
    for ll=1:numStats
        baselineStats(ll).prctile(ii) = quantile(Tboot(:,ll),1-1/100);
        baselineStats(ll).mean(ii) = mean(Tboot(:,ll));
        baselineStats(ll).stdError(ii) = std(Tboot(:,ll));
    end
end


% WALD TEST - VEP magnitude significantly greater in presence of a stimulus
%  than in the absence of a stimulus
significantStimuli = zeros(numChans,numStimuli,numRadii,numStats);
alpha = 0.05/(numStimuli*numRadii);
c = norminv(1-alpha,0,1);
for ii=1:numChans
    for jj=1:numStimuli
        for kk=1:numRadii
            for ll=1:numStats
                W = (dataStats(ll).mean(ii,jj,kk)-baselineStats(ll).mean(ii))...
                    /sqrt(dataStats(ll).stdError(ii,jj,kk)^2+baselineStats(ll).stdError(ii)^2);
                if W > c
                    significantStimuli(ii,jj,kk,ll) = dataStats(ll).mean(ii,jj,kk); % or equals W itself
                end
            end
        end
    end    
end

% GRAPHS OF RESPONSES
for ii=1:numChans
    h(ii) = figure();
    for jj=1:numStimuli
        subplot(1,2,jj);
        for kk=1:numStats
            cms = (Radii(jj,:)).*mmPerPixel./10;
            degrees = atan(cms./DistToScreen).*180/pi;
            errorbar(log2(degrees),squeeze(dataStats(kk).mean(ii,jj,:)),...
                squeeze(dataStats(kk).stdError(ii,jj,:)),'LineWidth',2);
            hold on;
        end
        title(Stimulus(jj).name);
        legend('VEP Max','VEP Min','VEP Max-Min');
        xlabel('Stimulus Radius (log2[degrees of visual space])');
        ylabel('VEP Magnitude(\muVolts)')
    end
end
savefig(h,sprintf('TargetSRP%d_%d.fig',Date,AnimalName));
save(sprintf('TargetSRPConv%d_%d.mat',Date,AnimalName),'Response',...
    'dataStats','baselineStats','significantStimuli','Stimulus','Radii',...
    'mmPerPixel','DistToScreen');
end


