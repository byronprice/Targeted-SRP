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
% Updated: 2016/08/29
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
    n = 30;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    ChanData(:,ii) = filter(blo,1,voltage);
end

timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeTimes = tsevs{1,strobeStart};
stimLen = round(0.3*sampleFreq); % 200 milliseconds
minWin = round(0.05*sampleFreq):1:round(0.15*sampleFreq);
maxWin = round(.1*sampleFreq):1:round(0.25*sampleFreq);
smoothKernel = 4;

reps = reps-1;
% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
numPhases = round((2*pi)/phaseShift);
Response = zeros(numChans,numStimuli*numRadii,reps,stimLen);
meanResponse = zeros(numChans,numStimuli*numRadii,stimLen);
for ii=1:numChans
    for jj=1:numStimuli*numRadii
        stimStrobes = strobeTimes(svStrobed == jj);
        for ll=1:reps
            stimOnset = stimStrobes(ll);
            [~,index] = min(abs(timeStamps-stimOnset));
            temp = ChanData(index:index+stimLen-1,ii);
            Response(ii,jj,ll,:) = temp;
            clear temp;
        end
        temp = mean(squeeze(Response(ii,jj,:,:)),1);
        meanResponse(ii,jj,:) = smooth(temp,smoothKernel);
    end
end

% STATISTICS OF INTEREST are T1 = min(meanResponse), T2 =
% max(meanResponse), T3 = max-min (meanResponse)
% in the interval from 0 to ~ 0.2 seconds after an image is flashed on the 
% screen, this is a measure of the size of a VEP
numStats = 3;
dataStats = struct;
dataStats(1).name = 'VEP Positivity';dataStats(2).name = 'VEP Negativity';dataStats(3).name = 'VEP Total';
dataStats(1).specs = '--^k';dataStats(2).specs = ':vr';dataStats(3).specs = '-+c';
for ii=1:numStats
    dataStats(ii).stdError = zeros(numChans,numStimuli*numRadii);
    dataStats(ii).latencySEM = zeros(numChans,numStimuli*numRadii);
end
[dataStats(1).mean,dataStats(1).latency] = max(meanResponse(:,:,maxWin),[],3);
[mins,dataStats(2).latency] = min(meanResponse(:,:,minWin),[],3);
dataStats(2).mean = -mins;
dataStats(3).mean = max(meanResponse(:,:,maxWin),[],3)-min(meanResponse(:,:,minWin),[],3);

% BOOTSTRAP FOR STANDARD ERROR OF STATISTICS IN PRESENCE OF VISUAL STIMULI
N = 1000;
for ii=1:numChans
    for jj=1:numStimuli*numRadii
        Tboot = zeros(N,numStats);
        Latency = zeros(N,numStats);
        for mm=1:N
            indeces = random('Discrete Uniform',reps,[reps,1]);
            group = squeeze(Response(ii,jj,indeces,:));
            meanGroup = mean(group,1);
            [Tboot(mm,1),Latency(mm,1)] = max(meanGroup(maxWin));
            [mins,Latency(mm,2)] = min(meanGroup(minWin));
            Tboot(mm,2) = -mins;
            Tboot(mm,3) = max(meanGroup(maxWin))-min(meanGroup(minWin));
        end
        for nn=1:numStats
            dataStats(nn).stdError(ii,jj) = std(Tboot(:,nn));
            dataStats(nn).latencySEM(ii,jj) = std(Latency(:,nn));
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
    baselineStats(ii).latency = zeros(numChans,1);
    baselineStats(ii).latencySEM = zeros(numChans,1);
end

pauseOnset = strobeTimes(svStrobed == 0);nums = length(pauseOnset);
for ii=1:numChans
    Tboot = zeros(N,numStats);
    Latency = zeros(N,numStats);
    for jj=1:N
        indeces = random('Discrete Uniform',noStimLen,[reps,1]);
        temp = zeros(reps,stimLen);
        num = random('Discrete Uniform',nums);
        [~,index] = min(abs(timeStamps-pauseOnset(num)));
        index = index+stimLen;
        for kk=1:reps
            temp(kk,:) = ChanData(index+indeces(kk):index+indeces(kk)+stimLen-1,ii);
        end
        meanTrace = mean(temp,1);
        [Tboot(jj,1),Latency(jj,1)] = max(meanTrace(maxWin));
        [mins,Latency(jj,2)] = min(meanTrace(minWin));
        Tboot(jj,2) = -mins;
        Tboot(jj,3) = max(meanTrace(maxWin))-min(meanTrace(minWin));
    end
    %figure();histogram(Tboot);
    for ll=1:numStats
        baselineStats(ll).prctile(ii) = quantile(Tboot(:,ll),1-1/100);
        baselineStats(ll).mean(ii) = mean(Tboot(:,ll));
        baselineStats(ll).stdError(ii) = std(Tboot(:,ll));
        baselineStats(ll).latency(ii) = mean(Latency(:,ll));
        baselineStats(ll).latencySEM(ii) = std(Latency(:,ll));
    end
end


% WALD TEST - VEP magnitude significantly greater in presence of a stimulus
%  than in the absence of a stimulus
significantStimuli = zeros(numChans,numStimuli*numRadii,numStats);
alpha = 0.05/(numStimuli*numRadii);
c = norminv(1-alpha,0,1);
for ii=1:numChans
    for jj=1:numStimuli*numRadii
        for ll=1:numStats
            W = (dataStats(ll).mean(ii,jj)-baselineStats(ll).mean(ii))...
                /sqrt(dataStats(ll).stdError(ii,jj)^2+baselineStats(ll).stdError(ii)^2);
            if W > c
                significantStimuli(ii,jj,ll) = dataStats(ll).mean(ii,jj); % or equals W itself
            end
        end
    end    
end

if Channel == 1
    chans = {'Target','Off-Target'};
else
    chans = {'Off-Target','Target'};
end

maxVEP = max(max(dataStats(3).mean));
% GRAPHS OF RESPONSES
h = figure();
count = 1;order = [1,2,5,6,3,4,7,8];
for ii=1:numChans
    for jj=1:numStimuli
        figure(h);subplot(4,2,order(count));
        members = Radii(:,2) == jj;
        for kk=1:numStats
            cms = (Radii(members,1)).*mmPerPixel./10;
            degrees = atan(cms./DistToScreen).*180/pi;
            degrees = degrees';
            [sortedDegs,indeces] = sort(degrees);
            Stats = squeeze(dataStats(kk).mean(ii,members));
            sortedStats = Stats(indeces);
            Errors = squeeze(dataStats(kk).stdError(ii,members));
            sortedErrors = Errors(indeces);
            errorbar(log2(sortedDegs),sortedStats,...
                sortedErrors,dataStats(kk).specs,...
                'LineWidth',2);
            hold on;
        end
        title(sprintf('%s, Channel: %s, Animal: %d',Stimulus(jj).name,chans{ii},AnimalName));
        legend(dataStats(1).name,dataStats(2).name,dataStats(3).name);
        xlabel('Radius (degrees of visual space [log scale])');
        ylabel('VEP Magnitude(\muVolts)')
        axis([(log2(sortedDegs(1))-1) (log2(sortedDegs(end))+1) -50 500]);
        ax = gca;
        ax.XTick = [log2(sortedDegs(1))-1,log2(sortedDegs),log2(sortedDegs(end))+1];
        ax.XTickLabel = [2^(log2(sortedDegs(1))-1),sortedDegs,2^(log2(sortedDegs(end))+1)];
        count = count+1;
   
        figure(h);subplot(4,2,order(count));
        VEP = squeeze(meanResponse(ii,members,:));
        sortedVEP = VEP(indeces,:);
        shift = 250;
        for kk=1:length(indeces)
            plot(1:stimLen,sortedVEP(kk,:)+shift*(kk),'LineWidth',2);
            hold on;
        end
        title(sprintf('%s, Channel: %s, Animal: %d',Stimulus(jj).name,chans{ii},AnimalName));
        xlabel('Time from phase shift (milliseconds)');
        ylabel('Radius (degrees of visual space [log scale])');
        axis([0 stimLen 0 (shift*length(indeces)+shift)]);
        ax = gca;
        ax.YTick = linspace(0,shift*length(indeces)+shift,length(indeces)+2);
        ax.YTickLabel = [2^(log2(sortedDegs(1))-1),sortedDegs,2^(log2(sortedDegs(end))+1)];
        count = count+1;
        
        clear indeces;
    end
end
savefig(h,sprintf('TargetSRP%d_%d.fig',Date,AnimalName));
save(sprintf('TargetSRPConv%d_%d.mat',Date,AnimalName),'Response',...
    'dataStats','baselineStats','significantStimuli','Stimulus','Radii',...
    'mmPerPixel','DistToScreen','meanResponse');
end


