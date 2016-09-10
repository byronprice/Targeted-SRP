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
% Updated: 2016/09/08
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
LatWin = round(0.04*sampleFreq):round(0.15*sampleFreq);
minWin = round(0.04*sampleFreq):round(0.19*sampleFreq);
maxWin = round(.1*sampleFreq):1:round(0.25*sampleFreq);
smoothKernel = 4;
alpha = 0.05;

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

LatencyStats = struct;
[~,temp] = min(Response(:,:,:,LatWin),[],4);
LatencyStats.mag = mean(Response(:,:,:,round(0.02*sampleFreq):round(0.04*sampleFreq)),4)-min(Response(:,:,:,LatWin),[],4);
LatencyStats.trials = (temp+LatWin(1)-1)./sampleFreq;
LatencyStats.mean = zeros(numChans,numStimuli*numRadii);
LatencyStats.sem = zeros(numChans,numStimuli*numRadii);
LatencyStats.ci = zeros(numChans,numStimuli*numRadii,2);
LatFun = @(x,win) max(abs(eCDF(x.*sampleFreq,win)-linspace(0,1,length(win))')); %mad(x,1);

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
    dataStats(ii).ci = zeros(numChans,numStimuli*numRadii,2);
end
[dataStats(1).mean,maxLats] = max(meanResponse(:,:,maxWin),[],3);
dataStats(1).latency = (maxLats+maxWin(1)-1)./sampleFreq;
[mins,minLats] = min(meanResponse(:,:,minWin),[],3);
dataStats(2).latency = (minLats+minWin(1)-1)./sampleFreq;
dataStats(2).mean = -mins;
dataStats(3).mean = max(meanResponse(:,:,maxWin),[],3)-min(meanResponse(:,:,minWin),[],3);
dataStats(3).latency = dataStats(1).latency-dataStats(2).latency;

% BOOTSTRAP FOR STANDARD ERROR OF STATISTICS IN PRESENCE OF VISUAL STIMULI
N = 2000;
for ii=1:numChans
    for jj=1:numStimuli*numRadii
        Tboot = zeros(N,numStats);
        Latency = zeros(N,numStats);
        
        LatBoot = zeros(N,1);
        LatencyStats.mean(ii,jj) = LatFun(squeeze(LatencyStats.trials(ii,jj,:)),LatWin);
        for mm=1:N
            indeces = random('Discrete Uniform',reps,[reps,1]);
            LatGroup = squeeze(LatencyStats.trials(ii,jj,indeces));
            LatBoot(mm) = LatFun(LatGroup,LatWin);
            group = squeeze(Response(ii,jj,indeces,:));
            meanGroup = mean(group,1);
            [Tboot(mm,1),maxLats] = max(meanGroup(maxWin));
            Latency(mm,1) = (maxLats+maxWin(1)-1)./sampleFreq;
            [mins,minLats] = min(meanGroup(minWin));
            Latency(mm,2) = (minLats+minWin(1)-1)./sampleFreq;
            Tboot(mm,2) = -mins;
            Tboot(mm,3) = max(meanGroup(maxWin))-min(meanGroup(minWin));
            Latency(mm,3) = Latency(mm,1)-Latency(mm,2);
        end
%         figure();subplot(2,1,1);histogram(squeeze(LatencyStats.trials(ii,jj,:)));subplot(2,1,2);plot(squeeze(meanResponse(ii,jj,:)));
        LatencyStats.sem(ii,jj) = std(LatBoot);
        LatencyStats.ci(ii,jj,:) = [quantile(LatBoot,alpha/2),quantile(LatBoot,1-alpha/2)];
        for nn=1:numStats
            dataStats(nn).stdError(ii,jj) = std(Tboot(:,nn));
            dataStats(nn).latencySEM(ii,jj) = std(Latency(:,nn));
            dataStats(nn).ci(ii,jj,:) = [quantile(Tboot(:,nn),alpha/2),quantile(Tboot(:,nn),1-alpha/2)];
        end
    end
end

% BOOTSTRAP FOR 95% CONFIDENCE INTERVALS OF STATISTIC IN ABSENCE OF VISUAL STIMULI
%  PLUS MEAN AND STANDARD ERRORS
%  interspersed stimulus repetitions with holdTime seconds of a blank grey
%  screen
noStimLen = holdTime*sampleFreq-stimLen*2;

baseLatStats = struct;
baseLatStats.mean = zeros(numChans,1);
baseLatStats.sem = zeros(numChans,1);
baseLatStats.ci = zeros(numChans,2);

baseStats = struct;
for ii=1:numStats
    baseStats(ii).ci = zeros(numChans,2);
    baseStats(ii).mean = zeros(numChans,1);
    baseStats(ii).stdError = zeros(numChans,1);
    baseStats(ii).latency = zeros(numChans,1);
    baseStats(ii).latencySEM = zeros(numChans,1);
end

pauseOnset = strobeTimes(svStrobed == 0);nums = length(pauseOnset);
for ii=1:numChans
    Tboot = zeros(N,numStats);
    Latency = zeros(N,numStats);
    LatBoot = zeros(N,1);
    for jj=1:N
        indeces = random('Discrete Uniform',noStimLen,[reps,1]);
        temp = zeros(reps,stimLen);
        num = random('Discrete Uniform',nums);
        [~,index] = min(abs(timeStamps-pauseOnset(num)));
        index = index+stimLen;
        for kk=1:reps
            temp(kk,:) = ChanData(index+indeces(kk):index+indeces(kk)+stimLen-1,ii);
        end
        [~,minLats] = min(temp(:,LatWin),[],2);
        minLats = (minLats+LatWin(1)-1)./sampleFreq;
        LatBoot(jj) = LatFun(minLats,LatWin);
        meanTrace = mean(temp,1);
        [Tboot(jj,1),maxLats] = max(meanTrace(maxWin));
        Latency(jj,1) = (maxLats+maxWin(1)-1)./sampleFreq;
        [mins,minLats] = min(meanTrace(minWin));
        Latency(jj,2) = (minLats+minWin(1)-1)./sampleFreq;
        Tboot(jj,2) = -mins;
        Tboot(jj,3) = max(meanTrace(maxWin))-min(meanTrace(minWin));
    end
    baseLatStats.mean(ii) = mean(LatBoot);
    baseLatStats.sem(ii) = std(LatBoot);
    baseLatStats.ci(ii,:) = [quantile(LatBoot,alpha/2),quantile(LatBoot,1-alpha/2)];
    %figure();histogram(Tboot);
    for ll=1:numStats
        baseStats(ll).ci(ii,:) = [quantile(Tboot(:,ll),alpha/2),quantile(Tboot(:,ll),1-alpha/2)];
        baseStats(ll).mean(ii) = mean(Tboot(:,ll));
        baseStats(ll).stdError(ii) = std(Tboot(:,ll));
        baseStats(ll).latency(ii) = mean(Latency(:,ll));
        baseStats(ll).latencySEM(ii) = std(Latency(:,ll));
    end
end


% CI test (something like a non-parametric Wald test) ... for the deviation
%  from a uniform distribution (of latency to peak negativity times)

% CI test ... for VEP magnitude, given 95% confidence intervals
significantStimuli = zeros(numChans,numStimuli*numRadii,2);
c = norminv(1-alpha,0,1);
for ii=1:numChans
    for jj=1:numStimuli*numRadii
        if LatencyStats.ci(ii,jj,1) > baseLatStats.ci(ii,2)
            significantStimuli(ii,jj,1) = 1;
        end
%         W = (dataStats(2).mean(ii,jj)-baseStats(2).mean(ii))/...
%             sqrt(dataStats(2).stdError(ii,jj)^2+baseStats(2).stdError(ii)^2);
%         if W > c
%             significantStimuli(ii,jj,2) = 1;
%         end
        if dataStats(3).ci(ii,jj,1) > baseStats(3).ci(ii,2) && dataStats(1).latency(ii,jj) > dataStats(2).latency(ii,jj)
            significantStimuli(ii,jj,2) = 1;
        end
    end    
end


if Channel == 1
    chans = {'Target','Off-Target'};
else
    chans = {'Off-Target','Target'};
end

% GRAPHS OF RESPONSES
h = figure();
count = 1;order = [1,2,3,7,8,9,4,5,6,10,11,12];
for ii=1:numChans
    for jj=1:numStimuli
        figure(h);subplot(4,3,order(count));
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
   
        figure(h);subplot(4,3,order(count));
        VEP = squeeze(meanResponse(ii,members,:));
        LAT = squeeze(dataStats(2).latency(ii,members));
        MIN = squeeze(dataStats(2).mean(ii,members));
        SIGNIFS = squeeze(significantStimuli(ii,members,:));
        
        VEP = VEP(indeces,:);
        LAT = LAT(indeces);
        MIN = -MIN(indeces);
        SIGNIFS = SIGNIFS(indeces,:);
        shift = 250;
        for kk=1:length(indeces)
            plot(1:stimLen,VEP(kk,:)+shift*kk,'LineWidth',2);
            hold on;
            if SIGNIFS(kk,1) == 1 
                plot(LAT(kk)*sampleFreq,MIN(kk)+shift*kk+100,'m*');
                hold on;
            end
            if SIGNIFS(kk,2) == 1
                plot(LAT(kk)*sampleFreq,MIN(kk)+shift*kk+50,'c*');
            end
        end
        title(sprintf('%s, Channel: %s, Animal: %d',Stimulus(jj).name,chans{ii},AnimalName));
        xlabel('Time from phase shift (milliseconds)');
        ylabel('Radius (degrees of visual space [log scale])');
        axis([0 stimLen 0 (shift*length(indeces)+shift)]);
        ax = gca;
        ax.YTick = linspace(0,shift*length(indeces)+shift,length(indeces)+2);
        ax.YTickLabel = [2^(log2(sortedDegs(1))-1),sortedDegs,2^(log2(sortedDegs(end))+1)];
        count = count+1;
        
        figure(h);subplot(4,3,order(count));
        
        for kk=1:numStats
            Stats = squeeze(dataStats(kk).latency(ii,members)).*sampleFreq;
            sortedStats = Stats(indeces);
            Errors = squeeze(dataStats(kk).latencySEM(ii,members)).*sampleFreq;
            sortedErrors = Errors(indeces);
            errorbar(log2(sortedDegs),sortedStats,...
                sortedErrors,dataStats(kk).specs,...
                'LineWidth',2);
            hold on;
        end
        title(sprintf('%s, Channel: %s, Animal: %d',Stimulus(jj).name,chans{ii},AnimalName));
        xlabel('Radius (degrees of visual space [log scale])');
        ylabel('Latency to Peak (milliseconds)')
        legend(dataStats(1).name,dataStats(2).name,'Positivity - Negativity');
        axis([(log2(sortedDegs(1))-1) (log2(sortedDegs(end))+1) 0 max(maxWin)+50]);
        ax = gca;
        ax.XTick = [log2(sortedDegs(1))-1,log2(sortedDegs),log2(sortedDegs(end))+1];
        ax.XTickLabel = [2^(log2(sortedDegs(1))-1),sortedDegs,2^(log2(sortedDegs(end))+1)];
        count = count+1;
        
        clear indeces;
    end
end
savefig(h,sprintf('TargetSRP%d_%d.fig',Date,AnimalName));
save(sprintf('TargetSRPConv%d_%d.mat',Date,AnimalName),'Response',...
    'dataStats','baseStats','LatencyStats','baseLatStats','significantStimuli','Stimulus','Radii',...
    'mmPerPixel','DistToScreen','meanResponse','Channel','numChans','numStimuli','numRadii','numStats');

% fileList = dir('TargetSRPConv*.mat');
% Latency = zeros(23,20);
% sizes = zeros(23,20);
% SizeErr = zeros(23,20);
% LatErr = zeros(23,20);
% for ii=1:23
% load(fileList(ii).name);
% Latency(ii,:) = reshape(dataStats(2).latency,[1,20]);
% sizes(ii,:) = reshape(dataStats(2).mean,[1,20]);
% SizeErr(ii,:) = reshape(dataStats(2).stdError,[1,20]);
% LatErr(ii,:) = reshape(dataStats(2).latencySEM,[1,20]);
% end
end

function [Fx] = eCDF(Data,win)
%eCDF.m
%   Creation of the empirical distribution function (empirical cumulative
%   distribution function) for an array of Data
%
%INPUT: Data - the data as a vector
%       OPTIONAL:
%       alpha - confidence level for nonparametric 1-alpha confidence
%         bands, defaults to 0.05 for a 95% confidence band
%OUTPUT: Fx - the empirical distribution function
%        x - the points at which Fx is calculated

%Created: 2016/07/09
%  Byron Price
%Updated: 2016/09/07
%By: Byron Price

n = length(Data);
Fx = zeros(length(win),1);
Data = sort(Data);

count = 1;
for ii=win
    Fx(count) = sum(Data<=ii)/n;
    count = count+1;
end

end




