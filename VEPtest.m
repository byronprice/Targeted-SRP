function [] = VEPtest(AnimalName,Date)
% VEPtest.m
%
%  Show procedure for determination of whether or not a VEP occurred.
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%
%OUTPUT: saved file, figures
%
% Created: 2016/08/15, 24 Cummington Mall, Boston, MA
%  Byron Price
% Updated: 2016/09/07
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
stimLen = round(0.3*sampleFreq); % 300 milliseconds
minWin = round(0.05*sampleFreq):round(0.2*sampleFreq);
maxWin = round(.1*sampleFreq):1:round(0.3*sampleFreq);
smoothKernel = 4;
alpha = 0.05;

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

for ii=1:numStimuli*numRadii
    h(ii) = figure;
end

LatencyStats = struct;
[~,temp] = min(Response(:,:,:,minWin),[],4);
LatencyStats.trials = (temp+minWin(1)-1)./sampleFreq;
LatencyStats.mad = zeros(numChans,numStimuli*numRadii);
LatencyStats.sem = zeros(numChans,numStimuli*numRadii);
LatencyStats.ci = zeros(numChans,numStimuli*numRadii,2);
LatFun = @(x,win) max(eCDF(x.*sampleFreq,win)-linspace(0,1,length(win))'); %mad(x,1);

for ii=1:numChans
    for stim = 1:numStimuli*numRadii
        count = 1+4*(ii-1);
        figure(h(stim));subplot(numChans,4,count);
        hold on;
        for jj=1:reps
            plot(squeeze(Response(ii,stim,jj,:)));
        end
        plot(mean(squeeze(Response(ii,stim,:,:)),1),'LineWidth',3)
        axis([0 stimLen -1000 1000]);
        hold off;
        count = count+1;
        subplot(numChans,4,count);
        histogram(squeeze(LatencyStats.trials(ii,stim,:)));
        axis([0 0.5 0 100]);
        x = squeeze(LatencyStats.trials(ii,stim,:));
        width = LatFun(squeeze(LatencyStats.trials(ii,stim,:)),minWin);
        legend(sprintf('KS stat = %3.3f',width))
        count = count+1;
        
        subplot(numChans,4,count);
        plot(eCDF(x.*sampleFreq,minWin));hold on;plot(linspace(0,1,length(minWin)));
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
[dataStats(1).mean,maxLats] = max(meanResponse(:,:,maxWin),[],3);
dataStats(1).latency = (maxLats+maxWin(1)-1)./sampleFreq;
[mins,minLats] = min(meanResponse(:,:,minWin),[],3);
dataStats(2).latency = (minLats+minWin(1)-1)./sampleFreq;
dataStats(2).mean = -mins;
dataStats(3).mean = max(meanResponse(:,:,maxWin),[],3)-min(meanResponse(:,:,minWin),[],3);

% BOOTSTRAP FOR STANDARD ERROR OF STATISTICS IN PRESENCE OF VISUAL STIMULI
N = 2000;
for ii=1:numChans
    for jj=1:numStimuli*numRadii
        count = 4+4*(ii-1);
        Tboot = zeros(N,numStats);
        Latency = zeros(N,numStats);
        
        LatBoot = zeros(N,1);
        LatencyStats.mad(ii,jj) = LatFun(squeeze(LatencyStats.trials(ii,jj,:)),minWin);
        for mm=1:N
            indeces = random('Discrete Uniform',reps,[reps,1]);
            LatGroup = squeeze(LatencyStats.trials(ii,jj,indeces));
            LatBoot(mm) = LatFun(LatGroup,minWin);
            group = squeeze(Response(ii,jj,indeces,:));
            meanGroup = mean(group,1);
            [Tboot(mm,1),maxLats] = max(meanGroup(maxWin));
            Latency(mm,1) = (maxLats+maxWin(1)-1)./sampleFreq;
            [mins,minLats] = min(meanGroup(minWin));
            Latency(mm,2) = (minLats+minWin(1)-1)./sampleFreq;
            Tboot(mm,2) = -mins;
            Tboot(mm,3) = max(meanGroup(maxWin))-min(meanGroup(minWin));
        end
        figure(h(jj));subplot(numChans,4,count);histogram(LatBoot);
        axis([0 1 0 500]);
        LatencyStats.sem(ii,jj) = std(LatBoot);
        LatencyStats.ci(ii,jj,:) = [quantile(LatBoot,alpha/2),quantile(LatBoot,1-alpha/2)];
        legend(sprintf('95%% CI: [%3.4f,%3.4f]',LatencyStats.ci(ii,jj,1),LatencyStats.ci(ii,jj,2)))

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
% noStimLen = holdTime*sampleFreq-stimLen*2;
% 
% baseLatStats = struct;
% baseLatStats.mad = zeros(numChans,1);
% baseLatStats.sem = zeros(numChans,1);
% baseLatStats.ci = zeros(numChans,2);
% 
% baseStats = struct;
% for ii=1:numStats
%     baseStats(ii).ci = zeros(numChans,2);
%     baseStats(ii).mean = zeros(numChans,1);
%     baseStats(ii).stdError = zeros(numChans,1);
%     baseStats(ii).latency = zeros(numChans,1);
%     baseStats(ii).latencySEM = zeros(numChans,1);
% end
% 
% pauseOnset = strobeTimes(svStrobed == 0);nums = length(pauseOnset);
% for ii=1:numChans
%     Tboot = zeros(N,numStats);
%     Latency = zeros(N,numStats);
%     LatBoot = zeros(N,1);
%     
%     count = numStimuli*numRadii*4+1;
%     for jj=1:N
%         indeces = random('Discrete Uniform',noStimLen,[reps,1]);
%         temp = zeros(reps,stimLen);
%         num = random('Discrete Uniform',nums);
%         [~,index] = min(abs(timeStamps-pauseOnset(num)));
%         index = index+stimLen;
%         for kk=1:reps
%             temp(kk,:) = ChanData(index+indeces(kk):index+indeces(kk)+stimLen-1,ii);
%         end
%         if jj==100
%             figure(h(ii));subplot(numRows,4,count);
%             hold on;
%             for kk=1:reps
%                 plot(squeeze(temp(kk,:)));
%             end
%             plot(mean(temp,1),'LineWidth',3);
%             hold off;
%             count = count+1;
%         end
%         [~,minLats] = min(temp(:,minWin),[],2);
%         minLats = (minLats+minWin(1)-1)./sampleFreq;
%         if jj==100
%             figure(h(ii));subplot(numRows,4,count);
%             histogram(minLats);axis([0 .5 0 100])
%             width = LatFun(minLats,minWin);
%             legend(sprintf('KS stat = %3.3f',width))
%             count = count+1;
%             subplot(numRows,4,count);
%             plot(eCDF(minLats.*sampleFreq,minWin));hold on;plot(linspace(0,1,length(minWin)));
%             count = count+1;
%         end
%         LatBoot(jj) = LatFun(minLats,minWin);
%         meanTrace = mean(temp,1);
%         [Tboot(jj,1),maxLats] = max(meanTrace(maxWin));
%         Latency(jj,1) = (maxLats+maxWin(1)-1)./sampleFreq;
%         [mins,minLats] = min(meanTrace(minWin));
%         Latency(jj,2) = (minLats+minWin(1)-1)./sampleFreq;
%         Tboot(jj,2) = -mins;
%         Tboot(jj,3) = max(meanTrace(maxWin))-min(meanTrace(minWin));
%     end
%     figure(h(ii));subplot(numRows,4,count);histogram(LatBoot);
%     axis([0 1 0 500])
%     baseLatStats.mad(ii) = mean(LatBoot);
%     baseLatStats.sem(ii) = std(LatBoot);
%     baseLatStats.ci(ii,:) = [quantile(LatBoot,alpha/2),quantile(LatBoot,1-alpha/2)];
%     legend(sprintf('95%% CI: [%3.4f,%3.4f]',baseLatStats.ci(ii,1),baseLatStats.ci(ii,2)))
%     for ll=1:numStats
%         baseStats(ll).ci(ii,:) = [quantile(Tboot(:,ll),alpha/2),quantile(Tboot(:,ll),1-alpha/2)];
%         baseStats(ll).mean(ii) = mean(Tboot(:,ll));
%         baseStats(ll).stdError(ii) = std(Tboot(:,ll));
%         baseStats(ll).latency(ii) = mean(Latency(:,ll));
%         baseStats(ll).latencySEM(ii) = std(Latency(:,ll));
%     end
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



