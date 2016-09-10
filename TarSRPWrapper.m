function [] = TarSRPWrapper(AnimalName)
% TarSRPWrapper.m
%
%  To be used with the output of TarSRPAnalysis.m ... the file named 
%   TargetSRPConv ... Will take the data from multiple days for review
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%
%OUTPUT: saved file, figures
%
% Created: 2016/09/08, 24 Cummington Mall, Boston, MA
%  Byron Price
% Updated: 2016/09/08
%  By: Byron Price

cd ('~/CloudStation/ByronExp/SRP');
set(0,'DefaultFigureWindowStyle','docked');

% find not unique numbers
% y = sort(vec);
% nonUnique = unique(y(diff(sort(vec))==0));

fileList = dir(sprintf('TargetSRPConv*_%d.mat',AnimalName));

numDates = length(fileList);
Dates = zeros(numDates,1);

VEPsize(1:numDates) = struct;
VEPlat(1:numDates) = struct;
stimDegree(1:numDates) = struct;
LatDist(1:numDates) = struct;

for ii=1:numDates
    index = regexp(fileList(ii).name,'_');
    temp = fileList(ii).name(index-8:index-1);
    Dates(ii) = str2double(temp);
    ConvFileName = sprintf('TargetSRPConv%d_%d.mat',Dates(ii),AnimalName);
    load(ConvFileName);
    
    VEPsize(ii).total = zeros(numChans,numStimuli,numRadii,numStats);
    VEPsize(ii).sem = zeros(numChans,numStimuli,numRadii,numStats);
    VEPlat(ii).total = zeros(numChans,numStimuli,numRadii,numStats);
    VEPlat(ii).sem = zeros(numChans,numStimuli,numRadii,numStats);
    LatDist(ii).total = zeros(numChans,numStimuli,numRadii);
    LatDist(ii).sem = zeros(numChans,numStimuli,numRadii);
    stimDegree(ii).vals = zeros(numRadii,1);
    for jj=1:numChans
        for kk=1:numStimuli
            members = Radii(:,2) == kk;
            cms = (Radii(members,1)).*mmPerPixel./10;
            degrees = atan(cms./DistToScreen).*180/pi;
            degrees = degrees';
            [sortedDegs,indeces] = sort(degrees);
            stimDegree(ii).vals = sortedDegs;
            for ll=1:numStats
                Stats = squeeze(dataStats(ll).mean(jj,members));
                sortedStats = Stats(indeces);
                Errors = squeeze(dataStats(ll).stdError(jj,members));
                sortedErrors = Errors(indeces);
                
                Lats = squeeze(dataStats(ll).latency(jj,members));
                LatErrs = squeeze(dataStats(ll).latencySEM(jj,members));
                sortedLats = Lats(indeces);
                sortedLatErrs = LatErrs(indeces);
                
                VEPsize(ii).total(jj,kk,:,ll) = sortedStats;
                VEPsize(ii).sem(jj,kk,:,ll) = sortedErrors;
                VEPlat(ii).total(jj,kk,:,ll) = sortedLats;
                VEPlat(ii).sem(jj,kk,:,ll) = sortedLatErrs';
            end
            Lats = squeeze(LatencyStats.mean(jj,members));
            LatErrs = squeeze(LatencyStats.sem(jj,members));
            
            Lats = Lats(indeces);
            LatErrs = LatErrs(indeces);
            
            LatDist(ii).total(jj,kk,:) = Lats;
            LatDist(ii).sem(jj,kk,:) = LatErrs;
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
order = [1,2,3,7,8,9,4,5,6,10,11,12];

for ii=1:length(order)
    ax(ii) = subplot(4,3,order(ii));
    hold(ax(ii),'on');
end

count = 1;
for jj=1:numChans
    for kk=1:numStimuli
        vepsizes = zeros(numDates,numRadii);
        vepsem = zeros(numDates,numRadii);
        for mm=1:numDates
            vepsizes(mm,:) = squeeze(VEPsize(mm).total(jj,kk,:,3));
            vepsem(mm,:) = squeeze(VEPsize(mm).sem(jj,kk,:,3));
        end
        for ll=1:numRadii
            errorbar(ax(order(count)),1:numDates,vepsizes(:,ll),...
                vepsem(:,ll),...
                'LineWidth',2);
        end
        title(ax(order(count)),sprintf('%s, Channel: %s, Animal: %d',Stimulus(kk).name,chans{jj},AnimalName));
        xlabel(ax(order(count)),'Experimental Day');
        legend(ax(order(count)),'4','8','16','32','64');
        ylabel(ax(order(count)),'Total VEP Magnitude(\muVolts)')
        axis(ax(order(count)),[0 numDates+1 0 500]);
        
        count = count+1;
        
        veplats = zeros(numDates,numRadii);
        veplatsem = zeros(numDates,numRadii);
        for mm=1:numDates
            veplats(mm,:) = squeeze(VEPlat(mm).total(jj,kk,:,2));
            veplatsem(mm,:) = squeeze(VEPlat(mm).sem(jj,kk,:,2));
        end
        for ll=1:numRadii
            errorbar(ax(order(count)),1:numDates,veplats(:,ll).*1000,...
                veplatsem(:,ll).*1000,...
                'LineWidth',2);
        end
        title(ax(order(count)),sprintf('%s, Channel: %s, Animal: %d',Stimulus(kk).name,chans{jj},AnimalName));
        xlabel(ax(order(count)),'Experimental Day');
        legend(ax(order(count)),'4','8','16','32','64');
        ylabel(ax(order(count)),'Latency to Peak Negativity (milliseconds)')
        axis(ax(order(count)),[0 numDates+1 0 200]);
        count = count+1;
        
        latdist = zeros(numDates,numRadii);
        latdistsem = zeros(numDates,numRadii);
        for mm=1:numDates
           latdist(mm,:) = squeeze(LatDist(mm).total(jj,kk,:));
           latdistsem(mm,:) = squeeze(LatDist(mm).sem(jj,kk,:));
        end
        for ll=1:numRadii
            errorbar(ax(order(count)),1:numDates,latdist(:,ll),...
                latdistsem(:,ll),...
                'LineWidth',2);
        end
        title(ax(order(count)),sprintf('%s, Channel: %s, Animal: %d',Stimulus(kk).name,chans{jj},AnimalName));
        xlabel(ax(order(count)),'Experimental Day');
        legend(ax(order(count)),'4','8','16','32','64');
        ylabel(ax(order(count)),'KS Stat for Peak Negativity')
        axis(ax(order(count)),[0 numDates+1 0 0.5]);
        count = count+1;
    end
end

savefig(h,sprintf('TargetSRPDays_%d.fig',AnimalName));


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


