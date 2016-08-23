function [] = TargetSRP(AnimalName,holdTime)
%TargetSRP.m
%  Will display a series of  small, circular phase-reversing sinusoidal 
%   gratings with an alpha mask that decays over a varying range of radii. 
%   For each radius, there will be an ON center stimulus (you see the 
%   grating within a circle of that radius and grey beyond) and an OFF center 
%   stimulus (you see the grating beyond the circle and grey within). Must
%   have previously used the Retinotopy.m and MapRetinotopy.m files to
%   determine the retinotopic location of the LFP recording electrode being
%   used.
%
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%
%        Optional- 
%        holdTime - amount of time to wait between blocks of about 50 stimuli
% 
%        see file SRPVars.mat for other changeable presets
%
% OUTPUT: a file with stimulus parameters named TargetSRPDate_AnimalName
%           e.g. TargetSRP20160708_12345.mat , to be saved on CloudStation
%           under '~/CloudStation/ByronExp/SRP'
%
% Created: 2016/08/11 at 24 Cummington, Boston, MA
%  Byron Price
% Updated: 2016/08/18
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');
load(sprintf('RetinoMap%d.mat',AnimalName));

centerMass = MapParams.centerMass;
Channel = MapParams.Channel;

cd('~/CloudStation/ByronExp/SRP');
load('SRPvars.mat');

directory = '~/Documents/MATLAB/Byron/Targeted-SRP';

if nargin < 2
    holdTime = 30;
end

reps = reps-mod(reps,blocks);
orientation = orientation*pi/180;
phaseShift = phaseShift*pi/180;

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

usb = usb1208FSPlusClass;
display(usb);
WaitSecs(10);

% % Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));
% 
% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black;255 = white
background = 127;
[win,~] = Screen('OpenWindow', screenid,background);

gammaTable = makeGrayscaleGammaTable(gama,0,255);
Screen('LoadNormalizedGammaTable',win,gammaTable);

% Switch color specification to use the 0.0 - 1.0 range
Screen('ColorRange', win, 1);

% % Query window size in pixels
[w_pixels, h_pixels] = Screen('WindowSize', win);
% 
% % Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

dgshader = [directory '/TargetSRP.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/TargetSRP.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;
mmPerPixel = conv_factor;
conv_factor = 1/conv_factor;

degreeRadii = zeros(numStimuli,numRadii);
for ii=1:numRadii
    degreeRadii(:,ii) = 2^(ii);
end
% perform unit conversions
Radii = (tan(degreeRadii.*pi./180).*(DistToScreen*10)).*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy

temp = (tan((1/spatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;clear temp;

centerVals = zeros(2,1);
centerVals(1) = centerMass.x(Channel);centerVals(2) = centerMass.y(Channel);

alpha = ones(numStimuli,numRadii);

% 0 for ON center, 1 for OFF center
onOFF = [0,1];

estimatedTime = ((stimTime*reps/blocks+holdTime)*blocks)*numRadii*numStimuli/60;
display(sprintf('\nEstimated time: %3.2f minutes',estimatedTime));

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

% Perform initial flip to gray background and sync us to the retrace:
Priority(9);

% we have event words from 0 to 255, so 3 stimuli and greater than 10
% radii would break the code
stimNums = zeros(numStimuli,numRadii);
for ii=1:numStimuli
    for jj=1:numRadii
    mystr = sprintf('%d%d',ii,jj);
    stimNums(ii,jj) = str2double(mystr);
    end
end

usb.startRecording;
WaitSecs(1);
usb.strobeEventWord(0);
WaitSecs(holdTime);

% Animation loop
vbl = Screen('Flip',win);
for ii=1:numStimuli
    for jj=1:numRadii
        for kk=1:blocks
%             vbl = Screen('Flip',win,vbl+ifi/2);
            phase = 0;
            for ll = 1:reps/blocks
                phase = phase+phaseShift;
                % Draw the procedural texture as any other texture via 'DrawTexture'
                Screen('DrawTexture', win,gratingTex,[],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[alpha(ii,jj),phase,...
                    Radii(ii,jj),centerVals(1),centerVals(2),spatFreq,orientation,onOFF(ii)]);
                % Request stimulus onset
                vbl = Screen('Flip', win,vbl-ifi/2+stimTime);
                usb.strobeEventWord(stimNums(ii,jj));
            end
            vbl = Screen('Flip',win,vbl+ifi/2);usb.strobeEventWord(0);
            vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
        end
    end
end
WaitSecs(2);
usb.stopRecording;
Priority(0);

Stimulus = struct('name',cell(numStimuli,1));
Stimulus(1).name = 'Gabor patch';
Stimulus(2).name = 'Full-Field Grating with Grey Center';

cd('~/CloudStation/ByronExp/SRP');
fileName = sprintf('TargetSRP%d_%d.mat',Date,AnimalName);
save(fileName,'centerVals','Radii','reps','stimTime','numStimuli',...
    'w_pixels','h_pixels','spatFreq','mmPerPixel','holdTime','Channel',...
    'DistToScreen','orientation','phaseShift','numRadii','Stimulus')
% Close window
Screen('CloseAll');
end

function gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
% Generates a 256x3 gamma lookup table suitable for use with the
% psychtoolbox Screen('LoadNormalizedGammaTable',win,gammaTable) command
% 
% gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
%
%   gamma defines the level of gamma correction (1.8 or 2.2 common)
%   blackSetPoint should be the highest value that results in a non-unique
%   luminance value on the monitor being used (sometimes values 0,1,2, all
%   produce the same black pixel value; set to zero if this is not a
%   concern)
%   whiteSetPoint should be the lowest value that returns a non-unique
%   luminance value (deal with any saturation at the high end)
% 
%   Both black and white set points should be defined on a 0:255 scale

gamma = max([gamma 1e-4]); % handle zero gamma case
gammaVals = linspace(blackSetPoint/255,whiteSetPoint/255,256).^(1./gamma);
gammaTable = repmat(gammaVals(:),1,3);
end
