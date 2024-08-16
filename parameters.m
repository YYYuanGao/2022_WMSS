% BY Ke Jia: Version 1_20220623 15:31
% list all the parameters used in this experiment

Param = struct;

%% Screen Settings
Param.Settings.ViewDistance      = 550;              % 1100 mm 
Param.Settings.ScrnResolution    = [0 0 1600 1200];   % rect_computer = [0 0 40 30];
Param.Settings.SquareSize        = 1600;             % 400 mm 
Param.Settings.SquareLength      = 350;              % 154 mm  
Param.Settings.PixelPerDegree    = 2*Param.Settings.ViewDistance*tan(1/2/180*pi)*Param.Settings.SquareSize/Param.Settings.SquareLength;       

%% Keys for response
Param.Keys.Space    = 32;  
Param.Keys.Esc      = 27;
Param.Keys.Left     = 37;  
Param.Keys.Right    = 39;
Param.Keys.Down     = 40;
Param.Keys.Trigger1 = 83;  % 's'       

%% parameters for stimuli
% Locations
Param.Stimuli.Eccentricity   = 5; % degree
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2 - Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2] + sqrt(2)/2*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree;
Param.Stimuli.Locations(4,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Lower Right','Centre'};

% Gratings
Param.Stimuli.OuterRadius      = 3;
% Param.Stimuli.OuterRadius_edge = Param.Stimuli.OuterRadius - 0.5;
Param.Stimuli.InnerRadius      = 1;
Param.Stimuli.OuterSize        = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
Param.Stimuli.InnerSize        = round(Param.Stimuli.InnerRadius*Param.Settings.PixelPerDegree);   % radius pixels

Param.Stimuli.GratingOri       = (0:15:165)+7.5;
Param.Stimuli.GratingOriDiff   = [0 10 20 30 40 50]; %[15,30,45,60]
Param.Stimuli.OriJitter        = 3; 
Param.Stimuli.ResponseJitter   = 15; 

Param.Stimuli.GratingContrast  = 0.8;
Param.Stimuli.Spatial_freq     = 1/Param.Settings.PixelPerDegree;       % 0.02 % 6 bars in total
% Param.Stimuli.SmoothSD         = Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree/3;          % 0.5degree*1/3
 
%% parameters for fixation
Param.Fixation.CrossColor    = [1,1,1]*255;
Param.Fixation.CrossSize     = 0.3*Param.Settings.PixelPerDegree;
Param.Fixation.CrossWidth    = 0.1*Param.Settings.PixelPerDegree;
% +
Param.Fixation.CrossLoc      = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(3)/2;...
                                     Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize] + [offset',offset',offset',offset'];
% X
Param.Fixation.CrossLoc2     = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2;...
                                     Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2]+[offset',offset',offset',offset'];  

Param.Fixation.OvalSize      = 1*Param.Settings.PixelPerDegree;
Param.Fixation.OvalColor     = [0,0,0]; 
Param.Fixation.Offset        = 0.5*Param.Settings.PixelPerDegree; 
Param.Fixation.FontColor     = [1,1,1]*255;

%% parameters for trials
Param.Trial.ITI              = 0.4;     
Param.Trial.ISI              = 0.4;    % 12 frames= 600 ms
Param.Trial.StiDura          = 0.2; 
Param.Trial.Delay            = 3; 
% Param.Trial.MaxRT            = 1.5; 
% Param.Trial.Duration         = 3;

Param.Trial.TotalNum         = 150; 
Param.Trial.minirun          = 30;
Param.Trial.ResponsePrec     = 1;

%% open window 
screens = Screen('Screens');
screenNumber = max(screens);	
Screen('Preference', 'SkipSyncTests', 0);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);
wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
Screen('BlendFunction',wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('TextFont',wnd,'Arial');
Screen('TextSize',wnd,30);

RefreshDur = Screen('GetFlipInterval',wnd);
RefreshRate = 1./RefreshDur;
Slack = RefreshDur/2;

if abs(RefreshRate-100)>1
    disp(' ');
    disp('Please reset your refreshrate!');
    disp(' ');
    abort;
end
