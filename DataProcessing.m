%% Description:

%       Detailed pipeline of data processing to identify timing and degree
%       of Au NR 3D rotation.
%       (1) Short-time Fourier transform (STFT) using MATLAB built-in fft function
%       (2) Outlier filtering using MATLAB built-in smooth function
%       (3) Change point detection using MATLAB built-in findchangepts function
%       (4) Robust local regression for phi and theta using weighted linear least squares
%       (5) Calculate instantaneous diffusion coefficient
%
% Lead Contact: Daeha Seo <livewire@dgist.ac.kr>
% Credit: Wonhee John Lee <whlee1009@dgist.ac.kr>
% Reference: Siwoo Jin et al., Adv. Sci., 2024

%% Initialization

clear; 
clc;
close all;

% Read x(Column 1), y(Column 2) coordinate of Au NR 
% and its oscillating intensity(Column 3) 
% This values were acquired by tracing the dark-field image sequence
% using the TrackMate plugin of ImageJ/Fiji 
% (Ref: Dmitry Ershov et al., Nature Methods 19., 829-832., 2022) 
load('RawData_Rotating.mat');

% Read oscillating intensity of fixed Au NR for comparison
% This value was measured by Tracking the intensity value using ROI manager in Fiji algorithm  
load('RawData_Fixed.mat');

% Define temporal information
Ts = 1/85.17; % Sampling time interval
Fs = 1/Ts; % Sampling frequency
Time = (0:length(Intensity)-1)'*Ts; % Whole time
PolarizerRotationPeriod = 0.5; % Rotation period of polarizer

%% (1) Short-time Fourier transform (STFT) using MATLAB built-in fft function

FourierWindow = round(PolarizerRotationPeriod/Ts); % Define window length for STFT
Length_STFT = length(Time)-FourierWindow+1; % Length for time dependent angular coordinate

% Time domain for STFT
if mod(FourierWindow,2) == 0
    STFT_Time = Time((1:Length_STFT)+FourierWindow/2-1);
else
    STFT_Time = Time((1:Length_STFT)+(FourierWindow+1)/2-1);
end

% Initialization
STFT_Fixed_Amp = zeros(Length_STFT,1);
STFT_Fixed_Phase = zeros(Length_STFT,1);
STFT_Amp = zeros(Length_STFT,1);
STFT_Phase = zeros(Length_STFT,1);

% For loop for STFT
for i = 1:Length_STFT
    TempIndex = (i:i+FourierWindow-1)'; % Index for short-time masking
    
    % STFT for signal from fixed Au NR
    TempSignal = Intensity_Fix(TempIndex); % Masked signal for one cycle length
    
    % Perform FFT
    X = TempSignal;
    L = length(X);
    Y = fft(X);
    P2 = abs(Y/L);
    if mod(L,2) == 0
        P1 = P2(1:(L/2+1));
        P1(2:end-1) = 2*P1(2:end-1);
    else
        P1 = P2(1:(L+1)/2);
        P1(2:end) = 2*P1(2:end);
    end
    
    % Save amplitude and phase of Fourier component corresponding to k = 1 (index = 2)
    % (i.e., the oscillating signal caused by rotating polarizer)
    STFT_Fixed_Amp(i) = P1(2);
    STFT_Fixed_Phase(i) = angle(Y(2))/2;
    
    % STFT for signal from rotating Au NR
    TempSignal = Intensity(TempIndex);
    
    % Perform FFT
    X = TempSignal;
    L = length(X);
    Y = fft(X);
    P2 = abs(Y/L);
    if mod(L,2) == 0
        P1 = P2(1:(L/2+1));
        P1(2:end-1) = 2*P1(2:end-1);
    else
        P1 = P2(1:(L+1)/2);
        P1(2:end) = 2*P1(2:end);
    end
    
    % Save amplitude and phase of Fourier component corresponding to k = 1 (index = 2)
    % (i.e., the oscillating signal caused by rotating polarizer)
    STFT_Amp(i) = P1(2);
    STFT_Phase(i) = angle(Y(2))/2;    
end

% Unwrapping phase for subtraction
Unwrap_Fixed_Phase = unwrap(2*STFT_Fixed_Phase)/2;
Unwrap_Phase = unwrap(2*STFT_Phase)/2;

% Phi of rotating Au NR
STFT_Phi_Pre = Unwrap_Phase - Unwrap_Fixed_Phase;

% Calculate Min-Max of amplitude
MaxAmp = max(STFT_Amp);
MinAmp = min(STFT_Amp);

% Theta of rotating Au NR
STFT_Theta_Pre = 4*asin((STFT_Amp-MinAmp)./(MaxAmp-MinAmp));
% Multiplied by 4 to match range of theta (0 ~ pi/2) and phi (-pi ~ pi), it will be divided by 4 later(*)
% Range mismatch yields different change point detection sensitivity for theta and phi

%% (2) Outlier filtering using MATLAB built-in smooth function

% Smoothing with rloess option
SW = 21; % Window length for smoothing
STFT_Phi_Smoothing = smooth(STFT_Phi_Pre,SW,'rloess');
STFT_Theta_Smoothing = smooth(STFT_Theta_Pre,SW,'rloess');

%% (3) Change point detection using MATLAB built-in findchangepts function

% Grouping
STFT_Theta_Phi_Smoothing = [STFT_Theta_Smoothing, STFT_Phi_Smoothing];

% Change point detection
MT = 20; % Penalized contrast (MinThreshold)
MD = 50; % Minimum number of samples between changepoints (MinDistance)
ipt = findchangepts(STFT_Theta_Phi_Smoothing','MinThreshold',MT,'MinDistance',MD)';

% Discretized spherical coordinate (theta, phi)
x = STFT_Theta_Phi_Smoothing; % x: Continous data
icp = ipt; % icp: Index of detected change point
istart = [1; icp(:)]; % List of start index of each segment
iend = [icp(:)-1; length(x)]; % List of end index of each segment
nseg = length(icp)+1; % The number of segments
y = zeros(nseg,2); % y: Discretized data
for s = 1:nseg
    ix = (istart(s):iend(s))';
    y(ix,1) = mean(x(ix,1));
    y(ix,2) = mean(x(ix,2));
end

% Divide 4*theta by 4(*)
STFT_Theta_Pre = STFT_Theta_Pre/4;
STFT_Theta_Smoothing = STFT_Theta_Smoothing/4;
STFT_Theta_Phi_Smoothing(:,1) = STFT_Theta_Phi_Smoothing(:,1)/4;
y(:,1) = y(:,1)/4;

%% (4) Robust local regression for phi and theta using weighted linear least squares

% Renormalize the range of theta (0 ~ pi/2) based on the discretized theta
% This process reduces error from spiky noise
MaxState = max(y(:,1)); 
MinState = min(y(:,1));
Final_Theta = (STFT_Theta_Smoothing-MinState)./(MaxState-MinState)*pi/2;
Final_Theta_State = (y(:,1)-MinState)./(MaxState-MinState)*pi/2;

% Reset range of unwrapped phi (-pi ~ pi) for visualization
Final_Phi = zeros(length(STFT_Phi_Smoothing),1);
Final_Phi_State = zeros(length(STFT_Phi_Smoothing),1);
for i = 1:length(STFT_Phi_Smoothing)
    Temp_Phi_State = mod(y(i,2)+pi/2,pi)-pi/2;
    Final_Phi(i) = STFT_Phi_Smoothing(i)-y(i,2)+Temp_Phi_State;
    Final_Phi_State(i) = Temp_Phi_State;
end

%% (5) Calculate instantaneous diffusion coefficient

LengthPerPixel = 0.27; % Length/Pixel (um)
AnomalousFactor = 0.25; %ratio to determine movement within the given range(Ex. Higher AnomalousFactor values lead to considering shorter time intervals)
TimeLag = 4; % lag time(τ) points for linear fitting in calculation of MSD analysis (〈r^2(τ)〉=4Dτ^α) 
D_Cutoff = 0.009; % Cutoff for diffusion coefficient
D_Percent = 0.95; % Ratio of exceeding the D cut off between changepoints

D = zeros(Length_STFT,1);

% Calculation of mean square displacement (MSD) and alpha value
for i = 1:Length_STFT
    TempIndex = (i:i+FourierWindow-1)';
    
    TrackData = Track(TempIndex,1:2)*LengthPerPixel;

    Range = floor((size(TrackData,1)-1)*AnomalousFactor);

    msd = zeros(1,Range); 
    for j = 1:Range
        DeltaCoords = TrackData(1+j:end,1:2) - TrackData(1:end-j,1:2);
        SquaredDisplacement = sum(DeltaCoords.^2,2);
        msd(j) = mean(SquaredDisplacement);
    end

    t = (1:Range)*Ts;

    if TimeLag == 1
        d = (msd(1))/(4*0.001*Ts);
    else
        param_d = polyfit(t(1:TimeLag),msd(1:TimeLag),1);
        d = param_d(1)/4;
    end

    D(i) = d;

    if Range == 1
        alpha = 1;
    else
        param_a = polyfit(log(t),log(msd),1);
        alpha = param_a(1);
    end

    ALPHA(i) = alpha;
end

%marking of change point for 3D rotation
STFT_Full_Theta_Phi = [Final_Theta, Final_Phi];
    icp = ipt; % icp: Index of detected change point
    Newipt = ipt;
 
D_PauseZeroMoveOne = zeros(length(Newipt)-1,1);

%marking of diffusion coefficient in the region between change points using D cut off value,
%and distinguish whether it is in a 'Run' or 'Pause' state.
for i = 1:(length(Newipt)-1)
    TempIndex = Newipt(i):Newipt(i+1);
    TempPercent = zeros(length(TempIndex),1);
    
    cnt = 1;
    for j = TempIndex
       if D(j) >= D_Cutoff
           TempPercent(cnt) = 1;
       end
       cnt = cnt+1;
    end
    
    if sum(TempPercent)/length(TempPercent) >= D_Percent
        D_PauseZeroMoveOne(i) = 1;
    end
end


