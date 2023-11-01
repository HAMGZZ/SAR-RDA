% Lewis Hamilton 28/10/2023
% With help from Vanishka Kapoor and Michael Jae Arbotante
%
% An RDA implementation for point target simulation and RADARSAT echo data processing

clc; clear all; close all;

% Simulation? 1 for sim, 0 for RADARSAT
simulation = 0;

% Fast? Fase for black and white images, slow for nice colour ones
fast = 1;

% Load the 'echo.mat' data
if simulation == 0
    load('echo.mat');
end

% How many sim point targets?
NumberofSimTargets = 3;

% Constants
C = 2.9979e8; % Speed of Light

R0  = 0 ; % Center Slant Range U: m
Vr  = 0 ; % Radar Velocity U: m/s
Tr  = 0 ; % Pulse Duration U: s
Kr  = 0 ; % Pulse Rate U: Hz/s
f0  = 0 ; % Carrier (radar) Frequency U: Hz
Fr  = 0 ; % Smapling Rate U: Hz
Fa  = 0 ; % Pulse Repetition Frequecny U: Hz
Naz = 0 ; % Range lines
Nrg = 0 ; % Smaples per Range
fc  = 0 ; % Doppler centroid U: Hz

% Radar Params for RADARSAT: 
if simulation == 0
    R0  = 988647.462; % Center Slant Range U: m
    Vr  = 7062;       % Radar Velocity U: m/s
    Tr  = 41.74e-6;   % Pulse Duration U: s
    Kr  = -0.72135e12;% Pulse Rate U: Hz/s
    f0  = 5.3e9;      % Carrier (radar) Frequency U: Hz
    Fr  = 32.317e6;   % Smapling Rate U: Hz
    Fa  = 1256.98;    % Pulse Repetition Frequecny U: Hz
    Naz = 4096;       % Range lines
    Nrg = 4096;       % Smaples per Range
    fc  = -596.271;   % Doppler centroid U: Hz
end

% Radar Params for Simulation: 
if simulation == 1
    R0  = 20e3;       % Center Slant Range U: m
    Vr  = 150;        % Radar Velocity U: m/s
    Tr  = 2.5e-6;     % Pulse Duration U: s
    Kr  = 20e12;      % Pulse Rate U: Hz/s
    f0  = 5.3e9;      % Carrier (radar) Frequency U: Hz
    Fr  = 60e6;       % Smapling Rate U: Hz
    Fa  = 100;        % Pulse Repetition Frequecny U: Hz
    Naz = 256;        % Range lines
    Nrg = 320;        % Smaples per Range
    fc  = 0;          % Doppler centroid U: Hz
end

% Calculated Values
lambda = C / f0; % Wavelength
t0 = 2 * R0 / C;  % Starting time of the data window

% Freq and Time Axis
tr = t0 + (Nrg-1) / Fr;
Frg = ((0:Nrg-1) - Nrg/2) / Nrg * Fr;
ta = ((0:Naz-1) - Naz/2) * 1/Fa;
Faz = fc + ((0:Naz-1) - Naz/2) / Naz * Fa;

% Var to hold echo data, be it sim or not.
data = 0;

if simulation == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate random target positions
    for i = 1:NumberofSimTargets
        xmin = 0;
        xmax = Nrg;
        ymin = -Naz;
        ymax = Naz;
        x =  R0 + (xmax-xmin).*rand() + xmin;
        y = (ymax-ymin).*rand() + ymax;
        xs(i) = x;
        ys(i) = y;
        btc = (y-x)/Vr;
        btcs(i) = btc;
    end

    % if it should be made more accurate in terms of physical antenna size
    La_real = 8; % Antenna length in azimuth
    bBandwidth = 0.886 * lambda / La_real; % Radar 3dB beamwidth

    % simulation freq and time Axis
    s_tr = 2 * x(1) / C + (-Nrg / 2 : (Nrg / 2 - 1)) / Fr;
    s_ta = (-Naz / 2 : Naz / 2 - 1) / Fa;

    % Generate azimuth and frequency matrix.
    trAxis = ones(Naz,1) * s_tr;
    taAxis = s_ta.' * ones(1,Nrg);

    % Generate simulated data for each target
    s_k = zeros(Naz, Nrg, NumberofSimTargets);
    for k = 1:NumberofSimTargets
        Rn = sqrt((xs(k) .* ones(Naz,Nrg)) .^ 2 + (Vr .* taAxis - ys(k) .* ones(Naz,Nrg)) .^ 2);
        inSwath = ((abs(trAxis - 2 .* Rn ./ C)) <= ((Tr / 2) .* ones(Naz,Nrg)));

        % sinc-squared function, Formula 4.31 in textbook
        s = atan(Vr .* (taAxis - btcs(k) .* ones(Naz,Nrg)) / xs(k)); 
        azimuth = (sinc(0.886 .* s ./ bBandwidth)) .^ 2;
    
        % reflection of LFMW, Formula 4.32 in textbook
        s_k(:,:,k) = inSwath .* azimuth .* exp(-(1j * 4 * pi * f0) .* Rn ./ C) .* exp((1j * pi * Kr) .* (trAxis - 2 .* Rn ./ C) .^ 2);
    end

    % Load the simlatued point target data into the data matrix for RDA
    data = sum(s_k, 3);

    % Display each simulated signal and the sum.
    plotSimSignals(s_k, NumberofSimTargets);
else
    data = double(echo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Raw Data
plotData(data, fast, 'Raw');

% Align data to doppler centroid
disp("Allign to doppler")
data = data .* exp(-1j * 2 * pi * fc * (ta' * ones(1, Nrg)));
disp("Allign to doppler done")

% Range compression
disp("Range Compression")
GFilter = exp(1j * pi * ((Frg.^2) / Kr));
data = ifty(fty(data) .* GFilter);
disp("Range Compression done")

% Display Range Compressed
plotData(data, fast, 'Range Compressed');

% Azimuth FFT
disp("Azimuth FFT")
data = ftx(data);
disp("Azimuth FFT done")

% Display Azimuth FFT
plotData(data, fast, 'Azimuth FFT');

% RCMC
disp("RCMC")
dR = lambda^2 * R0 .* Faz.^2 / (8 * Vr^2);
[Frg_2D, dR2D] = meshgrid(Frg, dR);
G = exp(1j * 4 * pi * Frg_2D .* dR2D / C);
data = data .* G;
disp("RCMC done")

% Display RCMC
plotData(data, fast, 'RCMC');

% Azimuth compression
disp("Azimuth compression")
Ka = 2 * Vr^2 / lambda / R0;
H = exp(-1j * pi * Faz.^2 ./ Ka);
H_2D = repmat(H.', [1, size(data, 2)]);  % Replicate H across the range dimension
data = data .* H_2D;
disp("Azimuth compression done")

% Display H filter
figure; 
plot(Faz,H);
title('H Filter');
shg

plotData(H_2D, fast, 'H Filter Matrix');

% Display after Azimuth Compression
plotData(data, fast, 'Azimuth Compression');

% Inverse FFT to get the image
disp("Final FFT")
data = iftx(data);
data = abs(data);
disp("Final FFT done")

% Normalize the SAR image with the maximum pixel amplitude
SARImage = data / max(data(:));

% Display the final SAR image
if fast == 0
    figure; 
    pcolor(SARImage); 
    shading flat
    title('Final Raw (unedited) SAR Image');
    h = colorbar;
    ylabel(h, 'Normalized Amplitude')
    colormap 'jet'
    shg
else  
    figure; imshow(abs(SARImage), []); title('Final Raw (unedited) SAR Image');
end

% Image processing
if simulation == 0
    % Constrast increase
    disp("Image contrast increase")
    SARImageEdit = imadjust(SARImage, [0.00015 0.06]);
    disp("Image contrast increase done")
    if fast == 0
        figure; 
        pcolor(SARImageEdit); 
        shading flat
        title('Constrast Increased');
        h = colorbar;
        ylabel(h, 'Normalized Amplitude')
        colormap 'jet'
        shg
    else  
        figure; imshow(abs(SARImageEdit), []); title('Constrast Increased');
    end

    % Circshift
    disp("image shift")
    SARImageShifted = circshift(SARImageEdit, -295, 1);
    SARImageShifted = circshift(SARImageShifted, -580, 2);
    disp("image shift done")
    if fast == 0
        figure; 
        pcolor(SARImageShifted); 
        shading flat
        title('Image shifted');
        h = colorbar;
        ylabel(h, 'Normalized Amplitude')
        colormap 'jet'
        shg
    else  
        figure; imshow(abs(SARImageShifted), []); title('Image shifted');
    end

    % Despekle
    disp("image despeckle")
    SARImageSpek = specklefilt(SARImageShifted,DegreeOfSmoothing=0.2,NumIterations=50);
    disp("image despeckle done")
    if fast == 0
        figure; 
        pcolor(SARImageSpek); 
        shading flat
        title('Speckle Filter');
        h = colorbar;
        ylabel(h, 'Normalized Amplitude')
        colormap 'jet'
        shg
    else  
        figure; imshow(abs(SARImageSpek), []); title('Speckle Filter');
    end
end

function plotSimSignals(d, NumberofSims)
    number = NumberofSims;

    %Row Major Order
    r = ceil(sqrt(number));
    c = ceil(number / r);
    figure;
    for i = 1:number
        data = d(:,:,i);
        data = squeeze(data);
        subplot(r,c,i)
        pcolor(abs(data));
        shading flat
        title(string(append('A s: ', int2str(i))));
        xlabel('Range');
        ylabel('Azimuth');
        % No colour bar
        colormap 'jet'
    end
    shg
end

function plotData(data, speed, t)
    if speed == 0
        figure; 
        subplot(2,2,1)
        pcolor(real(data)); 
        shading flat
        title(string(append('Real ', t, ' Data')));
        xlabel('Range');
        ylabel('Azimuth');
        h = colorbar;
        ylabel(h, 'Normalized Amplitude')
        colormap 'jet'
    
        subplot(2,2,2)
        pcolor(imag(data)); 
        shading flat
        title(append('Imaginary ', t, ' Data'));
        xlabel('Range');
        ylabel('Azimuth');
        h = colorbar;
        ylabel(h, 'Normalized Amplitude')
        colormap 'jet'
    
        subplot(2,2,3)
        pcolor(abs(data)); 
        shading flat
        title(append('Amplitude of ', t, ' Data'));
        xlabel('Range');
        ylabel('Azimuth');
        h = colorbar;
        ylabel(h, 'Normalized Amplitude')
        colormap 'jet'
    
        subplot(2,2,4)
        pcolor(angle(data)); 
        shading flat
        title(append('Phase of ', t, ' Data'));
        xlabel('Range');
        ylabel('Azimuth');
        h = colorbar;
        ylabel(h, 'Normalized Amplitude')
        colormap 'jet'
    
        shg
    else  
        figure; imshow(abs(data), []); title(string(append(t, ' Data (Amplitude)')));
    end
end

% Functions used
function fs=ftx(s)
    fs=fftshift(fft(fftshift(s)));
end

function fs=fty(s)
    fs=fftshift(fft(fftshift(s.'))).';
end
   
function s=iftx(fs)
    s=fftshift(ifft(fftshift(fs)));
end

function fs=ifty(fs)
    fs=fftshift(ifft(fftshift(fs.'))).';
end
