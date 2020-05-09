clear all;
clc;
close all;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
% Speed of light = 3e8
%%%%%%%%%%%%%%%%%%%%%%%%%%%
carrierFreq = 77e9;
maxRange = 200;
rangeResolution = 1;
maxV = 70;
minV = -70;
velResolution = 3;
c = 3e8;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
initTargetRange = 100;
initTargetVelocity = 10;


%% FMCW Waveform Generation

% *%TODO* :
% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

sweepBandwidth = c/(2*rangeResolution);
Tchirp = 5.5*2*maxRange/c;
sweepSlope = sweepBandwidth/Tchirp;

                                                          
%% Set up loop sequence                              
% The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
% for Doppler Estimation. 
Nd = 128;                   % #of doppler cells OR #of sent periods % number of chirps

% The number of samples on each chirp. 
Nr = 1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t = linspace(0, Nd*Tchirp, Nr*Nd); %total time for samples


% Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t)); %transmitted signal
Rx = zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

% Similar vectors for range_covered and time delay.
rangeTarget = zeros(1,length(t));
timeDelay = zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i = 1 : length(t)         
 
    % *%TODO* :
    % For each time stamp update the Range of the Target for constant velocity.
    rangeTarget(i) = initTargetRange + initTargetVelocity*t(i);
    timeDelay(i) = 2*rangeTarget(i)/c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(carrierFreq*t(i)+sweepSlope*t(i)^2/2));
    Rx(i) = cos(2*pi*(carrierFreq*(t(i)-timeDelay(i)) + (sweepSlope*(t(i)-timeDelay(i))^2)/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = times(Tx(i), Rx(i));
    
end

%% RANGE MEASUREMENT
% *%TODO* :
% run the FFT on the beat signal along the range bins dimension (Nr) and
% normalize.
MixFFT = fft(Mix,Nr)./Nr;

% *%TODO* :
% Take the absolute value of FFT output
MixFFT = abs(MixFFT);

% *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
MixFFT  = MixFFT(1:Nr/2);

% plotting the range
figure ('Name', 'Range from First FFT');

% *%TODO* :
% plot FFT output 
plot(MixFFT);
pbaspect([2 1 1]);
xlabel('Range');
ylabel('dB');
axis ([0 maxRange 0 1]);

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

% use the surf function to plot the output of 2DFFT and to show axis in both
% dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
xlabel('Doppler Velocity');
ylabel('Range');

%% CFAR implementation

% Slide Window through the complete Range Doppler Map

% *%TODO* :
% Select the number of Training Cells in both the dimensions.
Tr = 2;
Td = 1;

% *%TODO* :
% Select the number of Guard Cells in both dimensions around the Cell under 
% test (CUT) for accurate estimation
Gr = 2;
Gd = 1;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 1.5;

% *%TODO* :
% Create a vector to store noise_level for each iteration on training cells
windowSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
numTrainingCells = windowSize-(2*Gr+1)*(2*Gd+1);
mapCFAR = zeros(size(RDM));
[RDMr, RDMd] = size(RDM);

% *%TODO* :
% design a loop such that it slides the CUT across range doppler map by
% giving margins at the edges for Training and Guard Cells.
% For every iteration sum the signal level within all the training
% cells. To sum convert the value from logarithmic to linear using db2pow
% function. Average the summed values for all of the training
% cells used. After averaging convert it back to logarithimic using pow2db.
% Further add the offset to it to determine the threshold. Next, compare the
% signal under CUT with this threshold. If the CUT level > threshold assign
% it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR
for i = 1 : RDMd-2*(Td+Gd)
    for j = 1 : RDMr-2*(Tr+Gr)
        
        % Crop the RDM with patch size
        trainingCellsPatch = db2pow( RDM(j : j+2*(Tr+Gr), i : i+2*(Td+Gd)) );
        % set guard cells and CUT to 0
        trainingCellsPatch(Tr+1 : end-Tr, Td+1 : end-Td) = 0;
        
        noiseLevel = pow2db(sum(sum(trainingCellsPatch))/numTrainingCells);
        noiseLevel = noiseLevel*offset; % scale noise level with offset
        
        if RDM(j+(Td+Gr),i+(Td+Gd)) > noiseLevel
            mapCFAR(j+(Td+Gr), i+(Td+Gd)) = 1;
        else
            mapCFAR(j+(Td+Gr), i+(Td+Gd)) = 0;
        end
           
    end
end

% *%TODO* :
% display the CFAR output using the Surf function like we did for Range
% Doppler Response output.
figure('Name','CA-CFAR Filtered RDM'),surf(doppler_axis, range_axis, mapCFAR);
xlabel('Doppler Velocity');
ylabel('Range');
colorbar;
