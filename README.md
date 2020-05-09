# SFND_Radar_Target_Generation_and_Detection

## Writeup for submission
### Radar Specifications
```
Frequency of operation = 77GHz
Max Range = 200m
Range Resolution = 1 m
Max Velocity = 100 m/s
Speed of light = 3e8

Target initial status:
Initial Range = 100 m
Initial Velocity = 10 m/s;
```
### FMCW Waveform Design
```
sweepBandwidth = c/(2*rangeResolution);
Tchirp = 5.5*2*maxRange/c;
sweepSlope = sweepBandwidth/Tchirp;
```

### Simulation Loop
Set the vectors for looping : 
```
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
```

Then run the loop:
```
for i = 1 : length(t)      

    rangeTarget(i) = initTargetRange + initTargetVelocity*t(i);
    timeDelay(i) = 2*rangeTarget(i)/c;
    Tx(i) = cos(2*pi*(carrierFreq*t(i)+sweepSlope*t(i)^2/2));
    Rx(i) = cos(2*pi*(carrierFreq*(t(i)-timeDelay(i)) + (sweepSlope*(t(i)-timeDelay(i))^2)/2));
    Mix(i) = times(Tx(i), Rx(i));
    
end
```
### Range FFT
<p align="center">
  <img  src="https://github.com/paulyehtw/SFND_Radar_Target_Generation_and_Detection/blob/master/results/1D_FFT.png">
</p>

Implement Fast Fourier Transform and do some post-processing: 
```
MixFFT = fft(Mix,Nr)./Nr;
MixFFT = abs(MixFFT);
MixFFT  = MixFFT(1:Nr/2);
```

In the image above, a range of **100m** is detected.
### 2D CFAR
The output from 2D-FFT(Range Doppler Map, RDM) as follows:
<p align="center">
  <img  src="https://github.com/paulyehtw/SFND_Radar_Target_Generation_and_Detection/blob/master/results/2D_FFT.png">
</p>

First set numbers of training and guard cells for both directions, and set offset for threshold:
```
% Select the number of Training Cells in both the dimensions.
Tr = 2;
Td = 1;

% Select the number of Guard Cells in both dimensions around the Cell under 
% test (CUT) for accurate estimation
Gr = 2;
Gd = 1;

% offset the threshold by SNR value in dB
offset = 1.5;
```

Then calculate the number of training cells in a window based on parameters above:
```
windowSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
numTrainingCells = windowSize-(2*Gr+1)*(2*Gd+1);
```

Then slide the window over the Range Doppler Map : 
```
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
```

The output of 2D CA-CFAR : 
<p align="center">
  <img  src="https://github.com/paulyehtw/SFND_Radar_Target_Generation_and_Detection/blob/master/results/2D_CA-CFAR.png">
</p>

As shown in the image, a range of **100m** and a doppler velocity of **10m/s** are detected.
