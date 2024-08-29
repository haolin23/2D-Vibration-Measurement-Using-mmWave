function [T1_Rx_chirp,T1_Rx_FFT,T2_Rx_chirp,T2_Rx_FFT,T3_Rx_chirp,T3_Rx_FFT] = get_chirp_data(filename,cut_range_bin)
%
format long;

numADCSamples = 256; % number of ADC samples per chirp
numADCBits = 16; % number of ADC bits per sample
numRX = 4; % number of receivers
numLanes = 2; % do not change. number of lanes is always 2
isReal = 0; % set to 1 if real only data, 0 if complex data0
numTx = 3;

%%% read file
fid = fopen(filename);
adcData_raw = fread(fid, 'int16');

%%%if 12 or 14 bits ADC per sample compensate for sign extension
if numADCBits ~= 16
  l_max = 2^(numADCBits - 1) - 1;
  adcData_raw(adcData_raw > l_max) = adcData_raw(adcData_raw > l_max) - 2^numADCBits;
end
fclose(fid);
fileSize = size(adcData_raw, 1);  %%% total number of samples

%%%real data reshape, filesize = numADCSamples * numChirps
if isReal
    numChirps = fileSize/numADCSamples/numRX;
    LVDS = zeros(1, fileSize);
    %create column for each chirp
    LVDS = reshape(adcData_raw, numADCSamples*numRX, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
else
    % for complex data
    numChirps = fileSize/2/numADCSamples/numRX; %% total number of chirps
    
    %%if numChirps is not integer times of numTx, cut the end samples 
    if mod(numChirps, numTx) ~= 0
      numChirps = numChirps - mod(numChirps, numTx);
      fileSize = 2 * numADCSamples * numChirps * numRX;
      adcData = adcData_raw(1:fileSize);
    else
      adcData = adcData_raw;
    end   
    
    %%%%% concantenate samples into complex numbers
    LVDS = zeros(1, fileSize/2);
    counter = 1;
    %disp(fileSize);
    for i=1:4:fileSize-1
      LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2);
      LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3);
      counter = counter + 2;
    end
    %disp(length(LVDS));
end
%%% row is total samples of 4 Rx for each chirp  
LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
LVDS = LVDS.'; %%% transpose

%seperate adc_data per Tx
adcData1 = zeros(numRX,numChirps*numADCSamples/numTx);
adcData2 = zeros(numRX,numChirps*numADCSamples/numTx);
adcData3 = zeros(numRX,numChirps*numADCSamples/numTx); 
for row = 1 : numRX
    num = 1;
    for i = 1 : numChirps
        if mod(i,numTx) == 1
            %adcData1(row, (i-1)*numADCSamples+1:i*numADCSamples) = LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
            adcData1(row, (num-1) * numADCSamples + 1 : num * numADCSamples) = LVDS(i, (row-1) * numADCSamples + 1 : row * numADCSamples);
        end
        if mod(i,numTx) == 2
            adcData2(row, (num-1) * numADCSamples + 1 : num * numADCSamples) = LVDS(i, (row-1) * numADCSamples + 1 : row * numADCSamples);          
        end
        if mod(i,numTx) == 0
            adcData3(row, (num-1) * numADCSamples + 1 : num * numADCSamples) = LVDS(i, (row-1) * numADCSamples + 1 : row * numADCSamples);
            num = num + 1;
        end
    end
end

numChirps_per_Tx = numChirps/numTx;

%%%%% initialize the Tx matrix, each matrix is 3-dimension, 3rd dimension is the Rx 
T1_Rx_chirp = zeros(numChirps_per_Tx,numADCSamples);
T2_Rx_chirp = zeros(numChirps_per_Tx,numADCSamples);
T3_Rx_chirp = zeros(numChirps_per_Tx,numADCSamples);

T1_Rx_FFT = zeros(numChirps_per_Tx,numADCSamples);
T2_Rx_FFT = zeros(numChirps_per_Tx,numADCSamples);
T3_Rx_FFT = zeros(numChirps_per_Tx,numADCSamples);
for i = 2:4
   T1_Rx_chirp(:,:,i) = zeros(numChirps_per_Tx,numADCSamples);
   T2_Rx_chirp(:,:,i) = zeros(numChirps_per_Tx,numADCSamples);
   T3_Rx_chirp(:,:,i) = zeros(numChirps_per_Tx,numADCSamples);
   
   T1_Rx_FFT(:,:,i) = zeros(numChirps_per_Tx,numADCSamples);
   T2_Rx_FFT(:,:,i) = zeros(numChirps_per_Tx,numADCSamples);
   T3_Rx_FFT(:,:,i) = zeros(numChirps_per_Tx,numADCSamples);
end

%%%% fill in the value for each point
for j = 1: 4
    for i = 1: numChirps_per_Tx
        T1_Rx_chirp(i,1:numADCSamples,j) = adcData1(j,(i-1)*numADCSamples+1:i*numADCSamples);
        T2_Rx_chirp(i,1:numADCSamples,j) = adcData2(j,(i-1)*numADCSamples+1:i*numADCSamples);
        T3_Rx_chirp(i,1:numADCSamples,j) = adcData3(j,(i-1)*numADCSamples+1:i*numADCSamples);
    end
end
disp('raw data acquired.....');

%% Range-FFT with Hanning window
Hann_win = hann(numADCSamples);

for j = 1: 4
    T1_Rx_FFT(:,:,j) = fft(T1_Rx_chirp(:,:,j).* Hann_win',[],2);
    T2_Rx_FFT(:,:,j) = fft(T2_Rx_chirp(:,:,j).* Hann_win',[],2);
    T3_Rx_FFT(:,:,j) = fft(T3_Rx_chirp(:,:,j).* Hann_win',[],2);
end

%  figure; mesh(abs(T1_Rx_FFT(:,1:cut_range_bin,1))'); title("T1 Rx1 FFT 1:cutRangebin");view(2);
%  shading interp; colormap(parula); view(2); ylim([1,cut_range_bin]);