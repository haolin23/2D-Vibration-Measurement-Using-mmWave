function Y_beamformed_fft = get_beamformed(T1_Rx_chirp,T3_Rx_chirp,cut_range_bin,low_cut_range_bin,chirps_num)

Slope = 50.018e12;
c = 3e8;
Chirp_fs = 1e7;
numADCSamples = 256;
angle_num = 181;
num_antenna = 8;

chirp_fft_resol = Chirp_fs / numADCSamples;
range = ((((0:255))*chirp_fft_resol)*c) / (2*Slope);

beam_seg_start = 1;
% beam_seg_end = floor(length(T1_Rx_chirp)/2) * 2;  %length(T1_Rx_chirp);5460，2730
beam_seg_end = chirps_num;
beam_seg_len = beam_seg_end - beam_seg_start + 1; %% only use part of signal for beamforming to save computation

Y_beamformed = zeros(beam_seg_len,numADCSamples);
for theta = 2:angle_num
    Y_beamformed(:,:,theta) = zeros(beam_seg_len,numADCSamples);
end
disp('start beamforming...');

%%%%% steering vector construction
omega_vector = zeros(angle_num,num_antenna);
phi = zeros(angle_num,1);
for theta = 1:angle_num
    for i = 1:num_antenna
        phi(theta) = (i-1) * pi * sin(((theta-1)-((angle_num-1)/2)) * pi / 180); %
        omega_vector(theta,i) = exp(-1j * phi(theta));
    end
end

for theta = 1:angle_num 
    for t = beam_seg_start:beam_seg_end
        for theta_ant = 1:num_antenna
            if theta_ant <= 4
                Y_beamformed(t-beam_seg_start+1,:,theta) = Y_beamformed(t-beam_seg_start+1,:,theta) + omega_vector(theta,theta_ant) * T1_Rx_chirp(t,:,theta_ant);
            else
                Y_beamformed(t-beam_seg_start+1,:,theta) = Y_beamformed(t-beam_seg_start+1,:,theta) + omega_vector(theta,theta_ant) * T3_Rx_chirp(t,:,theta_ant-4);
            end
        end
    end
end

%%%% FFT on beamformed output
Y_beamformed_fft = zeros(beam_seg_len,numADCSamples);
for theta = 2:angle_num
    Y_beamformed_fft(:,:,theta) = zeros(beam_seg_len,numADCSamples);
end

Hann_win = hann(numADCSamples);
for theta = 1:angle_num
    Y_beamformed_fft(:,:,theta) = fft(Y_beamformed(:,:,theta).* Hann_win',[],2);
end

%% use this plot
angle = -90:1:90;
for i = 1:50:105
     Y_2d = squeeze(Y_beamformed_fft(i,low_cut_range_bin:cut_range_bin,:)); 
     figure; 
     h = surf(abs(Y_2d)); view(2); h.YData = range(low_cut_range_bin:cut_range_bin); h.XData = angle;
     shading interp; colormap('jet'); 
     colorbar; xlabel('Angle'); ylabel('Distance(m)'); xlim([-90,90]); %ylim([1,25]);
     title(sprintf("range - angle %d",i));
end

% old version
% for i = 1:20:400
%  %     Y_2d = squeeze(Y_beamformed_fft(i,1:20,:)); figure; mesh(abs(Y_2d));
%  %     view(2);
%      Y_2d = squeeze(Y_beamformed_fft(i,5:cut_range_bin,:)); 
%      figure; 
%      h = surf(abs(Y_2d)); view(2); h.YData = range(5:cut_range_bin); h.XData = angle;
%      shading interp; colormap('jet'); 
%      colorbar; xlabel('Angle'); ylabel('Distance(m)'); xlim([-90,90]); %ylim([1,25]);
%      title(sprintf("range - angle %d",i));
% end
 
 %% add
%  for i = 1:50:50
%  %     Y_2d = squeeze(Y_beamformed_fft(i,1:20,:)); figure; mesh(abs(Y_2d));
%  %     view(2);
%      Y_2d = squeeze(Y_beamformed_fft(i,1:15,:)); 
%      figure; 
%      h = surf(abs(Y_2d)); view(2); h.YData = range(1:15); h.XData = angle;
%      shading interp; colormap('jet'); 
%      colorbar; xlabel('Angle'); ylabel('Distance(m)'); xlim([-90,90]); %ylim([1,25]);
%      title('range - angle');
%  end

 disp('end beamforming...')
