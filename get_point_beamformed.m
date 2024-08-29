function Y_beamformed_fft = get_point_beamformed(T1_Rx_chirp,T3_Rx_chirp,cut_range_bin,cand_bin,AoA_theta)

Slope = 50.018e12;
c = 3e8;
Chirp_fs = 1e7;
numADCSamples = 256;
angle_num = 181;
num_antenna = 8;

chirp_fft_resol = Chirp_fs / numADCSamples;
range = ((((0:255))*chirp_fft_resol)*c) / (2*Slope);

beam_seg_start = 1;
beam_seg_end = floor(length(T1_Rx_chirp)/2) * 2;  %length(T1_Rx_chirp);5460，2730
beam_seg_len = beam_seg_end - beam_seg_start + 1; %% only use part of signal for beamforming to save computation

%% 
Y_beamformed_fft = zeros(beam_seg_len,1);

%%%%% steering vector construction
omega_vector = zeros(angle_num,num_antenna);
phi = zeros(angle_num,1);
for theta = 1:angle_num
    for i = 1:num_antenna
        phi(theta) = (i-1) * pi * sin(((theta-1)-((angle_num-1)/2)) * pi / 180); %
        omega_vector(theta,i) = exp(-1j * phi(theta));
    end
end

theta = AoA_theta + 91;
    for t = beam_seg_start:beam_seg_end
        for theta_ant = 1:num_antenna
            if theta_ant <= 4
                Y_beamformed_fft(t-beam_seg_start+1) = Y_beamformed_fft(t-beam_seg_start+1) + omega_vector(theta,theta_ant) * T1_Rx_chirp(t,cand_bin,theta_ant);
            else
                Y_beamformed_fft(t-beam_seg_start+1) = Y_beamformed_fft(t-beam_seg_start+1) + omega_vector(theta,theta_ant) * T3_Rx_chirp(t,cand_bin,theta_ant-4);
            end
        end
    end

