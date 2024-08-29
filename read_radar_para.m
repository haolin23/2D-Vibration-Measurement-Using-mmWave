%%
adc_data_20240603_20mm_150cm_para_50_radar_disp = []; % 1*56
adc_data_20240603_20mm_150cm_para_50_radar_freq = [];
distance = 1900;

%% Get chirp data
%             close all;
%             clear;
%             clearvars -except adc_data20240504_20mm_los2m_perp_exp36789_radar_disp adc_data20240504_20mm_los2m_perp_exp36789_radar_freq;

for wii = 1:7
    clearvars -except wii adc_data_20240603_20mm_150cm_para_50_radar_freq adc_data_20240603_20mm_150cm_para_50_radar_disp distance;
                freq_0 = 77e9; 
                lambda = 3.89e-3;
                Fs = 128/0.024/3; % periodicity = 24 ，slow sample rate 1777
                cut_range_bin = 25; % number of range bins shown in the result
                low_cut_range_bin = 14; %5
                Chirp_fs = 1e7; 
                numADCSamples = 256;
                c = 3e8; Slope = 50.018e12;
                chirp_fft_resol = Chirp_fs / numADCSamples;

                filename = sprintf('adc_data_20240603_20mm_260cm_para_50_%d.bin',wii);   
%                 filename = sprintf('adc_data_20240519_20mm_perp_170_5.bin');   

                [T1_Rx_chirp,T1_Rx_FFT,T2_Rx_chirp,T2_Rx_FFT,T3_Rx_chirp,T3_Rx_FFT] = get_chirp_data(filename,cut_range_bin);

    %% select T1_RX_FFT and T3_RX_FFT
    T1_Rx_FFT_all = T1_Rx_FFT; T2_Rx_FFT_all = T2_Rx_FFT; T3_Rx_FFT_all = T3_Rx_FFT;
    T1_Rx_chirp_all = T1_Rx_chirp; T2_Rx_chirp_all = T2_Rx_chirp; T3_Rx_chirp_all = T3_Rx_chirp;
    for ni = 1:8
        close all;
        clearvars -except wii ni T1_Rx_FFT_all T2_Rx_FFT_all T3_Rx_FFT_all T1_Rx_chirp_all T2_Rx_chirp_all T3_Rx_chirp_all adc_data_20240603_20mm_150cm_para_50_radar_freq adc_data_20240603_20mm_150cm_para_50_radar_disp distance;

        freq_0 = 77e9; 
        lambda = 3.89e-3;
        Fs = 1777; % periodicity = 24 ，slow sample rate
        cut_range_bin = 32; % number of range bins shown in the result
        low_cut_range_bin = 20; %5
        Chirp_fs = 1e7; 
        numADCSamples = 256;
        c = 3e8; Slope = 50.018e12;
        chirp_fft_resol = Chirp_fs / numADCSamples;

        if ni ~= 8  
            T1_Rx_FFT = T1_Rx_FFT_all(10922*(ni-1)+1:10922*ni,:,:);
            T2_Rx_FFT = T2_Rx_FFT_all(10922*(ni-1)+1:10922*ni,:,:);
            T3_Rx_FFT = T3_Rx_FFT_all(10922*(ni-1)+1:10922*ni,:,:);
            T1_Rx_chirp = T1_Rx_chirp_all(10922*(ni-1)+1:10922*ni,:,:);
            T2_Rx_chirp = T2_Rx_chirp_all(10922*(ni-1)+1:10922*ni,:,:);
            T3_Rx_chirp = T3_Rx_chirp_all(10922*(ni-1)+1:10922*ni,:,:);
        else % ni == 8
            if 10922*8 > length(T1_Rx_FFT_all) % length(8th) < 8890
                T1_Rx_FFT = T1_Rx_FFT_all(10922*(ni-1)+1:end,:,:);
                T2_Rx_FFT = T2_Rx_FFT_all(10922*(ni-1)+1:end,:,:);
                T3_Rx_FFT = T3_Rx_FFT_all(10922*(ni-1)+1:end,:,:);
                T1_Rx_chirp = T1_Rx_chirp_all(10922*(ni-1)+1:end,:,:);
                T2_Rx_chirp = T2_Rx_chirp_all(10922*(ni-1)+1:end,:,:);
                T3_Rx_chirp = T3_Rx_chirp_all(10922*(ni-1)+1:end,:,:);
            else
                T1_Rx_FFT = T1_Rx_FFT_all(10922*(ni-1)+1:10922*ni,:,:);
                T2_Rx_FFT = T2_Rx_FFT_all(10922*(ni-1)+1:10922*ni,:,:);
                T3_Rx_FFT = T3_Rx_FFT_all(10922*(ni-1)+1:10922*ni,:,:);
                T1_Rx_chirp = T1_Rx_chirp_all(10922*(ni-1)+1:10922*ni,:,:);
                T2_Rx_chirp = T2_Rx_chirp_all(10922*(ni-1)+1:10922*ni,:,:);
                T3_Rx_chirp = T3_Rx_chirp_all(10922*(ni-1)+1:10922*ni,:,:);
            end
        end


%                 figure;mesh(abs(T1_Rx_FFT(:,1:cut_range_bin,1)));view(2);title("T1 Rx1 FFT 1:cutrangebin");
    %             vfft = get_Doppler_FFT(T1_Rx_chirp);

                %% Beamforming for localization
    % %             Y_sbeamformed = get_beamformed(T1_Rx_chirp,T3_Rx_chirp,cut_range_bin);
                Y_beamformed = get_beamformed_jet(T1_Rx_chirp,T3_Rx_chirp,cut_range_bin,low_cut_range_bin,1100);

    chirps_num = floor(length(T1_Rx_FFT)/2) * 2;
    % chirps_num = 1100;

    %% Locate range bin and angle
    ep = 1; % expand point
    ram = squeeze(Y_beamformed(ep,low_cut_range_bin:cut_range_bin,-60+91:60+91));
    ram = abs(ram);
    R = max(ram.');
    [rpks,rlocs] = findpeaks(R);
%     figure;plot(R);title("R");
    rlocs = rlocs + low_cut_range_bin-1;
    [srpks,sr_indices] = sort(rpks,'descend');
    two_range_locs = rlocs(sr_indices(1:2));

    A = max(ram);
    [apks,alocs] = findpeaks(A);
%     figure;plot(-60:1:60,A);title("A");
    alocs = alocs - 91 + 30;
    [sapks,sa_indices] = sort(apks,'descend'); % indices(i)表示排序后的第i个元素在原向量中的索引
    two_angle_locs = alocs(sa_indices(1:2));

    [ascend_rangelocs,srangeindices] = sort(two_range_locs); % 
    ascend_anglelocs = two_angle_locs(srangeindices(1:2)); % ascend_anglelocs, ascend_rangelocs, 按距离增加的两个目标
    
    %% MUSIC
    % disp('start MUSIC...')
    num_antenna = 4;
    % %%% select out the candidate range bin and concantanate 8 antennas vectors
    sel_index = ascend_rangelocs(1); %16,24 %15，19 %17，14，6
    T1_sel_bin = squeeze(T1_Rx_FFT(:,sel_index,:)); T3_sel_bin = squeeze(T3_Rx_FFT(:,sel_index,:));
    Ant_8_sel_bin = [T1_sel_bin];
    P = 1; %3
    d=lambda/2;
    Cov = Ant_8_sel_bin' * Ant_8_sel_bin;
    [U,V] = eig(Cov);
    UU=U(:,1:num_antenna-P);% noise sub space
    theta=-90:0.5:90;
    for ii=1:length(theta)
        AA=zeros(1,length(num_antenna)); %%% AA is steeting vector
        for jj=0:num_antenna-1
            AA(1+jj)=exp(-1j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda); % -1j
        end
        WW=AA*UU*UU'*AA'; 
        Pmusic(ii)=abs(1/WW);
    end
    Pmusic=10*log10(Pmusic/max(Pmusic));%spatial spectrum
    Pmusic = fliplr(Pmusic);
%     figure;plot(theta,Pmusic,'-k','linewidth',2.0);
%     xlabel('DOA \theta/degree');
%     ylabel('Power Spectrum/dB');
%     title('MUSIC algorithm for r2');
%     grid on;
    [~,angle_r1] = max(Pmusic);
    angle_r1 = -90 + 0.5*(angle_r1-1);
    disp('end MUSIC...');
    %% MUSIC r2
    % disp('start MUSIC...')
    clear AA WW Pmusic;
    num_antenna = 4;
    % %%% select out the candidate range bin and concantanate 8 antennas vectors
    sel_index = ascend_rangelocs(2); %16,24 %15，19 %17，14，6
    T1_sel_bin = squeeze(T1_Rx_FFT(:,sel_index,:)); T3_sel_bin = squeeze(T3_Rx_FFT(:,sel_index,:));
    Ant_8_sel_bin = [T1_sel_bin];
    P = 1; %2
    d=lambda/2;
    Cov = Ant_8_sel_bin' * Ant_8_sel_bin;
    [U,V] = eig(Cov);
    UU=U(:,1:num_antenna-P);% noise sub space
    theta=-90:1:90;
    for ii=1:length(theta)
        AA=zeros(1,length(num_antenna)); %%% AA is steeting vector
        for jj=0:num_antenna-1
            AA(1+jj)=exp(-1j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda); % -1j
        end
        WW=AA*UU*UU'*AA'; 
        Pmusic(ii)=abs(1/WW);
    end
    Pmusic=10*log10(Pmusic/max(Pmusic));%spatial spectrum
    Pmusic = fliplr(Pmusic);
%     figure;plot(theta,Pmusic,'-k','linewidth',2.0);
%     xlabel('DOA \theta/degree');
%     ylabel('Power Spectrum/dB');
%     title('MUSIC algorithm for r2');
%     grid on;
    [~,angle_r2] = max(Pmusic);
    angle_r2 = -90 + 1*(angle_r2-1);
    disp('end MUSIC...');

    theta1 = angle_r2 - angle_r1;

    %% FFT and Chirp-Z on 51th chirp
    xx = T1_Rx_chirp(ep,1:numADCSamples,1);
    N = length(xx);
    fft_result = fft(xx, N); % 执行FFT 
    fft_frequencies = (0:N-1) * (Chirp_fs/N); % 计算频率   
    Y_single_side = fft_result(1:N/2+1);  % 仅保留单侧频谱  
    Y_single_side = Y_single_side / N;  % 调整幅度 
%     % figure;
%     % plot(fft_frequencies(1:N/2+1), abs(Y_single_side));
%     % title('fft on 851th chirp,x-axis is frequence');
%     figure;
%     plot(fft_frequencies(1:N/2+1)*c/(2*Slope), abs(Y_single_side));
%     title('fft on th chirp,x-axis is range(m)');
    
    %% czt 
     x = T1_Rx_chirp(ep,:,1); % 2730 * 256
    M = 256;
    czt_start_range = 2.3; % m
    czt_delta_range = 0.005; % m
    czt_delta_freq = 2 * Slope * czt_delta_range / c;
    czt_start_freq = 2 * Slope * czt_start_range / c;
    W = exp(-1i * 2 * pi * czt_delta_freq / Chirp_fs);
    A = exp(1i * 2 * pi * czt_start_freq / Chirp_fs);
    x_czt = czt(x,M,W,A);
    % x_czt = x_czt - mean(x_czt);
    f = czt_start_freq : czt_delta_freq : czt_start_freq+czt_delta_freq*(M-1);
    range_czt = f*c/(2*Slope);
    figure; plot(range_czt,abs(x_czt));title('czt');
    abs_czt = abs(x_czt);
    [peaks,peak_locs] = findpeaks(abs_czt);
    [sorted_peaks,sorted_indices]=sort(peaks,'descend');
    top_two_peak_locs = peak_locs(sorted_indices(1:2));
    top_two_ranges = [range_czt(top_two_peak_locs(1)),range_czt(top_two_peak_locs(2))];
    two_ranges = sort(top_two_ranges);
    % b = four_ranges(1)*2;
    r1_0 = two_ranges(1);
    r2_0 = two_ranges(2);

    tb = (r1_0^2-r2_0^2)/(2*r1_0*cosd(theta1)-2*r2_0);
    ta = r2_0 - tb;
    r12 = sqrt(r1_0^2+r2_0^2-2*r1_0*r2_0*((r1_0^2+ta^2-tb^2)/(2*r1_0*ta)));
    b = r12;

                    %% ------extract the phase of the beamformed points r1----------
                    cand_bin = ascend_rangelocs(1); %17
                    AoA_theta  = ascend_anglelocs(1); % -33 % -31

        Y_point_beamformed = get_point_beamformed(T1_Rx_FFT,T3_Rx_FFT,cut_range_bin,cand_bin,AoA_theta);
        focus_chirps = Y_point_beamformed;

                    %%-- plot raw phase --
%                     figure;
                    temp_phase = angle(focus_chirps);
                    focus_phase_extract = unwrap(temp_phase);
%                     plot(focus_phase_extract);title("raw phase");
%                     figure;plot(real(focus_chirps),imag(focus_chirps),'.');title("raw phase IQ r1");hold on;axis equal;

                    % calibrate phase
                    x_temp = real(focus_chirps);
                    y_temp = imag(focus_chirps);
                    xy = [x_temp(:), y_temp(:)]; 

                    ParIni = TaubinSVD(xy);
                    Par = LM(xy, ParIni,1); 
                    xCenter = Par(1);
                    yCenter = Par(2);
                    radius = Par(3);

                    ctheta = linspace(0, 2*pi, 100);    
                    cx = xCenter + radius * cos(ctheta);  
                    cy = yCenter + radius * sin(ctheta);  
%                     plot(cx, cy,'r');  

                    x_temp = x_temp(1:chirps_num);y_temp = y_temp(1:chirps_num);% 2730
                    temp_processed = (x_temp - xCenter) + 1j*(y_temp - yCenter);
                    temp_phase_processed = zeros(1,length(temp_processed));
                    for i = 1:length(temp_processed)
                        temp_phase_processed(i) = angle(temp_processed(i));
                    end
                    temp_phase_processed = unwrap(temp_phase_processed);

%                     figure;
%                     plot(real(temp_processed),imag(temp_processed),'.'); 
%                     title("phase IQ after fitcircle r1");
%                     axis equal;hold on;
                    cctheta = linspace(0, 2*pi, 100); 
                    ccx = radius * cos(cctheta); 
                    ccy = radius * sin(cctheta);  
%                     plot(ccx, ccy,'r');  

%         figure;plot(abs(temp_processed));title("abs r1");

        % get frequencies
        vib_fft = fft(temp_processed,chirps_num);
        frequencies = (0:chirps_num-1) * (Fs/chirps_num);
        Y_single_side = vib_fft(1:chirps_num/2+1);
        Y_single_side = Y_single_side / chirps_num;
        Y_single_side(1) = 0;
%         figure;
%         plot(frequencies(1:chirps_num/2+1),abs(Y_single_side));
%         title("FFT on focuschirptempprocessed1");

%                     figure;plot(temp_phase_processed);title("phase after fitcircle");

        %             first_phase = temp_phase_processed(1);
        %             phase_delta = temp_phase_processed - first_phase;
                    ep_phase = temp_phase_processed(ep);
                    phase_delta = temp_phase_processed - ep_phase;
                    phase_delta = phase_delta(ep:end);

                    delta_r1 = phase_delta * c / (4 * pi * freq_0);
%                     figure;plot(delta_r1);title("vibration r1");

                    %% average peak_to_peak value
                    [thp_val,thp_loc] = findpeaks(delta_r1,'MinPeakProminence',1e-4,'MinPeakDistance',distance); % 1e-3
                    [tlp_val,tlp_loc] = findpeaks(-delta_r1,'MinPeakProminence',1e-4,'MinPeakDistance',distance);

                    hp_val = thp_val;hp_loc = thp_loc;lp_val = tlp_val;lp_loc = tlp_loc;

                    ptp = [];
                    if length(hp_val) <= length(lp_val) 
                        pindex = length(hp_val);
                    else
                        pindex = length(lp_val);
                    end
                    for i = 1 : pindex % from 2 to pindex-1
                        ptp = [ptp,hp_val(i)+lp_val(i)];
                    end
                    ptp_avg = mean(ptp);
%                     figure; plot(delta_r1); title("vibration r1");hold on; plot(hp_loc,hp_val,'^'); plot(lp_loc,-lp_val,'o');
                    disp("average peak_to_peak of r1:(mm)");disp(ptp_avg*1000);

                    % get frequencies
                    chirps_num = chirps_num-ep+1;
                    vib_fft = fft(delta_r1,chirps_num);
                    frequencies = (0:chirps_num-1) * (Fs/chirps_num);
                    Y_single_side = vib_fft(1:chirps_num/2+1);
                    Y_single_side = Y_single_side / chirps_num;
                    Y_single_side(1) = 0;
%                     figure;
%                     plot(frequencies(1:chirps_num/2+1),abs(Y_single_side));
                    chirps_num = chirps_num+ep-1;
%                     title("FFT on r1");

          %% ------extract the phase of the beamformed points r2----------
                    cand_bin = ascend_rangelocs(2); 
                    AoA_theta  = ascend_anglelocs(2);
        %             focus_chirps = Y_besamformed(:,cand_bin,AoA_theta + 91);
        Y_point_beamformed = get_point_beamformed(T1_Rx_FFT,T3_Rx_FFT,cut_range_bin,cand_bin,AoA_theta);
        focus_chirps = Y_point_beamformed;
        % focus_chirps = T1_Rx_FFT(:,cand_bin,1);

                    %%-- plot raw phase --
%                     figure;
                    temp_phase = angle(focus_chirps);
                    focus_phase_extract = unwrap(unwrap(temp_phase));
%                     plot(focus_phase_extract);title("raw phase");

%                     figure;plot(real(focus_chirps),imag(focus_chirps),'.');title("raw phase IQ"); hold on;
%                     axis equal;

                    % calibrate phase
                    x_temp = real(focus_chirps);
                    y_temp = imag(focus_chirps);
                    xy = [x_temp(:), y_temp(:)]; 

                    ParIni = TaubinSVD(xy);
                    Par = LM(xy, ParIni,1); 
                    xCenter = Par(1);
                    yCenter = Par(2);
                    radius = Par(3);

                    ctheta = linspace(0, 2*pi, 100);    
                    cx = xCenter + radius * cos(ctheta);  
                    cy = yCenter + radius * sin(ctheta);  
%                     plot(cx, cy,'r');  

                    x_temp = x_temp(1:chirps_num);y_temp = y_temp(1:chirps_num);% 2730
                    temp_processed = (x_temp - xCenter) + 1j*(y_temp - yCenter);
                    temp_phase_processed = zeros(1,length(temp_processed));
                    for i = 1:length(temp_processed)
                        temp_phase_processed(i) = angle(temp_processed(i));
                    end
                    temp_phase_processed = unwrap(temp_phase_processed);

%                     figure;plot(real(temp_processed),imag(temp_processed),'.');
%                     title("phase IQ after fitcircle r2");
%                     axis equal;hold on;
                    cctheta = linspace(0, 2*pi, 100); 
                    ccx = radius * cos(cctheta); 
                    ccy = radius * sin(cctheta);  
%                     plot(ccx, ccy,'r');  

%                     figure;plot(temp_phase_processed);title("phase after fitcircle");

        %             first_phase = temp_phase_processed(1);
        %             phase_delta = temp_phase_processed - first_phase;
                    ep_phase = temp_phase_processed(ep);
                    phase_delta = temp_phase_processed - ep_phase;
                    phase_delta = phase_delta(ep:end);

                    delta_r2 = phase_delta * c / (4 * pi * freq_0);
%                     figure;plot(delta_r2);title("vibration r2");

                    % average peak_to_peak value
                    [thp_val,thp_loc] = findpeaks(delta_r2,'MinPeakProminence',1e-4,'MinPeakDistance',distance);
                    [tlp_val,tlp_loc] = findpeaks(-delta_r2,'MinPeakProminence',1e-4,'MinPeakDistance',distance);

                    hp_val = thp_val;hp_loc = thp_loc;lp_val = tlp_val;lp_loc = tlp_loc;

                    ptp = [];
                    if length(hp_val) <= length(lp_val) 
                        pindex = length(hp_val);
                    else
                        pindex = length(lp_val);
                    end
                    for i = 1 : pindex % from 2 to pindex-1
                        ptp = [ptp,hp_val(i)+lp_val(i)];
                    end
                    ptp_avg = mean(ptp);
%                     figure; plot(delta_r2); hold on; title("vibration r2"); plot(hp_loc,hp_val,'^'); plot(lp_loc,-lp_val,'o');
                    disp("average peak_to_peak of r2:(mm)");disp(ptp_avg*1000);

                    % get frequencies
                    chirps_num = chirps_num-ep+1;
                    vib_fft = fft(delta_r2,chirps_num);
                    frequencies = (0:chirps_num-1) * (Fs/chirps_num);
                    Y_single_side = vib_fft(1:chirps_num/2+1);
                    Y_single_side = Y_single_side / chirps_num;
                    Y_single_side(1) = 0;
%                     figure;
%                     plot(frequencies(1:chirps_num/2+1),abs(Y_single_side));
                    chirps_num = chirps_num+ep-1;
%                     title("FFT on r2");grid on;

        %% Calculate axial displacement
        syms r1 r2;
        f1(r1,r2) = 2*(1/b)*(1/4)*sqrt((r1+r2+b)*(r1+r2-b)*(r1+b-r2)*(r2+b-r1));
        f2(r1,r2) = sqrt(r1^2-f1^2);
        f3(r1,r2) = b - sqrt(r2^2-f1^2); % an alternative representation of y
        f1_r1 = diff(f1,r1);    f1_r2 = diff(f1,r2);
        f2_r1 = diff(f2,r1);    f2_r2 = diff(f2,r2);
        f3_r1 = diff(f3,r1);    f3_r2 = diff(f3,r2);

        f1_r1_at_point = subs(f1_r1,[r1,r2],[r1_0,r2_0]);
        f1_r2_at_point = subs(f1_r2,[r1,r2],[r1_0,r2_0]);
        f2_r1_at_point = subs(f2_r1,[r1,r2],[r1_0,r2_0]);
        f2_r2_at_point = subs(f2_r2,[r1,r2],[r1_0,r2_0]);
        f3_r1_at_point = subs(f3_r1,[r1,r2],[r1_0,r2_0]);
        f3_r2_at_point = subs(f3_r2,[r1,r2],[r1_0,r2_0]);

        % first-order expansion
        delta_f1 = f1_r1_at_point * delta_r1 + f1_r2_at_point * delta_r2;
        delta_f2 = f2_r1_at_point * delta_r1 + f2_r2_at_point * delta_r2;
        % delta_f3 = f3_r1_at_point * delta_r1 + f3_r2_at_point * delta_r2;

        delta_f1_double = double(delta_f1);
        delta_f2_double = double(delta_f2);
        % delta_f3_double = double(delta_f3);
        % delta_f23avg_double = (delta_f2_double - delta_f3_double)/2;

%         figure;plot(delta_f1_double);title("delta f1");
%         figure;plot(delta_f2_double);title("delta f2");
        % figure;plot(delta_f3_double);title("delta f3");
        % figure;plot(delta_f23avg_double);title("delta f23avg");

        %% average peak_to_peak value
        % delta_f1_double = double(delta_f1);
        [thp_val,thp_loc] = findpeaks(delta_f1_double,'MinPeakProminence',1e-4,'MinPeakDistance',distance);
        [tlp_val,tlp_loc] = findpeaks(-delta_f1_double,'MinPeakProminence',1e-4,'MinPeakDistance',distance);

        hp_val = thp_val;hp_loc = thp_loc;lp_val = tlp_val;lp_loc = tlp_loc;

        ptp = [];
        if length(hp_val) <= length(lp_val) 
            pindex = length(hp_val);
        else
            pindex = length(lp_val);
        end
        for i = 1 : pindex 
            ptp = [ptp,hp_val(i)+lp_val(i)];
        end
        ptp_avg = mean(ptp);
%         figure; plot(delta_f1_double); hold on; plot(hp_loc,hp_val,'^'); plot(lp_loc,-lp_val,'o');
%         title("vibration of parallel axis");
        disp("average peak_to_peak of parallel axis:(mm)");disp(ptp_avg*1000);
        adc_data_20240603_20mm_150cm_para_50_radar_disp = [adc_data_20240603_20mm_150cm_para_50_radar_disp,ptp_avg*1000]; %mm

        % delta_f2_double = double(delta_f2);
        [thp_val,thp_loc] = findpeaks(delta_f2_double,'MinPeakProminence',1e-4,'MinPeakDistance',distance);
        [tlp_val,tlp_loc] = findpeaks(-delta_f2_double,'MinPeakProminence',1e-4,'MinPeakDistance',distance);

        hp_val = thp_val;hp_loc = thp_loc;lp_val = tlp_val;lp_loc = tlp_loc;

        ptp = [];
        if length(hp_val) <= length(lp_val) 
            pindex = length(hp_val);
        else
            pindex = length(lp_val);
        end
        for i = 1 : pindex 
            ptp = [ptp,hp_val(i)+lp_val(i)];
        end
        ptp_avg = mean(ptp);
%         figure; plot(delta_f2_double); hold on; plot(hp_loc,hp_val,'^'); plot(lp_loc,-lp_val,'o');
%         title("vibration of perpendicular axis");
        disp("average peak_to_peak of perpendicular:(mm)");disp(ptp_avg*1000);
%         adc_data_20240603_20mm_150cm_para_50_radar_disp = [adc_data_20240603_20mm_150cm_para_50_radar_disp,ptp_avg*1000]; %mm

        % trace
        figure;
        plot(delta_f1_double, delta_f2_double, '.');  
        xlabel('parallel轴位移');  
        ylabel('perpendicular轴位移');  
        titleString = sprintf('物体轨迹 %d',(wii-1)*8+ni);
        title(titleString);  
        axis equal;
        disp((wii-1)*8+ni);

        %% FFT
        chirps_num = chirps_num-ep+1;
        % get frequencies of x
        vib_fft = fft(delta_f1_double,chirps_num);
        frequencies = (0:chirps_num-1) * (Fs/chirps_num);
        Y_single_side = vib_fft(1:chirps_num/2+1);
        Y_single_side = Y_single_side / chirps_num;
        Y_single_side(1) = 0;
        [maxY, indexOfMax] = max(Y_single_side);  
        correspondingFrequency = frequencies(indexOfMax);
%         adc_data_20240603_20mm_150cm_para_50_radar_freq = [adc_data_20240603_20mm_150cm_para_50_radar_freq,correspondingFrequency];
%         figure;
%         plot(frequencies(1:chirps_num/2+1),abs(Y_single_side));
%         title("FFT on parallel");
        M = chirps_num;
        startf = 0.5; intervalf=0.0005;
        W = exp(-1j*2*pi*intervalf/Fs);
        A = exp(1j*2*pi*startf/Fs);
        czt_result = czt(delta_f1_double,M,W,A);
        ff = startf:intervalf:startf+intervalf*(M-1);
%         figure;plot(ff,abs(czt_result));
        [maxY, indexOfMax] = max(czt_result);  
        corf = ff(indexOfMax);
        adc_data_20240603_20mm_150cm_para_50_radar_freq = [adc_data_20240603_20mm_150cm_para_50_radar_freq,corf];
%         title('czt');

        % get frequencies of y
        vib_fft = fft(delta_f2_double,chirps_num);
        frequencies = (0:chirps_num-1) * (Fs/chirps_num);
        Y_single_side = vib_fft(1:chirps_num/2+1);
        Y_single_side = Y_single_side / chirps_num;
        Y_single_side(1) = 0;
        [maxY, indexOfMax] = max(Y_single_side);  
        correspondingFrequency = frequencies(indexOfMax);
%         adc_data_20240603_20mm_150cm_para_50_radar_freq = [adc_data_20240603_20mm_150cm_para_50_radar_freq,correspondingFrequency];
%         figure;
%         plot(frequencies(1:chirps_num/2+1),abs(Y_single_side));
%         title("FFT on perpendicular");

    end
end

%%
% save('czt_0603_20mm_para_50_radar_freq_260cm.mat','adc_data_20240603_20mm_150cm_para_50_radar_freq');

figure;plot(adc_data_20240603_20mm_150cm_para_50_radar_freq,'.');





