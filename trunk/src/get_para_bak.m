function [paraEx paraSt] = get_para(network,freq,snapRate, snapNum, posBS,posMS,veloMS)
%GET_PARA Generate the external and stochastic parameters of the scenario
%Default call: [paraEx paraSt] = get_para(network,freq,snapRate, snapNum,
%posBS,posMS,veloMS)
%
%------
%Input:
%------
%network: network type, fours scenarios are available: 
%'macrocell', 'microcell', 'picocell', 'aalto'
%freq: frequency range, 2x1 vector, [freq_start freq_stop] [Hz]
%snapRate: channel snapshot rate, [s]
%snapeNum: channel snapshot number
%BSPos: position matrix of BSs, size = (number of BS, [x y z])  [m]
%MSPos: initial position matrix of MSs, size = (number of MS, [x y z])  [m]
%MSVelo: velocity vector of MSs, size = (number of MS, [x y z])  [m/s]
%
%------
%Output:
%------
%paraEx: external parameters
% .network: network type
% .freq_start: start frequency [Hz]
% .freq_stop: stop frequency [Hz]
% .freq: carrier frequency [Hz]
% .net_radii: network cell radius [m]
% .sample_rate: delay sampling rate [s]
% .c0: wave speed [m/s]
% .pos_BS: BS position [m]
% .pos_MS: MS position [m]
% .velo_MS: MS moving speed vector [m]
% .nfloor: average number of floors between BS and MS( picocell, aalto)
% .snap_rate: channel snapshot rate [s]
% .snap_num: channel snapshot number 
% .bandwidth: bandwidth [Hz]
% .delay_max: maximum delay, 5 times of net_radii [s]
%paraSt: stochastic parameters
% .r_c: visibility region radius for clusters [m]
% .l_c: transition region radius for clusters [m]
% .k_tau: cluster power attenuation coefficient [dB/mu_s]
% .tau_b: cut-off delay for cluster power attenuation [s]
% .r_l: visibility region radius for LOS [m]
% .l_l: transition region radius for LOS [m]
% .mu_k: LOS power factor mean 
% .sigma_k: LOS power factor dB scale std
% .n_c_local: average number of local clusters
% .k_sel: proportion of twin/single clusters
% .mu_n_c_far: average number of far clusters, mean
% .n_c_far: average number of far clusters in the scenario
% .n_mpc: number of MPCs per cluster
% .mu_tau: delay spread mean [s]
% .sigma_stau: delay spread dB scale std
% .mu_phi_BS: AoD mean [deg]
% .sigma_phi_BS: AoD dB std
% .mu_theta_BS: EoD mean [deg]
% .sigma_theta_BS: EoD dB scale std
% .mu_phi_MS: AoA mean [deg]
% .sigma_phi_MS: AoA dB scale std
% .mu_theta_MS: EoA mean [deg]
% .sigma_theta_MS: EoA dB scale std
% .sigma_sh: std of shadowing [dB]
% .rho: correlation matrix
% .corr_mat: Cholesky factorized correlation matrix
% .mu_phi_c: azimuth deviation of cluster to BS mean [deg]
% .sigma_phi_c: azimuth deviation of cluster to BS/VR std [deg]
% .mu_theta_c: azimuth deviation of cluster to BS/VR mean [deg]
% .sigma_theta_c: elevation deviation of cluster to BS/VR [deg]
% .mean_r_c: distance of cluster to BS/VR mean [m]
% .sigma_r_c: distance of cluster to BS/VR sigma [m]
% .BSCC: Ratio of BS-common cluster to total cluster ([2,3,...]BS-CC)
% .MSCC: average number of VRs in VR group
%
%See also: cost2100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng, ICTEAM, UCL, Belgium 
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paraEx.network = network; %Network type

switch network
    case 'macrocell'
        paraEx.freq_start = 900e6; %Reference start frequency
        paraEx.freq_stop = 2e9; %Reference stop frequency
        paraEx.h_BS = 50; %BS height
        paraEx.h_MS = 1.5; %MS height
        paraEx.net_radii = 1000; %cell radius
        paraEx.h_rooftop = 15; % rooftop height
        paraEx.w_r = 25; %road width 
        paraEx.w_b = 50; %Street width
        paraEx.phi_road = 45/180*pi; %road orientation
                
    case 'microcell'
        paraEx.freq_start = 1e9;%[Hz]
        paraEx.freq_stop = 5e9;%[Hz]
        paraEx.h_BS_low = 3; %BS height 3-10m
        paraEx.h_BS_high = 10; %BS height 3-10m
        paraEx.h_MS = 1.5; %MS height
        paraEx.net_radii = 1000*randn; %cell radius any
        paraEx.h_rooftop = 15; % rooftop height ???
        paraEx.w_r = 20; %road width ???
        paraEx.w_s = 50; %street width ???
        paraEx.phi_road = 45/180*pi; %road orientation ???
        
    case 'picocell'        
        paraEx.freq_start = 2e9;%[Hz]
        paraEx.freq_stop = 5e9;%[Hz]
        paraEx.h_MS = 2; %MS height
        paraEx.net_radii = 30; %cell radius
        paraEx.a_room = [2 3]; % room size [length width] ???
        paraEx.n_floor = 0;%round(4+(10-4)*rand); % floors between BS and MS ???
    case 'aalto' %OLOS scenario
        paraEx.freq_start = 5.3e9-60e6;%[Hz]
        paraEx.freq_stop = 5.3e9+60e6;%[Hz]
        paraEx.net_radii = 30; %cell radius [m]    
        paraEx.n_floor = 0;
    case 'test'
        paraEx.freq_start = 5.3e9-60e6;%[Hz] %Reference start frequency
        paraEx.freq_stop = 5.3e9+60e6;%[Hz] %Reference stop frequency
        paraEx.net_radii = 30; %cell radius [m]    
        paraEx.n_floor = 0; %Number of floors between BS and MS
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Common external parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.snap_rate = snapRate; %Channel snapshot rate [s]
paraEx.snap_num = snapNum; %Channel snapshot number
paraEx.c0 = 3e8; %Wave speed [m/s]

paraEx.freq = (freq(1)+freq(2))/2; %Carrier frequency [Hz]
paraEx.bandwidth = freq(2)-freq(1); %Bandwidth [Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BS/MS position & movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.num_BS = length(posBS(:,1));
paraEx.num_MS = length(posMS(:,1));
paraEx.pos_BS = posBS;
paraEx.pos_MS = posMS;
paraEx.velo_MS = veloMS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Delay resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.sample_rate = 1/paraEx.bandwidth; %Delay sample rate [s]
paraEx.delay_max = paraEx.net_radii*5/paraEx.c0; %Maximum delay 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get ST parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraSt = get_paraSt(paraEx);

end

function paraSt = get_paraSt(paraEx)
%Sub function to get the stochastic parameters

switch paraEx.network
    %%%%%%%%%%%%%%%%
    case 'macrocell'
    %%%%%%%%%%%%%%%%    
        paraSt.r_c = 100; %Visibility region size
        paraSt.l_c = 20; %Transition region size
        paraSt.k_tau = 1; %Cluster power
        paraSt.tau_b = 10e-6; %Cluster power
        
        paraSt.d_co = 500; %LOS cut-off distance
        paraSt.r_l = 30; %LOS visibility region
        paraSt.l_l = 20; %LOS transition region
        EPL = 0; %Excess pathloss
        paraSt.mu_k = ((26-EPL)/6); %LOS power factor, mean, dB
        paraSt.sigma_k = 6; %LOS power factor, std, dB
        
        paraSt.n_c_local = 1; %average number of local clusters
        paraSt.k_sel = 1; %Only single clusters
        paraSt.mu_n_c_far = 1.18; %average number of add. clusters(single interacting clusters)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        n_c_far = random('poiss',paraSt.mu_n_c_far); %Possion distributed number of clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n_c_far ==0
            n_c_far = 1;
        end       
        paraSt.n_c_far = n_c_far;

        paraSt.n_mpc = 20; %number of MPCs per cluster
        paraSt.k_mpc = 0; %Rice factor of additional clusters
        
        paraSt.mu_diff = 0.05; %Diffuse radiation mean
        paraSt.sigma_diff = 3.4; %Diffuse radiation std
        paraSt.pdp_tau = 0.5e-6; %PDP of the diffuse radiation, uniform in azimuth exp(-t/tau), tau = 0.5us
        
        paraSt.mu_tau = 0.4e-6; %Delay spread mean[s]
        paraSt.sigma_tau = 3; %Delay spread std, 3dB        
        paraSt.mu_phi_BS = 0.81; %AoD spread mean, in degree
        paraSt.sigma_phi_BS = 0.34; %AoD spread std, in degree
        paraSt.mu_theta_BS = 0.5; %EoD spread mean, in degree
        paraSt.sigma_theta_BS = 3; %EoD spread std, in degree
        paraSt.mu_phi_MS = 35; %AoA spread mean, in degree
        paraSt.sigma_phi_MS = 0; %AoA spread std, in degree        
        paraSt.mu_theta_MS = 0.81; %EoA spread mean, in degree
        paraSt.sigma_theta_MS = 0.34; %EoA spread std, in degree
        
%         paraSt.theta_MS_low = 0; 
%         paraSt.theta_MS_high = 45; %elevation spread at MS, uniform distributed, not used yet
        
        paraSt.sigma_sf = 6; %shadowing std[dB]
                
        %Cross-correlations
        paraSt.rho=[1 -0.6 -0.6 0   0   0;...        % shadow fading
                    -0.6 1   0.5 0   0   0;...        % delay spread
                    -0.6 0.5 1   0   0   0;...        % theta BS
                     0   0   0   1   0   0;...        % phi BS
                     0   0   0   0   1   0;...        % theta MS
                     0   0   0   0   0   1];          % phi MS   
        paraSt.corr_mat=chol(paraSt.rho);    % the correlation matrix is caculated using the Cholesky factorization
        
        %Autocorrelation distances, not used yet
        paraSt.l_s = 100;
        paraSt.l_tau = 100;
        paraSt.l_phi_BS = 100;
        paraSt.l_theta_BS = 100;
        paraSt.l_phi_MS = 100;
        paraSt.l_theta_MS = 100;
                
        %Other parameters
        paraSt.sigma_phi_c = 60/180*pi; %azimuth of cluster to visibility region , std
        paraSt.sigma_theta_c = 0; %elevation of cluster to visibility region , std
        paraSt.min_r_c = 400; %Minimum distance from cluster to VR
        paraSt.sigma_r_c = 50; %Distance from cluster to VR, std
        paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10);%LOS power factor
        
        %Polarization
        paraSt.mu_xpdv = 0;
        paraSt.sigma_xpdv = 0; %Mean and std for XPDV
        paraSt.mu_xpdh = 0; 
        paraSt.sigma_xpdh = 0; %Mean and std for XPDH
        paraSt.mu_cpr =0;
        paraSt.sigma_cpr = 0; %Mean and std for CPR
        
        %For dual link
        paraSt.ratio_common = 0; %Ratio of common cluster to total cluster
        paraSt.d_common = 0; %lifetime of common cluster, in distance (m)            
                
    %%%%%%%%%%%%%%%%    
    case 'microcell'
    %%%%%%%%%%%%%%%%
        paraSt.r_c = 50; %Visibility region size
        paraSt.l_c = 20; %Transition region size
        paraSt.k_tau = 40; %Cluster power
        paraSt.tau_b = 0.5; %Cluster power
        %Line of sight
        paraSt.d_co = 300; %LOS cut-off distance
        paraSt.r_l = 50;
        paraSt.l_l = 50;
        paraSt.mu_k = 7; %Power factor, mean, dB
        paraSt.sigma_k = 2.3; %Power factor, std, dB
        
        paraSt.n_c = 1; %average number of local clusters
        paraSt.k_sel = 0.5;        
        paraSt.n_c_add = 3; %average number of add. clusters(single interacting clusters)
        paraSt.n_mpc = 7; %number of MPCs per cluster
        paraSt.k_mpc = 2; %Rice factor of additional clusters
        
        paraSt.mu_diff = 0.05; %Diffuse radiation mean, not used yet
        paraSt.sigma_diff = 3.4; %Diffuse radiation std, not used yet
        paraSt.pdp_tau = 0.5e-6; %PDP of the diffuse radiation, uniform in azimuth exp(-t/tau), tau = 0.5us
        
        paraSt.mu_tau = 13e-6; %Delay spread mean[s]
        paraSt.sigma_tau = 14; %Delay spread std, dB
        paraSt.mu_phi_BS = 2.3; %AoD spread mean, in degree
        paraSt.sigma_phi_BS = 3.4; %AoD spread std, in degree
        paraSt.mu_theta_BS = 1.3; %EoD spread mean, in degree
        paraSt.sigma_theta_BS = 3.3; %EoD spread std, in degree
        
        paraSt.mu_phi_MS = 2.3; %AoA spread mean, in degree
        paraSt.sigma_phi_MS = 3.4; %AoA spread std, in degree
        paraSt.mu_theta_MS = 1.3; %EoA spread mean, in degree
        paraSt.sigma_theta_MS = 3.3; %EoA spread std, in degree
        
        paraSt.sigma_sh = 2.9; %shadowing[dB]
        
        %Autocorrelation distances, not used yet
        paraSt.l_s = 5;
        paraSt.l_tau = 5;
        paraSt.l_phi_BS = 50;
        paraSt.l_theta_BS = 50;
        paraSt.l_phi_MS = 25;
        paraSt.l_theta_MS = 25;
        
        %Cross-correlations
        paraSt.rho=[1    0.04 0   -0.2  0    0.2;...        % shadow fading
                     0.04 1    0    0.1  0    0.1;...        % delay spread
                     0    0    1    0    0    0;...        % theta BS
                    -0.2  0.1  0    1    0    0;...        % phi BS
                     0    0    0    0    1    0;...        % theta MS
                     0.2  0.1  0    0    0    1];          % phi MS   
        paraSt.corr_mat=chol(paraSt.rho);    % the correlation matrix is caculated using the Cholesky factorization        
        
        %Other parameters
        paraSt.sigma_phi_c = 60/180*pi; %angular deviation of cluster to visibility region 
        paraSt.sigma_r = 50; %radius deviation of cluster to V.R.
        paraSt.r_min = 400; %Cut-off distance for V.R. 
        paraSt.r_max = 400;
        
        paraSt.n_local_c = 1; %Number of local cluster
        
    %%%%%%%%%%%%%%%%    
    case 'picocell'
    %%%%%%%%%%%%%%%%
        paraSt.r_c = paraEx.net_radii; %Visibility region size, not available at the moment
        paraSt.l_c = paraEx.net_radii; %Transition region size, not available at the moment
        %We have to assume a visibility of all clusters.

        paraSt.k_tau_low = 50; %Cluster power
        paraSt.k_tau_high = 100; %Cluster power
        paraSt.k_tau = rand(1)*50+50;
        paraSt.tau_b = inf; %Cluster power, not available at the moment

        %Line of sight
%         paraSt.d_co = --; %LOS cut-off distance
%         paraSt.r_l = --;
%         paraSt.l_l = --;
%         paraSt.mu_k = --;
%         paraSt.delta_k = --;

        paraSt.n_c = 1; %Number of local clusters
        paraSt.k_sel = 0; %Only twin clusters    
        paraSt.n_c_add = 3; %average number of add. clusters(single interacting clusters)
        paraSt.n_mpc_min = 20; %Least number of MPCs per cluster
        paraSt.k_mpc = 0.5; %Rice factor of additional clusters
        
        paraSt.mu_diff = 0.05; %Diffuse radiation mean
        paraSt.delta_diff = 3.4; %Diffuse radiation std
        %paraSt.pdp_tau = --; %PDP of the diffuse radiation same as the corresponding clusters?
        
        %Delay spread
        paraSt.cluster_decay_low = 150; %Cluster decay, uniform distribution db/us
        paraSt.cluster_decay_high = 800; 
        
        paraSt.k_phi_BS_low = 30; %azimuth spread at BS, in degree, 5+k*exp(-(0.15*d)^2), k uniform dist.
        paraSt.k_phi_BS_high = 70; 
        paraSt.theta_BS_low = 5; %elevation spread at BS, in degree, uniform distributed
        paraSt.theta_BS_high = 10; 
        
        paraSt.phi_MS_low = 30; %azimuth spread at MS, mean, in degree, uniform dist.
        paraSt.phi_MS_high = 70;         
        paraSt.theta_MS_low = 5; %elevation spread at MS, mean, in degree, uniform dist.
        paraSt.theta_MS_high = 10; 
        
        paraSt.sigma_sh = 3; %shadowing, in dB
        
        %Autocorrelation distances, not available at the moment
%         paraSt.l_s = --;
%         paraSt.l_tau = --;
%         paraSt.l_phi_BS = --;
%         paraSt.l_theta_BS = --;
%         paraSt.l_phi_MS = --;
%         paraSt.l_theta_MS = --;
        
         %Cross-correlations, not available at the moment
        paraSt.rho=[1    0    0    0    0    0;...        % shadow fading
                     0    1    0    0    0    0;...        % delay spread
                     0    0    1    0    0    0;...        % theta BS
                     0    0    0    1    0    0;...        % phi BS
                     0    0    0    0    1    0;...        % theta MS
                     0    0    0    0    0    1];          % phi MS   
        paraSt.corr_mat=chol(paraSt.rho);    % the correlation matrix is caculated using the Cholesky factorization        
        
        %Other parameters
        paraSt.sigma_phi_c = 60/180*pi; %angular deviation of cluster to visibility region 
        paraSt.sigma_r = 50; %radius deviation of cluster to V.R.
        paraSt.r_min = 400; %Cut-off distance for V.R. 
        paraSt.r_max = 400;   
    %%%%%%%%%%%%%%%%        
    case 'aaltoOLOS' %Aalto Obstructive LOS measurements
    %%%%%%%%%%%%%%%%        
        paraSt.r_c = 2.72; %Visibility region size m
        paraSt.l_c = 1; %Transition region size m
        paraSt.k_tau = 16; %Cluster power decaying factor lin./us
        paraSt.tau_b = 30/3e8*1e6; %Cluster power decaying cut-off delay us (cell size)
        
        paraSt.d_co = 30; %LOS cut-off distance m
        paraSt.r_l = 0; %LOS VR size m (to be checked)
        paraSt.l_l = 0; %LOS TR size m
        paraSt.mu_k = -10; %LOS power factor, mean dB
        paraSt.sigma_k = 7.2; %LOS power factor, std dB
        
        paraSt.n_c_local = 2; %Average number of local cluster (1, local at MS, 2, local at BS and MS)
        paraSt.mu_n_c_far = 3.8; %Average number of far clusters, mean
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n_c_far = poissrnd(paraSt.mu_n_c_far); %Possion distributed number of clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n_c_far ==0
            n_c_far = 1;
        end
        paraSt.n_c_far = n_c_far; %Average number of far clusters
        
        paraSt.k_sel = 0; %K-selection factor (ratio of twin cluster)
        
        paraSt.mu_tau = 3e-9; %Delay spread mean s
        paraSt.sigma_tau = 0.6e-9; %Delay spread std s
        paraSt.mu_phi_BS = 3; %AoD spread mean deg.
        paraSt.sigma_phi_BS = 2; %AoD spread std deg.
        paraSt.mu_theta_BS = 1.5; %EoD spread mean deg.
        paraSt.sigma_theta_BS = 1.5; %EoD spread std deg.
        paraSt.mu_phi_MS = 4; %AoA spread mean deg.
        paraSt.sigma_phi_MS = 2; %AoA spread std deg.
        paraSt.mu_theta_MS = 1; %EoA spread mean deg.
        paraSt.sigma_theta_MS = 1; %EoA spread std deg.        
        paraSt.sigma_sf = 0; %shadow fading std dB
        
        paraSt.n_mpc = 3; %Number of MPCs per cluster
        
        %Cross-correlations
        paraSt.rho =[1    0    0    0    0    0;...        % shadow fading
                     0    1    0    0    0    0;...        % delay spread
                     0    0    1    0    0    0;...        % theta BS
                     0    0    0    1    0    0;...        % phi BS
                     0    0    0    0    1    0;...        % theta MS
                     0    0    0    0    0    1];          % phi MS   
        paraSt.corr_mat=chol(paraSt.rho);    %Cholesky factorization
               
        %Other parameters
		paraSt.mean_phi_c = 43; %azimuth of cluster to VR, mean deg.
        paraSt.sigma_phi_c = 39; %azimuth of cluster to VR, std deg.
        paraSt.sigma_theta_c = 0; %elevation of cluster to visibility region , std
        paraSt.min_r_c = 0; %Minimum distance from cluster to BS/MS
        paraSt.max_r_c = 30; %Distance from cluster to BS/MS, uniform distr. m
        paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10);%LOS power factor, dB?
        
        %Polarization
        paraSt.mu_xpdv = 0;
        paraSt.sigma_xpdv = 0; %Mean and std for XPDV
        paraSt.mu_xpdh = 0; 
        paraSt.sigma_xpdh = 0; %Mean and std for XPDH
        paraSt.mu_cpr =0;
        paraSt.sigma_cpr = 0; %Mean and std for CPR
        
        %For multi-link
        paraSt.ratio_common = 0.65; %Ratio of common cluster to total cluster
        paraSt.d_common = 2.72; %lifetime of common cluster, in distance (m)               
		
		%DMC		
    %%%%%%%%%%%%%%%%        
    case 'test' %Test scenarios
    %%%%%%%%%%%%%%%%        
        paraSt.r_c = 2.72; %Visibility region size [m]
        paraSt.l_c = 1; %Transition region size [m]
        paraSt.k_tau = 16; %Cluster power decaying factor lin. /us
        paraSt.tau_b = 30/3e8*1e6; %Cluster power decaying cut-off delay us (cell size)
        
        paraSt.d_co = 30; %LOS cut-off distance [m]
        paraSt.r_l = 0; %LOS VR size [m]
        paraSt.l_l = 0; %LOS TR size [m]
        paraSt.mu_k = 10^(-10/10); %LOS power factor, mean linear
        paraSt.sigma_k = 7.2; %LOS power factor, std dB
        paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); %LOS power factor        
                
        paraSt.BSLocal = true; %BS Local cluster activity true/false
        paraSt.MSLocal = true; %MS local cluster activity true/false
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        paraSt.mu_n_c_far = 3.8; %Average number of far clusters, mean
        n_c_far = poissrnd(paraSt.mu_n_c_far); %Possion distribution
        if n_c_far==0 n_c_far = 1; end %Guarantee min 1 far cluster
        paraSt.n_c_far = n_c_far; %Average number of far clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        paraSt.k_sel = 0; %K-selection factor (ratio of single cluster)
        
        E = 3e-9; %Delay spread mean [s]
        V = 0.6e-9; %Delay spread std [s]
        paraSt.mu_tau = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %Delay spread mean [s]
        paraSt.sigma_tau = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %Delay spread dB scale std        
        E=3; %AoD spread mean [deg]
        V=2; %AoD spread std [deg]		
        paraSt.mu_phi_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %AoD spread mean [deg]
        paraSt.sigma_phi_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %AoD spread dB scale std        
        E=1.5; %EoD spread mean [deg]
        V=1.5; %EoD spread std [deg]
        paraSt.mu_theta_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %EoD spread mean [deg]
        paraSt.sigma_theta_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %EoD dB scale spread
        E=4; %AoA spread mean [deg]
        V=2; %AoA spread std [deg]		
        paraSt.mu_phi_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %AoA spread mean [deg]      
        paraSt.sigma_phi_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %AoA dB scale spread std 
        E=1; %EoA spread mean [deg]
        V=1; %EoA spread mean [deg]
        paraSt.mu_theta_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %EoA spread mean [deg]
        paraSt.sigma_theta_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %EoA dB scale std       
        paraSt.sigma_sf = 0; %shadow fading std dB
        		
        paraSt.n_mpc = 3; %Number of MPCs per cluster
        
        %Cross-correlation matrix
        paraSt.rho =[1    0    0    0    0    0;...        % shadow fading
                     0    1    0    0    0    0;...        % delay spread
                     0    0    1    0    0    0;...        % theta BS
                     0    0    0    1    0    0;...        % phi BS
                     0    0    0    0    1    0;...        % theta MS
                     0    0    0    0    0    1];          % phi MS   
        paraSt.corr_mat=chol(paraSt.rho);    %Cholesky factorization
               
        %Cluster distribution & cluster link delay
		paraSt.phi_c = [43 39]; %mean/std azimuth of cluster to VR, mean [deg]        
		paraSt.pdf_phi_c = 'norm'; %phi_c normal distribution
		paraSt.theta_c = [0 0]; %mean/std elevation of cluster to visibility region [deg]    
		paraSt.pdf_theta_c = 'norm';  %theta_c normal distribution
        paraSt.para_r_c = [0 30]; %Min/max distance from cluster to BS/MS        
		paraSt.pdf_r_c = 'unif'; %r_c uniform distribution        
        paraSt.mu_tauCLink = 1e-9; %Cluster link delay mean [s]
        paraSt.sigma_tauCLink = 3.4; %Cluster link delay dB scale std        
        
        %Polarization
        paraSt.mu_xpdv = 0;
        paraSt.sigma_xpdv = 0; %Mean and std for XPDV
        paraSt.mu_xpdh = 0; 
        paraSt.sigma_xpdh = 0; %Mean and std for XPDH
        paraSt.mu_cpr =0;
        paraSt.sigma_cpr = 0; %Mean and std for CPR
        
        %For multi-link extension
        paraSt.BS_common = [0.35 0.65/2]; %BS-common cluster ratio to total number of clusters (1,2,3,...)
        paraSt.MS_common = 4; %Average number of VRs in VR group for one MS-common cluster
        
        %DMC parameters
        paraSt.dmc_spread_mean = 0; %dmc spatial spread [m]
        paraSt.dmc_spread_sigma = 0; %dmc spatial spread [m]
        paraSt.dmc_beta = 0; %dmc delay power decaying [dB]
end


end