function [paraEx paraSt link env BS MS ] = cost2100(network, freq, snapRate, snapNum, BSPos,MSPos,MSVelo, antBS, antMS)

%COST2100 The main body of COST 2100 channel model implementation,
%with multi-link extension implemented
%
%Default call:
%[paraEx paraSt link env BS MS ] = cost2100(network, freq, snapRate,
%snapNum, BSPos,MSPos,MSVelo)
%
%------
%Input:
%------
%network: network type, fours scenarios are available: 
%'macrocell', 'microcell', 'picocell', 'aalto'
%freq: frequency range, 2x1 vector, [freq_start freq_stop], Hz
%snapRate: channel snapshot rate, s
%snapeNum: channel snapshot number
%BSPos: position matrix of BSs, size = (number of BS, [x y z]), m
%MSPos: initial position matrix of MSs, size = (number of MS, [x y z]), m
%MSVelo: velocity vector of MSs, size = (number of MS, [x y z]), m/s
%
%------
%Output:
%------
%paraEx: external parameters
%paraSt: stochastic parameters
%BS(numBS): BS information array: 
% .idx: label of BS  (1,2,...)
% .pos: position of each BS [x y z], m
% .VRLOS: LOS VR for each BS
% .cluster_local: local clusters surrounding each BS
% .mpc_local: local mpcs surrounding each BS
% .dmc_local: local dmcs in BS local cluster
%MS(numMS): Array of Information related with MSs: 
% .idx: label of MS  (1,2,...)
% .pos: position of each MS, [x y z], m
% .velo: velocity vector of each MS, [x y z], m
% .cluster_local: local cluster surrounding each MS
% .mpc_local: local mpc surrounding each MS
% .dmc_local: local dmcs in BS local cluster
%link(numBS,numMS): Link information matrix
% .channel: the channel information
% .IR: impulse response of channel
% .MSRec: record of MS information at each snapshot
%env: radio propagation environment record
% .VRtable: VR assignment table 
% .VR: 
% .cluster:
% .mpc: 
% .dmc: 
%
%See also: get_para, get_VRtable, get_VR, get_cluster, get_mpc, get_dmc,
%get_VRLOS, get_cluster_local, get_channel, get_IR, update_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Lingfeng, ICTEAM, UCL, Belgium 
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
%
%--------------------------------------------------------------
%Author: LIU Lingfeng
%Organization: COST2100 & University catholique de Louvain, Belgium
%Address: Bâtiment Maxwell, Place du Levant 3, BE-1348, Louvain-la-Neuve
%Create date: 20-10-08
%Last Modified: 09-04-10
%Version histories: 
%                  1.00, 25/05/2005 by HH (ftw)                    
%                  1.10, 08/07/2006 by HH (Eurecom)
%                  1.2.0, 22/09/2008 by LIU (UCL) 
%                  1.2.1, 29/10/2008 by LIU (UCL)
%                  1.3    19/01/2010 by LIU (UCL)
%                  2.1    04/02/2010 by LIU (UCL)   
%                  2.2    08/04/2010 by LIU (UCL)   
%--------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(freq)~=2
    error('Please specify freq as [freq_start freq_stop].\n');
end
chkVal = size(BSPos);
if chkVal(2)~=3
    error('Please specify each BS position as [x y z].\n');
end
chkVal = size(MSPos);
if chkVal(2)~=3
    error('Please specify each MS position as [x y z].\n');
end 
if ~all(size(MSPos)==size(MSVelo))
    error('Please check MSPos and MSVelo.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numBS = length(BSPos(:,1)); %Number of BSs
numMS = length(MSPos(:,1)); %Number of MSs

[paraEx paraSt] = get_para(network, freq, snapRate, snapNum, BSPos, MSPos, MSVelo, antBS, antMS);
VRtable = get_VRtable(paraEx,paraSt);
VR = get_VR(VRtable, paraEx, paraSt);
cluster = get_cluster(VR, VRtable, paraEx, paraSt);
mpc = get_mpc(cluster, paraSt);
dmc = get_dmc(cluster, paraSt);

for m = 1:numBS    
    BS(m).idx = m; %Label of BS(m)
    BS(m).pos = BSPos(m,:); %POsition of BS
    BS(m).VRLOS = get_VRLOS(BS(m), paraEx, paraSt); %LOS VR of BS
    BS(m).cluster_local = get_cluster_local('BS', m, paraEx, paraSt); %BS local cluster
    BS(m).mpc_local = get_mpc( BS(m).cluster_local, paraSt); %MPC in BS local cluster
    BS(m).dmc_local = get_dmc( BS(m).cluster_local, paraSt); %DMC in BS local cluster         
end

for m = 1:numMS
    MS(m).idx = m; %Label of MS(m)       
    MS(m).pos = MSPos(m,:); %Position of MS
    MS(m).velo = MSVelo(m,:); %Velocity of MS
    
    MS(m).cluster_local = get_cluster_local('MS', m, paraEx, paraSt); %MS local cluster
    MS(m).mpc_local = get_mpc( MS(m).cluster_local, paraSt); %MPC in MS local cluster
    MS(m).dmc_local = get_dmc( MS(m).cluster_local, paraSt); %DMC in MS local cluster          
end 


for m = 1:snapNum
    %Get the channel of each link
    for nB = 1:numBS
        for nM = 1:numMS
            link(nB,nM).channel{m} = get_channel(BS(nB),MS(nM), VRtable, VR, cluster, mpc, dmc, paraEx, paraSt); %Get the unsampled channels            
            link(nB,nM).MS(m) = MS(nM); %Record of MS at snapshot m
            link(nB,nM).BS(m) = BS(nB); %Record of BS at snapshot m            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Get the channel IR covariant matrix
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load(antBS(nB).antFile);
            txFull=antFull;
            txRot=[antBS(nB).aziRot antBS(nB).eleRot];
            load(antMS(nM).antFile);
            rxFull=antFull;           
            rxRot=[antMS(nM).aziRot antMS(nM).eleRot];
            link(nB,nM).H{m} = get_H(link(nB,nM).channel{m},txFull,rxFull,txRot,rxRot,paraEx);                        
                        
            %link(nB,nM).IR{m} = get_IR(link(nB,nM).channel{m},paraEx,paraSt);  %impulse response
        end    
    end
    
    %Update the MS
    for nM = 1:numMS
        MS(nM) = update_chan( MS(nM), paraEx, paraSt );
    end
end

%Saving the environment information
env.VRtable = VRtable;
env.VR = VR;
env.cluster = cluster;
env.mpc = mpc;
env.dmc = dmc;