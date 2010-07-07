clear all; clc;

network = 'test';
freq = [-60e6 60e6]+5.3e9; %Hz
snapRate = 1; %/s
snapNum = 10;
BSPos  = [0    0   0;
		  10   10  0;];
MSPos  = [0    10  0;
          10   5   0];
		  
MSVelo = [1    0   0;
          0    1   0];
numBS=length(BSPos(:,1));
numMS=length(MSPos(:,1));

antSample=struct('antFile','antSample','aziRot',0,'eleRot',0);
antBS=repmat(antSample,1,numBS); %BS antenna array information
antMS=repmat(antSample,1,numMS); %MS antenna array information

[paraEx paraSt link env BS MS ] = cost2100(network, freq, snapRate, snapNum, BSPos,MSPos,MSVelo,antBS,antMS);

%visual_channel(link(1,1),env,paraEx,paraSt)
%visual_pddp(link(1,1).channel,paraEx,paraSt)