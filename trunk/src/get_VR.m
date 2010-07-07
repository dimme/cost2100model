function VR = get_VR( VRtable, paraEx, paraSt)
%GET_VR function to generate the VR
%Default call: VR = get_VR( VRtable, paraEx, paraSt)
%------
%Input:
%------
%paraEx,paraSt: external parameters and stochastic parameters
%VRtable: VR assignment table
%------
%Output:
%------
%VR: the VR distribution (numVR, [x y])
%
%See also: cost2100, get_para, get_VRtable

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

tic %For generating state of rand state
nVR = length(VRtable(1,:,1)); %VR number
posBS = paraEx.pos_BS(:,1:2); %BS position (x,y) [m]
net_radii = paraEx.net_radii; %Radius of cell
rangeX = [min(posBS(:,1))-net_radii max(posBS(:,1))+net_radii]; %VR distribution range
rangeY = [min(posBS(:,2))-net_radii max(posBS(:,2))+net_radii]; %VR distribution range


x = (rangeX(2)-rangeX(1))*rand(nVR,1)+rangeX(1); %Uniform distribution
y = (rangeY(2)-rangeY(1))*rand(nVR,1)+rangeY(1); %Uniform distribution
VR(:,1) = x;
VR(:,2) = y;

% for m = 1:nVR
%     S = toc*1e6;
%     rand('twister',S);   
%     x = (rangeX(2)-rangeX(1))*rand+rangeX(1); %Uniform distribution
%     y = (rangeY(2)-rangeY(1))*rand+rangeY(1); %Uniform distribution
%     while any(sqrt((posBS(:,1)-x).^2+(posBS(:,2)-y).^2)> paraEx.net_radii)
% %         S = toc*1e6;
% %         rand('twister',S);
%         x = (rangeX(2)-rangeX(1))*rand+rangeX(1); %Uniform distribution
%         y = (rangeY(2)-rangeY(1))*rand+rangeY(1); %Uniform distribution
%     end
%     VR(m,1) = x;
%     VR(m,2) = y;
% end

