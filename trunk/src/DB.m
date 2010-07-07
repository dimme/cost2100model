%DB function to transform the linear value into dB scale
%[ varout ] = DB( varin ) output the varin in dB scale. If varin is a
%complex number, then varout = 20*log10(abs(varin))+j*C*phase(varin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng ,Universit� catholique de Louvain, Belgium
%This file is part of cost2100_model.

%cost2100_model is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.

%cost2100_model is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.

%You should have received a copy of the GNU General Public License
%along with cost2100_model.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ varout ] = DB( varin )

varout = 20*log10(varin+eps);
%Real is the amplitude, imaginary is the phase
