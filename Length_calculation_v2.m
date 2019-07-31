%
% NAME:
%               Length_calculation
%
% PURPOSE:
%               To calculate the pixel length of specific ROIs with thinning algorithm.
%               
%               Require Single_particle_analysis_v2.m or Particle_measurement_v2.m, or Particle_measurement_AICL_v2.m
%               
% 
%               Written by Dr Aleks Ponjavic and Jason C Sang, University of Cambridge, 
%               2015-2016
%
%               Fixed a bug in the loop 26-28
%               Last updated on 2017/04/08
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l,skel] = Length_calculation_v2(R, m, n, mag)

xs = R(:,2);
ys = R(:,1);
rr = zeros(m*mag, n*mag);
%rr = zeros(ceil(max(max(R))));

for i = 1:length(ys)
    rr(ys(i),xs(i)) = 1;
end


skel = bwmorph(rr,'thin',inf);

measurements = regionprops(skel, 'Perimeter');
l(:,1) = [measurements.Perimeter];

zero = find(l(:,1) == 0);
l(zero,1)=1;

end