function [ dx, dy ] = latlon2xy( lats, lons, reflat, reflon, varargin )
%latlon2xy Converts a set of latitudes and longitudes to x-y coordinates
%   
% only want 3 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:latlon2xy:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = {true};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[useWGS] = optargs{:};

% now calculate dx and dy

if useWGS
    dy = distance(lats,reflon, reflat, reflon, referenceEllipsoid('wgs84'))/1000; 
    dx = distance(reflat,lons, reflat, reflon, referenceEllipsoid('wgs84'))/1000; 
else
    dy = deg2km(distance(lats,reflon, reflat, reflon));
    dx = deg2km(distance(reflat,lons, reflat, reflon));  
end
% dx and dy are positive only fix sign
i_lat = (lats<reflat)*(-1) + (lats>=reflat*(1)) ;% setm y distance to negative for lats < reflat
dy = dy.*i_lat;

i_lon = (lons<reflon)*(-1) + (lons>=reflon*(1)) ;% setm x distance to negative for lons < reflon
dx = dx.*i_lon;

end

