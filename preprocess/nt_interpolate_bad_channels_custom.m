function [toGood,fromGood]=nt_interpolate_bad_channels_custom(x,iBad,coordinates,n)
%y=interpolate_bad_channels(x,iBad,coordinates,n) - interpolate bad channels from good
%
%  y: interpolated data
% 
%  x: data to interpolate
%  iBad: indices of bad channels
%  coordinates: coordinate map (see nt_proximity)
%  n: number of neighboring channels to use [default: 3]
%
% NoiseTools;


if nargin<3; 
    error('Input arg number < 3'); 
end
if nargin<4; 
    n=3;
end

nchans=size(x,2);
toGood=eye(nchans);
toGood(:,iBad)=[];

[closest,d]=nt_proximity_custom(coordinates);
if size(closest,1)~=nchans; error('Length not equal with data nchan and given coordinate file.'); end

fromGood=eye(nchans);
for iChan=iBad
    iOthers=closest(iChan,:);
    iOthers=setdiff(iOthers, iBad, 'stable'); % don't include bad channels
    if numel(iOthers)<n; warning('Not enough good neighbors to interpolate bad channels!'); end
    iOthers=iOthers(1:n);
    w=1./(d(iChan,iOthers) + eps);
    w=w/sum(w);
    fromGood(iOthers,iChan)=w;
end
fromGood(iBad,:)=[];
   

