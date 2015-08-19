

% Monitor size and position variables
w = 50;  % width of screen, in cm
h = 30;  % height of screen, in cm
cx = w/2;   % eye x location, in cm
cy = h/2; % eye y location, in cm

% Distance to bottom of screen, along the horizontal eye line
zdistBottom = dist;     % in cm
zdistTop    = dist;     % in cm

% Alternatively, you can specify the angle of the screen
%screenAngle = 72.5;   % in degrees, measured from table surface in front of screen to plane of screen
%zdistTop = zdistBottom - (h*sin(deg2rad(90-screenAngle)));

pxXmax = 128; % number of pixels in an image that fills the whole screen, x
pxYmax = 72; % number of pixels in an image that fills the whole screen, y

% Internal conversions
top = h-cy;
bottom = -cy;
right = cx;
left = cx - w;

% Convert Cartesian to spherical coord
% In image space, x and y are width and height of monitor and z is the
% distance from the eye. I want Theta to correspond to azimuth and Phi to
% correspond to elevation, but these are measured from the x-axis and x-y
% plane, respectively. So I need to exchange the axes this way, prior to
% converting to spherical coordinates:
% orig (image) -> for conversion to spherical coords
% Z -> X
% X -> Y
% Y -> Z

[xi,yi] = meshgrid(1:pxXmax,1:pxYmax);
cart_pointsX = left + (w/pxXmax).*xi;
cart_pointsY = top - (h/pxYmax).*yi;
cart_pointsZ = zdistTop + ((zdistBottom-zdistTop)/pxYmax).*yi;
[sphr_pointsTh sphr_pointsPh sphr_pointsR] ...
            = cart2sph(cart_pointsZ,cart_pointsX,cart_pointsY);

% view results
figure
subplot(3,2,1)
imagesc(cart_pointsX)
colorbar
title('image/cart coords, x')
subplot(3,2,3)
imagesc(cart_pointsY)
colorbar
title('image/cart coords, y')
subplot(3,2,5)
imagesc(cart_pointsZ)
colorbar
title('image/cart coords, z')

subplot(3,2,2)
imagesc(rad2deg(sphr_pointsTh))
colorbar
title('mouse/sph coords, theta')
subplot(3,2,4)
imagesc(rad2deg(sphr_pointsPh))
colorbar
title('mouse/sph coords, phi')
subplot(3,2,6)
imagesc(sphr_pointsR)
colorbar
title('mouse/sph coords, radius')


%% try a distortion

% make source image
checkSize = 4; % pixels per side of each check
w = 64; % width, in pixels
h = 36; % height, in pixels


% alternate source image
%I = zeros(150*4,200*4);
%I(105*4:125*4,:)=0.2;
%I(20*4:40*4,:)=0.4;

% Rescale the Cartesian maps into dimensions of radians
ymaxRad = max(sphr_pointsPh(:));
xmaxRad = (w/h)*ymaxRad;


fx = xmaxRad/max(cart_pointsX(:));
fy = ymaxRad/max(cart_pointsY(:));

% Apply the distortion via interpolation

moviewarp = zeros(size(moviedata));
display('warping')
tic
for f = 1:length(moviedata)
if f/600 ==round(f/600);
    sprintf('%f done',f/length(moviedata))
end

moviewarp(:,:,f) = interp2(cart_pointsX.*fx,cart_pointsY.*fy,squeeze(double(moviedata(:,:,f)))',sphr_pointsTh,sphr_pointsPh)';

end
toc

figure
imagesc(max(moviewarp,[],3)');axis equal; colormap gray

moviedata = uint8(moviewarp);
