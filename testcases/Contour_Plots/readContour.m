function [u,v] = readContour(image, cmap, range, duct_radius, axialVel)

% Read in contour plot image
RGB = imread(image);

% Get pixel density
imgSize = size(RGB);

% Make copy
RGB_new = RGB;

% Array to store the indices of the pixels within the plot
plotPixels = [];
pix = 1;
count = 1;

% For tracking progress
totalPix = size(RGB_new,1)*size(RGB_new,2);
fprintf('Pixels in image: %d\n', totalPix)

% Loop through every pixel
for i = 1:imgSize(1)
    for j = 1:imgSize(2)
        % Extract pixel rgb
        color = squeeze(RGB(i,j,:));
         
        % Check if the pixel colour is grey scale - ie R,G,B array elements
        % are equal
        inPlot = true;
        if all(color==color(1))
            % If yes, turn the pixel black
            RGB_new(i,j,:) = [0,0,0];
            inPlot = false;
        end

        if inPlot
            plotPixels(count,:) = [i,j];
            count = count + 1;
        end
        
        pix = pix + 1;
        progress = pix/totalPix;
            
        % Progress reports
        if mod(progress,0.25) == 0
            fprintf('%.0f percent done\n', progress*100)
        end
            
    end
end

% Convert rgb to index
swirlangle = rgb2ind(RGB_new,cmap);

% Normalise and map to input range
swirlangle = double(swirlangle) / 255;
swirlangle = swirlangle*(range(2)-range(1)) + range(1);

% Make an array of pixel x,y coords (meshgrids) - origin in top left
[X,Y] = meshgrid(0:size(RGB,2)-1,0:size(RGB,1)-1);

% Convert 2d indices of pixels within plot to linear indices
plotPixels = sub2ind(size(RGB,1,2),plotPixels(:,1),plotPixels(:,2));

% Flatten coordinate and swirl angle arrays
x = reshape(X,[],1);
y = reshape(Y,[],1);
swirlangle = reshape(swirlangle,[],1);

% Discard the pixels which are outside the contour plot area
x = x(plotPixels);
y = y(plotPixels);
swirlangle = swirlangle(plotPixels);

% Get x,y coord of middle of contour plot
xc = (max(x)-min(x))/2 + min(x);
yc = (max(y)-min(y))/2 + min(y);

% Transform x,y coords to place origin in the middle of the contour plot
x = x-xc;
y = y-yc;

% Normalise and map to actual values using input duct radius
xmin = min(x);
xmax = max(x);
x = (x-xmin) ./ (xmax-xmin);
x = x .* 2*duct_radius - duct_radius;
ymin = min(y);
ymax = max(y);
y = (y-ymin) ./ (ymax-ymin);
y = y .* 2*duct_radius - duct_radius;

% Get polar coordinates relative to the centre of the contour plot
[theta,r] = cart2pol(x,y);

% Get tangential velocity
vel_theta = tan(deg2rad(swirlangle))*axialVel;

% Get theta dot
theta_dot = vel_theta./r;

% Calculate velocity in x and y direction - equations from solving a system
% of 3 equations using matlab symbolic toolbox
u = (r.*(theta_dot.*x.*cos(theta).^2 - r.*theta_dot.*cos(theta) + theta_dot.*x.*sin(theta).^2))./(y.*cos(theta) - x.*sin(theta));
v = (r.*(theta_dot.*y.*cos(theta).^2 - r.*theta_dot.*sin(theta) + theta_dot.*y.*sin(theta).^2))./(y.*cos(theta) - x.*sin(theta));

% Get dx and dy - uniformly spaced since pixel grid
dx = x(2)-x(1);
% Create new coord meshgrids for plotting
[X,Y] = meshgrid(min(x):dx:max(x),min(y):dx:max(y));

u(u==Inf) = 0;
v(v==Inf) = 0;

% Expand back variables of interest to grid so we can plot with the contourf function
swirl_grid = griddata(x,y,swirlangle,X,Y);
u_grid = griddata(x,y,u,X,Y);
v_grid = griddata(x,y,v,X,Y);

% Plot extracted swirl angle as a contour plot to see if it looks the same as the
% original image
figure()
contourf(X,Y,swirl_grid,256,'LineStyle','none')
colormap(cmap)
colorbar('southoutside')

% Plot streamlines
figure()
N = 1000;
xstart = max(x)*rand(N,1);
ystart = max(y)*rand(N,1);
streamline(X,Y,u_grid,v_grid,xstart,ystart)

end


