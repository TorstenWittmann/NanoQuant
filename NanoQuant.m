% NanoQuant loads an NIS Elements stack to measure sub-resolution nanocage fluorescence intensities. 
% Requires MATLAB Image processing and Statistics and Machine learning toolboxes
% and the Bio-Formats toolbox to load images (https://www.openmicroscopy.org/bio-formats/)

% Buttons on the figure allow user to step through the image stack
% Slider allows to adjust display brightness; default set to 20% of image maximum
% After selecting an ROI, local maxima are fitted with a 2D Gaussian 
% to measure the fluorescence intensity in a circle with a radius of two
% times the standard deviation, which encompasses ~86% of the signal, in both the original image and the fit.
% A new window shows the fit and if accepted results are appended to the
% results table, and the fit is overlaid on the data window. Fitting results are in pixel and intensity.
% Torsten Wittmann. University of California, San Francisco. 2024.

% DISCLAIMER:
% This software is provided "as-is," without any express or implied warranty. 
% In no event will the authors be held liable for any damages arising from 
% the use of this software.

% Permission is granted to anyone to use this software for any purpose.
% Redistributions of source code must retain the this notice and disclaimer.
% Altered versions must be plainly marked as such, and must not be 
% misrepresented as being the original software.

close all
clearvars

% Load the image file
[filename, filepath] = uigetfile('*.nd2');      % Let the user choose a .nd2 file
fullpath = fullfile(filepath, filename);
data = bfopen(fullpath);                        % Open the .nd2 file using Bio-Formats
imgStack = data{1};                             % This is the stack of images
img_raw = imgStack{1};                          % This is the first image in the stack

pixel_size = 27;                                % Image pixel size in nanometer, adjust as needed
rayleigh = 0.61 * 520 / 1.49;                   % Rayleigh resolution
sigma_guess = 0.5 * rayleigh / pixel_size;      % Initial guess of PSF Gaussian fit standard deviation

% Initialize results table
results = table([], [], [], [], [], [], [], 'VariableNames', {'CenterX', 'CenterY', 'Amplitude', 'Sigma', 'Offset', 'TotalIntensityImage', 'TotalIntensityFit'});

% Create figure and axis for the image
f = figure('Name','Image Stack');
set(f, 'WindowButtonMotionFcn', @(~,~)[]);
ax1 = axes(f);
max_val = double(max(img_raw(:)));
min_val = double(min(img_raw(:)));
imh = imshow(img_raw, 'Parent', ax1, 'DisplayRange', [min_val 0.2 * max_val]);
zoom on;
set(f, 'UserData', 1);                          % First image in stack is displayed

% Define slider; set to 20% of max image intensity
slider = uicontrol('Parent', f, 'Style', 'slider', 'Position', [10, 10, 200, 20],...
                   'value', 0.2 * max_val, 'min', min_val, 'max', max_val);

% Set the callback for the slider
set(slider, 'Callback', @(src, ~) adjustDisplay(src, imh, min_val));

% Add buttons for scrolling through the image stack:
btnPrev = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', 'Previous image',...
    'Position', [330 10 100 20], 'Callback', @(src, ~) prevImage(src, f, imh, imgStack));
btnNext = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', 'Next image',...
    'Position', [440 10 100 20], 'Callback', @(src, ~) nextImage(src, f, imh, imgStack));

% Button to perform fit at selected point
btn = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', 'Perform Fit',...
                'Position', [220 10 100 20], 'Callback', @(src, ~) performFit(src, f, imh, ax1, sigma_guess));

% Function to adjust display brightness
function adjustDisplay(src, imh, min_v)
    cmax = src.Value;
    clim(imh.Parent, [min_v - 1 cmax]); 
end

% Functions for the previous and next buttons
function prevImage(~, f, imh, imgStack)
    currImageIndex = get(f, 'UserData');        % Get the current image index from the figure UserData property
    if currImageIndex > 1                       % If the current image is not the first one in the stack
        currImageIndex = currImageIndex - 1;    % Go to the previous image
    end
    set(f, 'UserData', currImageIndex);         % Save the new image index in the UserData property
    img = imgStack{currImageIndex};             % Get the new image
    set(imh, 'CData', img);                     % Update the displayed image
end

function nextImage(~, f, imh, imgStack)
    currImageIndex = get(f, 'UserData'); 
    if currImageIndex < length(imgStack)        % If the current image is not the last one in the stack
        currImageIndex = currImageIndex + 1;    % Go to the next image
    end
    set(f, 'UserData', currImageIndex); 
    img = imgStack{currImageIndex};  
    set(imh, 'CData', img); 
end

% Function to perform 2D Circular Gaussian fit
function performFit(~, f, imh, ax, sg)
    set(f, 'Pointer', 'crosshair');            % Change the cursor to a crosshair
    w = round(5 * sg);                         % Window size for subsampling for Gaussian fitting
    waitforbuttonpress;                        % Wait for the user to click on the image
    point = get(ax, 'CurrentPoint');           % Get click coordinates   
    xClick = round(point(1, 1));
    yClick = round(point(1, 2));

    img_raw = get(imh, 'CData');               % Get the current image from the image handle
    img = double(img_raw);

    wr = w / 2;                                % Extract window around the clicked point
    xMin = max(1, xClick - wr);
    xMax = min(size(img, 2), xClick + wr);
    yMin = max(1, yClick - wr);
    yMax = min(size(img, 1), yClick + wr);
    imgWindow = img(yMin:yMax, xMin:xMax);

    % Create coordinate arrays for the window
    [X_window, Y_window] = meshgrid(1:size(imgWindow, 2), 1:size(imgWindow, 1));
    X_flat = X_window(:);                      % Flatten arrays for fitting function
    Y_flat = Y_window(:);
    img_flat = imgWindow(:);

    % Define 2D Gaussian fitting function with offset
    % Parameters: [Amplitude, x0, sigma, y0, offset]
    fitfcn = @(p, xy) p(1) * exp(-((xy(:,1) - p(2)).^2 + (xy(:,2) - p(4)).^2) / (2 * p(3)^2)) + p(5);

    % Define initial guess ([Amplitude, x0, sigma, y0, offset])
    [amplitude_guess, idx] = max(img_flat);    % Find the maximum pixel value and its position
    [y_guess, x_guess] = ind2sub(size(imgWindow), idx);     % Center window around maximum
    initial_guess = [amplitude_guess, x_guess, sg, y_guess, 0.5 * min(img_flat)];

    % Perform fit
    fit_model = fitnlm([X_flat, Y_flat], img_flat, fitfcn, initial_guess);

    f3 = figure;                               % Create a new figure for 3D data visualization

    % Calculate fitted Gaussian data
    fit_data = reshape(fitfcn(fit_model.Coefficients.Estimate, [X_flat, Y_flat]), size(X_window));
    
    
    hold on;                                   % Plot original image intensities and Gaussian fit
    scatter3(X_flat, Y_flat, img_flat, 'b');   % Original image data as points
    mesh(X_window, Y_window, fit_data);        % Gaussian fit as mesh
    xlim([1 w]);                               % Set x and y limits
    ylim([1 w]);
    hold off;
    view(3);
    title('3D plot of Original Data and Gaussian Fit');

    % Define "Accept" button
    uicontrol('Parent', f3, 'Style', 'pushbutton', 'String', 'Accept',...
    'Position', [20 20 60 20], 'Callback', @(src, ~) acceptCallback(src, f3, fit_model, imgWindow, ax, xMin, yMin, X_window, Y_window, fit_data));
    % Define "Cancel" button
    uicontrol('Parent', f3, 'Style', 'pushbutton', 'String', 'Cancel',...
    'Position', [90 20 60 20], 'Callback', 'close(gcbf)');

function acceptCallback(~, f3, fit_model, imgWindow, ax, xMin, yMin, X_window, Y_window, fit_data)
    % Calculate center coordinates with respect to the original image
    centerX = fit_model.Coefficients.Estimate(2) + xMin - 1;
    centerY = fit_model.Coefficients.Estimate(4) + yMin - 1;
    
    % Draw the fitted center on the original image
    hold(ax, 'on');
    plot(ax, centerX, centerY, 'r.','MarkerSize',20);
    
    % Draw a circle with a radius of 2x sigma on the original image
    theta = linspace(0, 2*pi, 100);
    r = 2*fit_model.Coefficients.Estimate(3);  % radius
    x_circle = r*cos(theta) + centerX;
    y_circle = r*sin(theta) + centerY;
    plot(ax, x_circle, y_circle, 'r', 'LineWidth', 1.5);
    hold(ax, 'off');

    % Create mask for the 2x sigma circle
    mask = (X_window - fit_model.Coefficients.Estimate(2)).^2 + (Y_window - fit_model.Coefficients.Estimate(4)).^2 <= r^2;

    % Compute total intensity in the 2x sigma circle from the original image data and fit data, remove offset
    total_intensity_image = sum(imgWindow(mask)) - fit_model.Coefficients.Estimate(5)*sum(mask(:));
    total_intensity_fit = sum(fit_data(mask)) - fit_model.Coefficients.Estimate(5)*sum(mask(:));

    % Add the results to a table
    t = table(centerX, centerY, fit_model.Coefficients.Estimate(1), fit_model.Coefficients.Estimate(3), fit_model.Coefficients.Estimate(5), total_intensity_image, total_intensity_fit,...
        'VariableNames', {'CenterX', 'CenterY', 'Amplitude', 'Sigma', 'Offset', 'TotalIntensityImage', 'TotalIntensityFit'});
            
    % Assign the updated results table to a variable in the base workspace
    % Append t to the results table
    tempresults = evalin('base', 'results');
    tempresults = [tempresults; t];
    assignin('base', 'results', tempresults);
    disp(t);
    close(f3);
end
    
    set(f, 'Pointer', 'arrow');             % Restore the cursor
end


