% define useful constants
wavelen_purple = 405;
wavelen_red = 650;
trim_wavelen_min = 420; % minimum wavelength to trim
trim_wavelen_max = 680; % maximum wavelength to trim
window_size = 200; % window size for smoothing

% demo data is used for red and purple laser
laser_intensity_red = csvread("data/demo_data/RedLaser.csv");
laser_intensity_purple = csvread("data/demo_data/PurpleLaser.csv");

% loading the collected led_intensity_white Data
led_intensity_white = csvread("data/blank-curvette.csv");

% find the peak locations for both red and purple lasers
[peak_red, peak_loc_red] = find_max_peak_loc(laser_intensity_red);
[peak_purple, peak_loc_purple] = find_max_peak_loc(laser_intensity_purple);

% Use mldivide to obtain grad (m) & offset (c) for the calibration curve
% Purple: x1 = peak_loc_purple; y1 = wavelen_purple
% Red: x2 = peak_loc_red; y2 = wavelen_red

A = [peak_loc_purple, 1; 
     peak_loc_red, 1];
b = [wavelen_purple; 
     wavelen_red];

% Assuming the matrix equation Ax = c were x = [m, c], we can use mldivide
x = mldivide(A, b);
m = x(1);
c = x(2);

% Calculate the wavelengths of different Lasers and WhiteLED
wv_white_led = (1:length(led_intensity_white)) * m + c;
wv_purple_laser = (1:length(laser_intensity_purple)) * m + c;
wv_red_laser = (1:length(laser_intensity_red)) * m + c;

% Plot all the absorbance curves
figure(1) 
plot(wv_white_led, led_intensity_white, 'c');
hold on; % hold on used to plot all in the same figure
plot(wv_purple_laser, laser_intensity_purple, 'b');
plot(wv_red_laser, laser_intensity_red, 'r');
grid on; 
legend('White LED', 'Purple Laser', 'Red Laser');
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');

% load blanks and dyes
blank_red = csvread("data/blank.csv");
blank_blue = csvread("data/blank2.csv");
dye_red = csvread("data/red.csv");
dye_blue = csvread("data/blue.csv");
   
% Calculate wavelengths and absorbance
absorbance_red = 1 - dye_red./blank_red;
absorbance_blue = 1 - dye_blue./blank_blue;
wv_red_cuvette = (1:length(dye_red)) * m + c;
wv_blue_cuvette = (1:length(dye_blue)) * m + c;

% trim the excessive noise 
[wv_blue_cuvette, absorbance_blue] = trim_spectra(wv_blue_cuvette, absorbance_blue, 420, 680);
[wv_red_cuvette, absorbance_red] = trim_spectra(wv_red_cuvette, absorbance_red, 420, 680);

% Smooth spectra using custom smoothing function

absorbance_blue = smooth_spectra(absorbance_blue, window_size);
absorbance_red = smooth_spectra(absorbance_red, window_size);

% Plot the absorbance values for both red and blue curvettes
figure(2)
plot(wv_blue_cuvette, absorbance_blue, 'b', 'LineWidth', 1.2);
hold on;
plot(wv_red_cuvette, absorbance_red, 'r', 'LineWidth', 1.2);
grid on; 
legend('Blue Cuvette', 'Red Cuvette');
xlabel('Wavelength (nm)');
ylabel('Absorbance (arbitrary units)');

% Function to find the location of the highest peak
function [max_peak, peak_location] = find_max_peak_loc(laser_spectra)
    % Find peaks
    [peaks, locations] = findpeaks(laser_spectra);
    % locate highest peak
    [max_peak, idx] = max(peaks);
    % Obtain x-axis value of the peak
    peak_location = locations(idx);
end

function [trimmed_wavelen, trimmed_absorbance] = trim_spectra(wavelen, absorbance, min_wavelen, max_wavelen)
    % Obtain indices to trim
    abs_diff_min = abs(wavelen - min_wavelen);
    [~, idx_1] = min(abs_diff_min);
    abs_diff_max = abs(wavelen - max_wavelen);
    [~, idx_2] = min(abs_diff_max);
    
    % Find which index is bigger / smaller as the wavelength might be
    % reversed
    min_index = min(idx_1, idx_2);
    max_index = max(idx_1, idx_2);
    
    trimmed_wavelen = wavelen(min_index:max_index);
    trimmed_absorbance = absorbance(min_index:max_index);
end

function smoothed = smooth_spectra(spectra, window_size)
    % Initialize an empty array and store len
    smoothed = zeros(size(spectra));
    len_spectra = length(spectra);
    
    % Compute the moving average
    for i = 1:length(spectra)
        % Determine the start and end indices for the window
        start_index = max(1, i - floor(window_size));
        end_index = min(len_spectra, i + floor(window_size));
        
        % Compute the average of the window
        smoothed(i) = mean(spectra(start_index:end_index));
    end
end