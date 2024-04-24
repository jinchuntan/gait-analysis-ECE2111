% Assignment 1
% Written by Tan Jin Chun - 32194471
% Last Modified: 20/10/2021
% File Name: section32.m

clear all;clc;close all;

%% Section 3.2
% Part 1

% Loading in the data
load assignment_data.mat

% Extracting the data from the data set provided
time = C1(:,1);
left_foot = C1(:,2);
right_foot = C1(:,3);

% Initialising the variable
Samples = [10 20 50];
interval = 30;
Fs = 120;
str = "";

% Getting the start time and the time interval
start_time = (11*Fs) + 1;
end_time = start_time + interval*Fs;
indices = start_time:end_time;
time_interval = time(indices);
left_interval = left_foot(indices);
right_interval = right_foot(indices);

% Plotting the spectrogram of the VGRF from either the right or the left
% Plotting the spectrogram of the VGRF for the left foot (for this case)
% foot for any 30 seconds of data
% Using a for loop to plot for each of the different samples
for i = 1:length(Samples)
    
    % Getting the Sample
    Sample = Samples(i);
    
    % Creating a new figure
    figure;
    
    % Creating the spectrogram
    % Spectrogram function will compute an FFT-based spectral estimate over
    % each sliding window, allowing us to visualize how the frequency
    % content of the signal changes over time.
    % y-axis to display frequency on the vertical axis and time on the 
    % horizontal axis
    spectrogram(left_interval, Sample,0, Fs, Fs, "yaxis");
    hold on;
    
    % Calling the spectrogram function again
    % P- Power Spectral Density
    [S, W, T, P] = spectrogram(left_interval, Sample ,0, Fs, Fs, "yaxis");
    
    % Sum of the power density
    sum_pd = sum(P);
    
    % Getting the relevant values to calculate the threshold of the phase
    max_pd = max(sum_pd);
    min_pd = min(sum_pd);
    phase_threshold = (max_pd - min_pd) / 2;
 
    % Getting the vector for the swing time
    swing_time = T(phase_threshold > sum_pd);
    
    % Getting the vector for the stance time
    stance_time =  T(sum_pd > phase_threshold);
    
    % Adding the title to the figure
    str = sprintf("Spectrogram of Left VGRF signal with Window Length %d", Sample);
    title(str);
    
    % Plotting the stance time and the swing time on the figure itself
    plot(stance_time, 5, "bd", "MarkerSize", 2, "LineWidth", 1);
    plot(swing_time, 2, "ko", "MarkerSize", 2, "LineWidth", 1);
    
    % Using the colormap jet option for a better visualization
    colormap jet;
    
    % Adding a legend in the figure
    legend("Stance Phase" , "Swing Phase");
    
    % Holding off the figure
    hold off;
    
end

%% Questions
% Question 1
% From visual observation of the spectrogram that we have plotted above,
% we can observe that when the window length is increased, the spectrogram will
% become blurrer. This is because as the number of the samples increase, the resolution 
% will decrease as the frequency content is an accumulation of each sub-phases of the 
% stance and the swing phase. 

% Question 2
% The main noise component can be pinpointed to the window length of 50
% just by plotting the spectrogram and relying on visual observation alone.
% The power noise at around the 47-55 Hz (Mainly at 50Hz). One major source
% of noise would be the power line interference.

