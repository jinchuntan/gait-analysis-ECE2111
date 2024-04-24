% Assignment 1
% Written by Tan Jin Chun - 32194471
% Last Modified: 20/10/2021
% File Name: section31.m

clear all;clc;close all;

%% Section 3.1
% Loading in the data
% Name of the data set is C1
load assignment_data.mat

% Part 1
% Plotting the VGRF (in Newtons) vs time (in seconds), for both feet
% Extracting the data from the data set provided
time = C1(:,1);
left_foot = C1(:,2);
right_foot = C1(:,3);

% Initialising the variable
interval = 5;
Fs = 120;

% Starting time at 11
start_time = (11*Fs) + 1;

% End time
end_time = start_time + interval*Fs;
indices = start_time:end_time;
time_interval = time(indices);
left_interval = left_foot(indices);
right_interval = right_foot(indices);

% Normal Plotting without taking into consideration the time interval
% Plotting for the VGRF left foot vs time graph
figure;
plot(time, left_foot);
xlabel("Time (s)");
ylabel("VGRF Data (N)");
title("VGRF Data (N) vs Time (s) - Without Interval");

% Plotting for the VGRF right foot vs time graph (For observation purpose)
hold on
plot(time, right_foot);
legend("Left Foot", "Right Foot");
xlim([1 5]);
hold off

% Plotting the graph that take into account the time interval
% Plotting for the VGRF left foot vs time graph
figure;
plot(time_interval, left_interval);
xlabel("Time (s)");
ylabel("VGRF Data (N)");
title("VGRF Data (N) vs Time (s)");

% Plotting for the VGRF right foot vs time graph
hold on
plot(time_interval, right_interval);
legend("Left Foot", "Right Foot");
hold off

%% Part 2
% Plotting the VGRF signal in the frequency domain
freq = 1./time;

% % Plotting for the VGRF left foot vs frequency graph
% figure;
% plot(freq, left_foot);
% xlabel("Frequency (Hz)");
% ylabel("VGRF Data (N)");
% title("VGRF Data (N) vs Frequency (Hz)");
% 
% % Plotting for the VGRF right foot vs frequency graph
% hold on
% plot(freq, right_foot);
% legend("Left Foot", "Right Foot");
% xlim([1 5]);
% hold off

% Plotting the left VGRF signal in frequncy domain (Proper)
% Initialising the variable
% Code obtained from the lab
N = length(left_foot);
L = (1/N) * fft(left_foot);
X = -floor(N/2):(N-1-floor(N/2));

% % Creating a new figure
% figure; 
% stem(X, fftshift(abs(L)));
% title("Frequency Domain of Left Foot VGRF signal");
% ylabel("Magnitude");
% xlabel("k (multiple of fundemental frequency)");

% Left foot signal in Hertz
L = fft(left_foot);
omega = (-floor(N/2):(N-1-floor(N/2))) * (Fs/N);
figure; 
plot(omega, fftshift(abs(L)));
title("Frequency Domain of Left Foot VGRF signal");
ylabel("Magnitude");
xlabel("Frequency (Hz)");

% % Left foot signal in radian per sample
% L = fft(left_foot);
% omega = (-floor(N/2):(N-1-floor(N/2))) * (2*pi/N);
% figure; 
% plot(omega, fftshift(abs(L)));
% title('Frequency Domain of Left Foot VGRF signal');
% ylabel('Magnitude');
% xlabel('Frequency (rad/sample)');

% Plotting the right VGRF signal in frequncy domain
% Initialising the variable
N = length(right_foot);
R = (1/N) * fft(right_foot);
X = -floor(N/2):(N-1-floor(N/2));

% figure; 
% stem(X, fftshift(abs(R)));
% title("Frequency Domain of Right Foot VGRF signal");
% ylabel("Magnitude");
% xlabel("k (multiple of fundemental frequency)");

% Right foot signal in Hertz
R = fft(right_foot);
omega = (-floor(N/2):(N-1-floor(N/2))) * (Fs/N);
figure; 
plot(omega, fftshift(abs(R)));
title("Frequency Domain of Right Foot VGRF signal");
ylabel("Magnitude");
xlabel("Frequency (Hz)");

% %Right foot signal in radians per sample
% R = fft(right_foot);
% omega = (-floor(N/2):(N-1-floor(N/2))) * (2*pi/N);
% figure; 
% plot(omega, fftshift(abs(R)));
% title("Frequency Domain of Right Foot VGRF signal");
% ylabel("Magnitude");
% xlabel("Frequency (rad/sample)");

%% Part 3
% No, we cannot find a simple threshold in the VGRF value. 
% This is because the signal given consists of noise that blurs the transition between the swing and the stance phase.
% This could be due to some part of the noise having the same amplitude when the leg is lifted off the ground. 
% This means that if we filter out the signal (noise), we could potentially filter off the signal that is essential 
% to the analysis of our signal, leading to inaccurate calculation of the gait cycle which would lead us to make 
% inaccurate predictions.



