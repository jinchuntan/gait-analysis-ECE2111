% Assignment 1
% Written by Tan Jin Chun
% Last Modified: 24/9/2021
% File Name: estimateGait.m

%% Section 3.3
function [STl, STr, SWl, SWr, Sl, Sr] = estimateGait(VGRF)

% Name: estimateGait

% Purpose: There are three main purpose of this function. The first purpose
% is to filter the noise frequency itself. The second purpose is to find
% the start time for stance and swing in each stride gait cycle and the
% last purpose is to find the duration of the swing/stance/stride phases
% from the swing and stance onset times.

% Input
% VGRF: the raw VGRF signals

% Outputs
% STl vector of the stance times (durations of stance phases) for the left foot
% STr vector of the stance times for the right foot
% SWl vector of the swing times (durations of swing phases) for the left foot
% SWr vector of the swing times for the right foot
% Sl vector of the stride times (durations of stride phases) for the left foot
% Sr vector of the stride times for the right foot

% Logic of the code
% The logic behind my code is that we begin with getting the data from the user. 
% The data will then be extracted into its respective categories. 
% A low pass filter will be created within the function itself and the data for the left foot and the right foot will undergo convolution for filtering. 
% We then get the starting index of the left foot. 
% The swing and the stance phase value will be obtained periodically, meaning that once the swing phase is obtained, we will start to find for the stance phase. 
% This is done by going through the filtered data and a for loop with an if-elseif condition. 
% After this, the swing and the stance onset time is obtained. 
% From here, we would need to get the index of all of the values obtained using the find MATLAB built-in function.
% Once the index is obtained, we then can find the stance duration, swing duration and the strides duration. 
% The stance duration would be the swing index minus the stance index divided by the sampling frequency of 120Hz. 
% The swing duration would be the stance index minus the swing index divided by the sampling frequency. 
% The stride duration can be directly obtained from the stance index since the stride durations is just the summation of the stance and the swing durations.

    % Initialising the variable
    Fs = 120;
    
    % Extracting the data
    time = VGRF(:,1);
    left_foot = VGRF(:,2);
    right_foot = VGRF(:,3);

    % Making a low pass filter
    cutoff_Hz = [10, 15];
    w = 2 * cutoff_Hz / Fs;
    L = 200;
    h = firpm(L, [0, w, 1], [1 1 0 0]); 
    
    % Convolution
    left_foot = conv(left_foot, h);
    right_foot = conv(right_foot, h);
    
    % Threshold value
    threshold = 30;
    
    % Preallocating the variables
    L_Swing_Phase = 0;
    L_Stance_Phase = 0;
    L_Stance_Onset = zeros(1, length(left_foot));
    L_Swing_Onset = zeros(1, length(left_foot));
    
    % The place to start
    L_Start = floor((L-1)/2 + length(L) - 1);
    
    % Using a for loop to loop through the data
    for i = L_Start:length(left_foot)
    % If the value of the left foot is bigger or equal than the threshold
    % and it is not the left stance phase
        if ((~L_Stance_Phase) && left_foot(i) >= threshold)
            
            % Set stance phase to one and swing phase to 0
            L_Stance_Phase = 1;
            L_Swing_Phase = 0;
            L_Stance_Onset(i) = 1;
    
    % If the threshold value is larger than the signal and it is not the
    % left swing phase
        elseif ((~L_Swing_Phase) && left_foot(i) < threshold)
            
            % Set stance phase to 0 and swing phase to one
            L_Stance_Phase = 0;
            L_Swing_Phase = 1;
            L_Swing_Onset(i) = 1;
            
        end
    end
    
    % Getting the starting indexes and starting value of the left foots for the
    % stance and the swing times
    % Left Foot Index
    L_Stance_Index = find(L_Stance_Onset);
    L_Swing_Index = find(L_Swing_Onset);
    
    % Duration for left foot
    L_Stance_Durations = (L_Swing_Index - L_Stance_Index) / Fs;
    L_Swing_Durations = (L_Stance_Index(2:end) - L_Swing_Index(1:end-1)) / Fs;
    L_Stride_Durations = (L_Stance_Index(2:end) - L_Stance_Index(1:end-1)) / Fs;
    
    % Updating the duration
    L_Stance_Durations = L_Stance_Durations(2:end-2);
    L_Swing_Durations = L_Swing_Durations(1:end-1);
    L_Stride_Durations = L_Stride_Durations(2:end-1);
    
    % Repeat the same process as above for the right foot as well
    % Right Foot
    R_Swing_Phase = 0;
    R_Stance_Phase = 0;
    R_Stance_Onset = zeros(1, length(right_foot));
    R_Swing_Onset = zeros(1, length(right_foot));
    
    % The place to start (For right foot)
    R_Start = floor((L-1)/2 + length(L) - 1);
    
    % Using a for loop to loop through the data
    for i = R_Start:length(right_foot)
        
    % If the value of the right foot is bigger or equal than the threshold
    % and it is not the right stance phase
        if ((~R_Stance_Phase) && right_foot(i) >= threshold)
            
            % Set stance phase to one and swing phase to 0
            R_Stance_Phase = 1;
            R_Swing_Phase = 0;
            R_Stance_Onset(i) = 1;
            
    % If the threshold value is larger than the data and it is not the
    % right swing phase
        elseif ((~R_Swing_Phase) && right_foot(i) < threshold)
            
            % Set stance phase to 0 and swing phase to one
            R_Stance_Phase = 0;
            R_Swing_Phase = 1;
            R_Swing_Onset(i) = 1;
            
        end

    end
    
    % Right Foot Index
    R_Stance_Index = find(R_Stance_Onset);
    R_Swing_Index = find(R_Swing_Onset);
    
    % Duration for right_
    R_Stance_Duration = (R_Swing_Index - R_Stance_Index) / Fs;
    R_Swing_Duration = (R_Stance_Index(2:end) - R_Swing_Index(1:end-1))/Fs;
    R_Stride_Duration = (R_Stance_Index(2:end) - R_Stance_Index(1:end-1))/Fs;
    
    % Updating the durations
    R_Stance_Duration = R_Stance_Duration(2:end-1);
    R_Swing_Duration = R_Swing_Duration(1:end-1);
    R_Stride_Duration = R_Stride_Duration(1:end-1);
    
    % Gait Parameter Calculations
    % Initialising the variable here
    Interval = 5;
    Start_Time = 10*Fs+1;
    End_Time = Start_Time + Interval * Fs;
    Indices = Start_Time:End_Time;
    Time_Interval = time(Indices);
    Left_Interval = left_foot(Indices);
    Right_Interval = right_foot(Indices);
    
    % Initialising the variable for the plots
    L_Stance_Val = left_foot(L_Stance_Index);
    L_Swing_Val = left_foot(L_Swing_Index);
    R_Stance_Val = right_foot(R_Stance_Index);
    R_Swing_Val = right_foot(R_Swing_Index);
    
    % For the left foot
    L_Stance_Plot = L_Stance_Index(L_Stance_Index >= min(Indices) & L_Stance_Index <= max(Indices)) / Fs - 1 / Fs;
    L_Swing_Plot = L_Swing_Index(L_Swing_Index >= min(Indices) & L_Swing_Index <= max(Indices)) / Fs - 1 / Fs;
    L_Stance_Plot_Val = L_Stance_Val(L_Stance_Index >= min(Indices) & L_Stance_Index <= max(Indices));
    L_Swing_Plot_Val = L_Swing_Val(L_Swing_Index >= min(Indices) & L_Swing_Index <= max(Indices));

    % For the right foot
    R_Stance_Plot = R_Stance_Index(R_Stance_Index >= min(Indices) & R_Stance_Index <= max(Indices)) / Fs - 1 / Fs;
    R_Swing_Plot = R_Swing_Index(R_Swing_Index >= min(Indices) & R_Swing_Index <= max(Indices)) / Fs - 1 / Fs;
    R_Stance_Plot_Val = R_Stance_Val(R_Stance_Index >= min(Indices) & R_Stance_Index <= max(Indices));
    R_Swing_Plot_Val = R_Swing_Val(R_Swing_Index >= min(Indices) & R_Swing_Index <= max(Indices));

    % Plotting the filtered VGRF
    figure;
    plot(Time_Interval, Left_Interval, Time_Interval, Right_Interval);
    hold on;
    plot(L_Stance_Plot, L_Stance_Plot_Val, 'bd', L_Swing_Plot, L_Swing_Plot_Val ,'bo');
    plot(R_Stance_Plot, R_Stance_Plot_Val, 'rd', R_Swing_Plot, R_Swing_Plot_Val, 'ro');
    
    % Labelling the plot
    xlabel("Time (s)");
    ylabel("VGRF Data (N)"); 
    title("Filtered VGRF (N) vs Time (s)");
    legend("left foot", "right foot", "left stance onset", "left swing onset", "right stance onset", "right swing onset");
    hold off;
    
    % Plotting the frequency domain of the left foot after filtering
    % Code obtained from the lab 
    % Initialising the variable
    Left_N = length(left_foot);
    Left_X = fft(left_foot);
    Left_Omega = (-floor(Left_N/2):(Left_N-1-floor(Left_N/2))) * (Fs/Left_N);
    
    % Plotting the figure
    figure; 
    plot(Left_Omega, fftshift(abs(Left_X)));
    
    % Labelling the figure
    xlabel("Frequency (Hz)");
    ylabel("Magnitude");
    title("Frequency Domain of Left Foot Filtered VGRF Signal");

    % Plotting the frequency domain of the right foot after filtering
    % Initialising the variable
    Right_N = length(right_foot);
    Right_X = fft(right_foot);
    Right_Omega = (-floor(Right_N/2):(Right_N-1-floor(Right_N/2))) * (Fs/Right_N);
    
    % Plotting the figure
    figure; 
    plot(Right_Omega, fftshift(abs(Right_X)));
    
    % Labelling the plot 
    xlabel("Frequency (Hz)");
    ylabel("Magnitude");
    title("Frequency Domain of Right Foot Filtered VGRF Signal");
    
    % Plotting 3 subplots
    % Plotting the subplot for the left and the right stance
    % Left Stance
    figure;
    subplot(2,1,1);
    stem(L_Stance_Durations);
    xlabel("Gait Cycle");
    ylabel("Stance Time (s)");
    title("Left Foot Stance Time Against Gait Cycle");
    
    % Right Stance
    subplot(2,1,2);
    stem(R_Stance_Duration);
    xlabel("Gait Cycle");
    ylabel("Stance Time (s)");
    title("Right Foot Stance Time Against Gait Cycle")
    
    % Plotting the subplot for the swing time
    % Left Swing
    figure;
    subplot(2,1,1);
    stem(L_Swing_Durations);
    xlabel("Gait Cycle");
    ylabel("Swing Time (s)");
    title("Left Foot Swing Time Against Gait Cycle");
    
    % Right Swing
    subplot(2,1,2);
    stem(R_Swing_Duration);
    xlabel("Gait Cycle");
    ylabel("Swing Time (s)");
    title("Right Foot Swing Time Against Gait Cycle")
    
    % Plotting the subplot for left and right stride
    % Left Stride
    figure;
    subplot(2,1,1);
    stem(L_Stride_Durations);
    xlabel("Gait Cycle");
    ylabel("Stride Time (s)");
    title("Left Foot Stride Time Against Gait Cycle");
    
    % Right Stride
    subplot(2,1,2);
    stem(R_Stride_Duration);
    xlabel("Gait Cycle");
    ylabel("Stride Time (s)");
    title("Right Foot Stride Time Against Gait Cycle")
    
    
    % Setting the values for the outputs
    STl = L_Stance_Durations;
    STr = R_Stance_Duration;
    SWl = L_Swing_Durations;
    SWr = R_Swing_Duration;
    Sl = L_Stride_Durations;
    Sr = R_Stride_Duration;
end