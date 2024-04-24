% Assignment 1
% Written by Tan Jin Chun - 32194471
% Last Modified: 21/10/2021
% File Name: section33.m

clear all; clc; close all;

% Loading the data
load ("assignment_data.mat");

% Calling the function estimateGait
% Putting the raw data into the function
[STl,STr,SWl,SWr,Sl,Sr] = estimateGait(C1);

% Printing out the values
fprintf("The mean for stride (left foot) is %.4f.\n", mean(Sl));
fprintf("The mean for stride (right foot) is %.4f.\n", mean(Sr));
fprintf("The mean for stance (left foot) is %.4f.\n", mean(STl));
fprintf("The mean for stance (right foot) is %.4f.\n", mean(STr));
fprintf("The mean for stance (right foot) is %.4f.\n", mean(SWl));
fprintf("The mean for stance (right foot) is %.4f.\n", mean(SWr));