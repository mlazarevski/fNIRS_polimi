%% Read DPF TRS data%%
% Caterina Amendola, Alessandro Lia, Eleonora Mastromatteo

%   CURVES = READDAT(SETTINGS) reads the data according to the specified
%   settings and stores the result in CURVES.
%   
%
%   Input:
%   SETTINGS is a structure that must have the following fields:
%   - inputFile: the path to the .dat file to read
%   - nChannels: the number of channels of the acquisition board
%   - nCurvesTot: the number of curves in the file


addpath('C:\Users\catea\OneDrive - Politecnico di Milano\Codici Analisi\')
Settings.inputFile=('C:\Users\catea\OneDrive\Desktop\CW-NIRS\MISURE\TRS\MOTORIO\DATA\DPFm0001.dat');
Settings.nChannels=4096;
Settings.nCurvesTot=20;

curves=readDat(Settings);