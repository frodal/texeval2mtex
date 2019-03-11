%% TexEval to MTEX to discrete orientations
%
% Creates a MTEX PoleFigure object from four pole figures
% Plots the ODF and extracts a set of discrete orientations that represent the ODF
% Writes the orientations to a file compatible with the Auswert program by Olaf Engler
%
% Author  : Bjørn Håkon Frodal
% Contact : bjorn.h.frodal@ntnu.no
%
% Working with:
%   MATLAB R2018a
%   MTEX-5.1.1
%
% Requires:
%   MTEX (Available here:
%   http://mtex-toolbox.github.io/download.html)
%   export_fig (Available here:
%   https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig)

clf
close all
clear
clc

%% User input

% path and filename prefix to the pole figure data
path = ['./example_data/'];
fnamesPrefix = 'polefigure_intensities';

% Levels of the ODF contour plot
levelsODF = [0,1,2,4,7,10,15,20,30,40,50,75,100];

% Number of discrete orientations to extract from ODF
Nori=1000;
% Number of iterations to find the best discrete orientations from ODF
Niter=100;

% Auswert .ori settings
GaussianSmoothing=7.0;
SeriesRank=23;

%% Plot ODF and extract orientations

[odf_experimental,odf_orientations,orientations]=PlotODFandExtractOrientations(path,fnamesPrefix,levelsODF,Nori,Niter,GaussianSmoothing,SeriesRank);
