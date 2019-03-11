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
%
%%
function [odf,odf_extract,ori]=PlotODFandExtractOrientations(path,fnamesPrefix,levelsODF,Nori,Niter,GaussianSmoothing,SeriesRank)
%% Import pole figure data and create PoleFigure object
cs = crystalSymmetry('m-3m',[4.04 4.04 4.04],'mineral','Al');

fnames = {
    [path fnamesPrefix '_pf111_uncorr.dat'],...
    [path fnamesPrefix '_pf200_uncorr.dat'],...
    [path fnamesPrefix '_pf220_uncorr.dat'],...
    [path fnamesPrefix '_pf311_uncorr.dat']};

% Specimen symmetry
ss = specimenSymmetry('1'); % Triclinic
ssO = specimenSymmetry('orthorhombic');

% Plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

% Set annotations to highlight spatial reference frame
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],...
    {'RD','TD'},'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);

h = {
    Miller(1,1,1,cs),...
    Miller(2,0,0,cs),...
    Miller(2,2,0,cs),...
    Miller(3,1,1,cs)};

% Load pole figures separately
columnNames = {'Polar Angle','Azimuth Angle','Intensity'};
pf1 = loadPoleFigure_generic(fnames{1},'ColumnNames',columnNames);
pf2 = loadPoleFigure_generic(fnames{2},'ColumnNames',columnNames);
pf3 = loadPoleFigure_generic(fnames{3},'ColumnNames',columnNames);
pf4 = loadPoleFigure_generic(fnames{4},'ColumnNames',columnNames);

% Construct pole figure object of the four pole figures
intensities = {
    pf1.intensities,...
    pf2.intensities,...
    pf3.intensities,...
    pf4.intensities};
pfs = PoleFigure(h,pf1.r,intensities,cs,ss);

%% Plot pole figures of raw, corrected data
% figure
% plot(pfs,'upper','projection','eangle','minmax')
% mtexColorbar('location','southOutside')

%% Calculate the ODF using default settings
disp('Calculating ODF')
odf = calcODF(pfs);

% Set correct specimen symmetry for calculation of texture strength
odf.SS = ssO;

% Calculate texture strength
% textureIndex = odf.textureindex
% entropy = odf.entropy
% odfMax = odf.max

%% Extract N orientations from ODF
disp(['Extracting ',num2str(Nori),' orientations from ODF'])

error=inf;
progress(0,Niter)
for k=1:Niter
    ori_tmp=calcOrientations(odf,Nori);

    % Calculating ODF from extracted orientations
%     psi_tmp = calcKernel(ori_tmp);
%     odf_extract_tmp = calcODF(ori_tmp,'kernel',psi_tmp);
    odf_extract_tmp = calcODF(ori_tmp);
    odf_extract_tmp.SS=ssO;
    error_tmp = calcError(odf_extract_tmp,odf,'resolution',5*degree);
    if error_tmp<error
        ori=ori_tmp;
%         psi=psi_tmp;
        odf_extract=odf_extract_tmp;
        error=error_tmp;
    end
    progress(k,Niter)
end

%% Write orientations to .ori file for Auswert
[phi1,Phi,phi2] = Euler(ori,'Bunge');
% converting from radians to degrees
phi1 = phi1/degree;
Phi = Phi/degree;
phi2 = phi2/degree;

ID=fopen([path,fnamesPrefix,'_mtex.ori'],'w');

fprintf(ID,'%s \n','TEX      EPS 0.000   texture from mtex');
fprintf(ID,'%s %2d \n','PHI2',SeriesRank);
fprintf(ID,'%d %3d %2.1f \n',Nori,0,GaussianSmoothing);
for i=1:Nori
    fprintf(ID,'%6.2f %6.2f %6.2f %8.6f %5.2f \n',phi1(i),Phi(i),phi2(i),1,0);
end
fcl=fclose(ID);

%% Plot ODF in Euler space phi2 sections
disp('Plotting ODFs')

% figure
% plot(odf,'phi2',[0 45 65]*degree,'contourf','minmax')
% mtexColorMap white2black

figure
plot(odf,'phi2',[0 45 65]*degree,'contourf',levelsODF,'minmax')
mtexColorMap white2black
mtexColorbar('location','southoutside')

export_fig([path fnamesPrefix '_odf.pdf'])

% figure
% plot(odf_extract,'phi2',[0 45 65]*degree,'contourf','minmax')
% mtexColorMap white2black

figure
plot(odf_extract,'phi2',[0 45 65]*degree,'contourf',levelsODF,'minmax')
mtexColorMap white2black
mtexColorbar('location','southoutside')

export_fig([path fnamesPrefix '_estimated_odf.pdf'])

% figure
% plot(odf,'sections',18,'contourf',levelsODF)
% mtexColorMap white2black

% figure
% plot(odf_extract,'sections',18,'contourf',levelsODF)
% mtexColorMap white2black
end