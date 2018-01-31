more off;
!synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0
!xinput set-prop 15 "Synaptics Two-Finger Scrolling" 1 0

%MATLAB_HOME='/home/user';
MATLAB_HOME='/home/johann/MATLAB';

path(path,[MATLAB_HOME filesep 'Fitting/ezyfit']);
path(path,[MATLAB_HOME filesep 'Graph/Analysis/Graph1/']);

path(path,'plot');
path(path,'graph');
path(path,'sim');
path(path,'data');
path(path,'analysis');

rand('seed',99);
randn('seed',17);
