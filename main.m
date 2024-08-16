% SS_main: this calls the main functions of the tDCS surround suppresion project
% subjID: e.g. S001

clear all;
clc;

SubjID   = 'jk';
sess_num = 1;
run_num  = 4; 

%%  
location_used = 4; 
offset = [0,0];
CurrDir = pwd;
SetupRand;
set_test_gamma;
HideCursor; 
parameters;

%%
WM_SS;

delete *.asv
