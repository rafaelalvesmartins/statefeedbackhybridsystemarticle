%% Unicamp - FEEC - 05/07/2018
close all;clear;clc;

%% ------------------Includes----------------------------------------------
addpath(genpath('../functions'));

%% ------------------System's Matrices ------------------------------------
functionName = 'genSaveDataEx01';
eval(['paraSys=' functionName '();']);

%% ------------------Call Master Function with File Name ------------------
opt = ['Mid' 'Int' 'Poly']; % Options = Mid, Poly and Int
outPut = masterFunction(paraSys,functionName,opt);
save([functionName '/dataOutPut']);

%% ------------------Remove Includes---------------------------------------
rmpath(genpath('functions'));