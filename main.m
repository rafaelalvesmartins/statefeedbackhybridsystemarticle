%% Unicamp - FEEC - 28/03/2020
close all;clear;clc;

%% ------------------Includes----------------------------------------------
addpath(genpath('funcoes'));

%% ------------------Matrizes do Sistema-----------------------------------
functionName = 'genSaveDataEx12NomCont';
eval(['paraSys=' functionName '();']);

%% ------------------Chama a Função com todos as chamadas------------------
opt = []; 
outPut = chamaTudo(paraSys,functionName,opt);
save([functionName '/dataOutPut']);

%% ------------------Remove os Includes------------------------------------                                                                                 
rmpath(genpath('functions'));