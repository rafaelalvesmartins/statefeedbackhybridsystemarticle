 %% Unicamp - FEEC - 15/02/2020
close all;clear;clc;

%% ------------------Includes----------------------------------------------
addpath(genpath('funcoes'));

%% ------------------Matrizes do Sistema-----------------------------------
functionName = 'genSaveDataEx12';
eval([functionName '();']);

%% ------------------Chama a Função com todos as chamadas------------------
opt = []; 
chamaTudo(functionName,opt);

%% ------------------Remove os Includes------------------------------------                                                                                   ----------------------------------
rmpath(genpath('functions'));