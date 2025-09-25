function setup_paths()
% SETUP_PATHS - Add all necessary directories to MATLAB path
%
% This function ensures that all LMI analysis functions are accessible
% by adding their directories to the MATLAB path
%
% Created: 2025-01-25
% Author: Claude Code Assistant

    % Get the base directory
    baseDir = pwd;

    fprintf('Configurando paths para LMI functions...\\n');

    % List of subdirectories to add to path
    paths_to_add = {
        fullfile(baseDir, 'funcoes'),
        fullfile(baseDir, 'funcoes', 'Diversas'),
        fullfile(baseDir, 'funcoes', 'Optimize'),
        fullfile(baseDir, 'HInf - Análise - Intervalar', 'funcoes'),
        fullfile(baseDir, 'HInf - Análise - Intervalar', 'funcoes', 'AnaliseSemInt'),
        fullfile(baseDir, 'HInf - Análise - Intervalar', 'funcoes', 'Diversas'),
        fullfile(baseDir, 'functions'),
        fullfile(baseDir, 'functions', 'Optimize')
    };

    % Add paths that exist
    paths_added = 0;
    for i = 1:length(paths_to_add)
        if exist(paths_to_add{i}, 'dir')
            addpath(paths_to_add{i});
            fprintf('  ✓ Adicionado: %s\\n', paths_to_add{i});
            paths_added = paths_added + 1;
        else
            fprintf('  ⚠ Não encontrado: %s\\n', paths_to_add{i});
        end
    end

    fprintf('Total de %d paths adicionados ao MATLAB.\\n', paths_added);

    % Verify critical functions are now available
    critical_functions = {
        'valEstHInfLMILab',
        'valEstHInfLMILabSemInt',
        'valEstHInfLMILabInt',
        'estHInfAnaPolyLMILab',
        'estHInfAnaPolyLMILabOriginal'
    };

    fprintf('\\nVerificando funções críticas:\\n');
    for i = 1:length(critical_functions)
        if exist(critical_functions{i}, 'file')
            fprintf('  ✓ %s\\n', critical_functions{i});
        else
            fprintf('  ✗ %s (não encontrada)\\n', critical_functions{i});
        end
    end

    fprintf('\\nSetup de paths concluído!\\n');

end