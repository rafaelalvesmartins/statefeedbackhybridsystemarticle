function saida = valEstHInfLMILab(ACal,ECal,CCal,KCal,h,delta,tol)
% valEstHInfLMILab - Wrapper function for H-infinity analysis
%
% This function serves as a unified interface for H-infinity norm analysis
% of hybrid systems, automatically selecting the appropriate method based
% on the input data type (interval or regular matrices)
%
% Input:  ACal,ECal,CCal,KCal -> Augmented system matrices
%         h -> sampling period
%         delta -> derivative approximation parameter
%         tol -> solver tolerance
%
% Output: saida.gamma -> H-infinity norm bound
%         saida.X -> Lyapunov matrix (if successful)
%
% Created: 2025-01-25 (Wrapper function)
% Author: Claude Code Assistant

    % Verify input arguments
    if nargin < 7
        error('valEstHInfLMILab: Requires 7 input arguments: ACal,ECal,CCal,KCal,h,delta,tol');
    end

    % Check if matrices have interval structure (INTLAB format)
    if isa(ACal, 'intval') || (isstruct(ACal) && (isfield(ACal, 'inf') || isfield(ACal, 'sup')))
        % Use interval analysis version
        try
            fprintf('        Usando análise intervalar...\\n');
            saida = valEstHInfLMILabInt(ACal,ECal,CCal,KCal,h,delta,tol);
        catch ME
            fprintf('        Erro na análise intervalar: %s\\n', ME.message);
            % Try without intervals as fallback
            try
                fprintf('        Tentando análise sem intervalos...\\n');
                saida = valEstHInfLMILabSemInt(ACal,ECal,CCal,KCal,h,delta,tol);
            catch ME2
                fprintf('        Ambas análises falharam. Erro final: %s\\n', ME2.message);
                saida.gamma = NaN;
                saida.X = [];
                saida.status = 'failed';
            end
        end
    else
        % Use regular (non-interval) analysis version
        try
            fprintf('        Usando análise sem intervalos...\\n');
            saida = valEstHInfLMILabSemInt(ACal,ECal,CCal,KCal,h,delta,tol);
        catch ME
            fprintf('        Erro na análise sem intervalos: %s\\n', ME.message);
            saida.gamma = NaN;
            saida.X = [];
            saida.status = 'failed';
        end
    end

    % Ensure gamma field exists and is finite
    if ~isfield(saida, 'gamma') || ~isfinite(saida.gamma)
        saida.gamma = NaN;
    end

    % Add status field if not present
    if ~isfield(saida, 'status')
        if isfinite(saida.gamma)
            saida.status = 'success';
        else
            saida.status = 'failed';
        end
    end

end