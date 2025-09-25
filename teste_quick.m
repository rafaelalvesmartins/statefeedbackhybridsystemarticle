% TESTE_QUICK - Quick test of fixed functions
fprintf('=== TESTE RÁPIDO DAS CORREÇÕES ===\n');

% Setup paths
setup_paths();

% Test if functions exist
fprintf('\nVerificando funções:\n');
functions_to_test = {'valEstHInfLMILab', 'estHInfAnaPolyLMILab'};

for i = 1:length(functions_to_test)
    if exist(functions_to_test{i}, 'file')
        fprintf('✓ %s encontrada\n', functions_to_test{i});
    else
        fprintf('✗ %s NÃO encontrada\n', functions_to_test{i});
    end
end

% Try a simple test simulation
fprintf('\n=== TESTE SINTÉTICO RÁPIDO ===\n');
try
    % Create simple 2x2 test system
    nx = 2; nu = 1; nw = 1;

    % Simple stable system
    A_test = [0 1; -2 -3];
    B_test = [0; 1];
    E_test = [1; 0];
    C_test = [1 0];
    D_test = 0;

    % Test parameters
    h = 0.1;
    K = [2 1];  % Simple stabilizing gain
    delta = 0.01;
    tol = 1e-6;

    % Create caligraphic matrices
    ACal = [A_test B_test; zeros(nu, nx) zeros(nu, nu)];
    ECal = [E_test; zeros(nu, nw)];
    CCal = [C_test D_test];
    KCal = [eye(nx) zeros(nx, nu); K zeros(nu)];

    fprintf('Testando valEstHInfLMILab...\n');
    saida = valEstHInfLMILab(ACal, ECal, CCal, KCal, h, delta, tol);

    if isfield(saida, 'gamma') && isfinite(saida.gamma)
        fprintf('✓ SUCCESS! gamma = %.4f\n', saida.gamma);
    else
        fprintf('⚠ LMI falhou, mas função foi chamada sem erro\n');
    end

catch ME
    fprintf('✗ ERRO no teste: %s\n', ME.message);
end

fprintf('\n=== TESTE CONCLUÍDO ===\n');