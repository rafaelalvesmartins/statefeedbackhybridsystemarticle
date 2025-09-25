function teste_simulacao_rapida()
% TESTE_SIMULACAO_RAPIDA - Quick test of the fixed simulation functions
%
% This function tests the corrected simulation with a small subset of data
% to verify everything works before running the full analysis

    fprintf('=== TESTE RÁPIDO DA SIMULAÇÃO CORRIGIDA ===\n\n');

    % Setup paths
    setup_paths();

    % Try loading real data first
    try
        fprintf('Carregando dados reais...\n');
        load('resultados_completos_massa_mola.mat', 'resultados');
        fprintf('✓ Dados carregados!\n');
        use_real_data = true;
    catch
        fprintf('⚠ Dados não encontrados, usando sintéticos...\n');
        use_real_data = false;
    end

    if use_real_data && exist('resultados', 'var')
        % Use only first controller for quick test
        controladores = fieldnames(resultados);
        controlador_teste = controladores{1};  % First controller

        fprintf('\nTestando controlador: %s\n', controlador_teste);

        try
            % Get controller data
            dados = resultados.(controlador_teste);

            % Test parameters - use smaller values for speed
            h = 0.1;
            delta = 0.01;
            tol = 1e-4;    % Reduced tolerance for speed
            imageName = 'teste_rapido';
            axisVector = [0 5 -0.5 0.5];  % Shorter time

            % Take only first polytope vertex for quick test
            combPoly = dados.combPoly;

            % Reduce to single vertex
            single_vertex_poly = struct();
            single_vertex_poly.A.alphaVecs = combPoly.A.alphaVecs(1);
            single_vertex_poly.A.polytopicMatrices = combPoly.A.polytopicMatrices(1);
            single_vertex_poly.B.alphaVecs = combPoly.B.alphaVecs(1);
            single_vertex_poly.B.polytopicMatrices = combPoly.B.polytopicMatrices(1);
            single_vertex_poly.E.alphaVecs = combPoly.E.alphaVecs(1);
            single_vertex_poly.E.polytopicMatrices = combPoly.E.polytopicMatrices(1);
            single_vertex_poly.C.alphaVecs = combPoly.C.alphaVecs(1);
            single_vertex_poly.C.polytopicMatrices = combPoly.C.polytopicMatrices(1);
            single_vertex_poly.D.alphaVecs = combPoly.D.alphaVecs(1);
            single_vertex_poly.D.polytopicMatrices = combPoly.D.polytopicMatrices(1);
            single_vertex_poly.D1.alphaVecs = combPoly.D1.alphaVecs(1);
            single_vertex_poly.D1.polytopicMatrices = combPoly.D1.polytopicMatrices(1);

            fprintf('  Sistema: %dx%d, 1 vértice (teste rápido)\n', ...
                size(combPoly.A.polytopicMatrices{1}, 1), ...
                size(combPoly.B.polytopicMatrices{1}, 2));

            % Test both analysis types
            fprintf('\n--- TESTANDO ANÁLISE INTERVALAR ---\n');
            tic;
            saida_int = simulatesSampledInputRobusto(single_vertex_poly, h, dados.K, imageName, axisVector, delta, tol, false);
            tempo_int = toc;

            fprintf('\n--- TESTANDO ANÁLISE POLITÓPICA ---\n');
            tic;
            saida_poly = simulatesSampledInputRobusto(single_vertex_poly, h, dados.K, imageName, axisVector, delta, tol, true);
            tempo_poly = toc;

            % Report results
            fprintf('\n=== RESULTADOS DO TESTE RÁPIDO ===\n');
            fprintf('Controlador testado: %s\n', controlador_teste);

            fprintf('\nAnálise Intervalar (%.2fs):\n', tempo_int);
            if isfield(saida_int, 'gamma') && isfield(saida_int.gamma, 'maxGamma')
                if isfinite(saida_int.gamma.maxGamma)
                    fprintf('  γ máximo: %.6f\n', saida_int.gamma.maxGamma);
                    fprintf('  Taxa sucesso LMI: %.1f%%\n', saida_int.lmi_success_rate*100);
                else
                    fprintf('  ⚠ LMI falhou\n');
                end
            end

            fprintf('\nAnálise Politópica (%.2fs):\n', tempo_poly);
            if isfield(saida_poly, 'gamma') && isfield(saida_poly.gamma, 'maxGamma')
                if isfinite(saida_poly.gamma.maxGamma)
                    fprintf('  γ máximo: %.6f\n', saida_poly.gamma.maxGamma);
                    fprintf('  Taxa sucesso LMI: %.1f%%\n', saida_poly.lmi_success_rate*100);
                else
                    fprintf('  ⚠ LMI falhou\n');
                end
            end

            % Check if at least one worked
            int_ok = isfield(saida_int, 'gamma') && isfinite(saida_int.gamma.maxGamma);
            poly_ok = isfield(saida_poly, 'gamma') && isfinite(saida_poly.gamma.maxGamma);

            if int_ok || poly_ok
                fprintf('\n✅ TESTE APROVADO! Pelo menos uma análise funcionou.\n');
                if int_ok && poly_ok
                    fprintf('   Ambas análises funcionaram - sistema pronto!\n');
                end
            else
                fprintf('\n⚠ TESTE PARCIAL - LMIs falharam mas funções executaram sem erro.\n');
            end

        catch ME
            fprintf('\n✗ ERRO no teste: %s\n', ME.message);
            if ~isempty(ME.stack)
                fprintf('   Em: %s (linha %d)\n', ME.stack(1).name, ME.stack(1).line);
            end
        end

    else
        fprintf('\n⚠ Dados reais não disponíveis. Execute primeiro:\n');
        fprintf('   main_massa_mola_completo()\n');
    end

    fprintf('\n=== TESTE RÁPIDO CONCLUÍDO ===\n');

end