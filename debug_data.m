% DEBUG_DATA - Check the structure of loaded data
fprintf('=== DEBUG: ESTRUTURA DOS DADOS ===\n');

try
    load('resultados_completos_massa_mola.mat', 'resultados');
    fprintf('Dados carregados com sucesso!\n');

    if exist('resultados', 'var')
        fprintf('Tipo: %s\n', class(resultados));

        if isstruct(resultados)
            campos = fieldnames(resultados);
            fprintf('Campos disponíveis (%d):\n', length(campos));
            for i = 1:length(campos)
                fprintf('  %d. %s\n', i, campos{i});

                % Mostrar estrutura do primeiro campo
                if i == 1
                    fprintf('     Estrutura de "%s":\n', campos{i});
                    if isstruct(resultados.(campos{i}))
                        subcampos = fieldnames(resultados.(campos{i}));
                        for j = 1:min(length(subcampos), 10)  % Max 10 subcampos
                            fprintf('       - %s\n', subcampos{j});
                        end
                        if length(subcampos) > 10
                            fprintf('       ... e mais %d campos\n', length(subcampos) - 10);
                        end
                    else
                        fprintf('       Tipo: %s\n', class(resultados.(campos{i})));
                    end
                end
            end
        else
            fprintf('Resultado não é struct, é: %s\n', class(resultados));
        end
    else
        fprintf('Variável "resultados" não encontrada no arquivo\n');
    end

catch ME
    fprintf('ERRO ao carregar: %s\n', ME.message);
end

fprintf('\n=== FIM DEBUG ===\n');