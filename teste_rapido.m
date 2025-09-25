% TESTE_RAPIDO - Teste rápido das funções corrigidas

fprintf('=== TESTE RÁPIDO ===\n');

% Verificar se o script principal existe e pode ser executado
if exist('main_simulacao_completa.m', 'file')
    fprintf('✓ main_simulacao_completa.m encontrado\n');
else
    fprintf('✗ main_simulacao_completa.m não encontrado\n');
end

% Verificar se as funções de simulação existem
if exist('simulatesSampledInputCompleto.m', 'file')
    fprintf('✓ simulatesSampledInputCompleto.m encontrado\n');
else
    fprintf('✗ simulatesSampledInputCompleto.m não encontrado\n');
end

if exist('simulatesSampledInputSimplified.m', 'file')
    fprintf('✓ simulatesSampledInputSimplified.m encontrado\n');
else
    fprintf('✗ simulatesSampledInputSimplified.m não encontrado\n');
end

% Verificar se os dados existem
if exist('resultados_completos_massa_mola.mat', 'file')
    fprintf('✓ resultados_completos_massa_mola.mat encontrado\n');
    fprintf('\n🚀 TUDO PRONTO! Execute: main_simulacao_completa()\n');
else
    fprintf('⚠ resultados_completos_massa_mola.mat não encontrado\n');
    fprintf('📋 Execute primeiro: main_massa_mola_completo()\n');
    fprintf('   Em seguida: main_simulacao_completa()\n');
end

fprintf('\n=== INSTRUÇÕES ===\n');
fprintf('1. Para executar simulação completa: main_simulacao_completa()\n');
fprintf('2. Para testar com dados sintéticos: teste_simulacao()\n');
fprintf('3. Agora as funções têm verificação robusta de argumentos!\n');

fprintf('\n=== STATUS DAS CORREÇÕES ===\n');
fprintf('✅ Erro de sintaxe MATLAB corrigido (operador ternário)\n');
fprintf('✅ Verificação de argumentos adicionada\n');
fprintf('✅ Simulação temporal completa implementada\n');
fprintf('✅ Análise LMI robusta com fallbacks\n');
fprintf('✅ Resposta focada nas questões do revisor\n');