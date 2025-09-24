function saida = verificaEstSisHib(ACal, KCal, h)
    mat = expm(h*ACal)*KCal;
    saida = eigs(mat, 1, 'lm');
end