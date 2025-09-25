function saida = normaSistemaContinuo(A,B,C,D)
    sys = ss(A,B,C,D);
    [hinfJ,freqNorm] = hinfnorm(sys);
    saida.hInfNorm = hinfJ;
    saida.freqNorm = freqNorm;
end