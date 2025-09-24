function saida = normaSistemaContinuo(ACal,ECal,CCal,DCal)
    sys = ss(ACal,ECal,CCal,DCal);
    [hinfJ,freqNorm] = hinfnorm(sys);
    saida.hInfNorm = hinfJ;
    saida.freqNorm = freqNorm;
end