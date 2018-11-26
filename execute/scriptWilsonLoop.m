clear; close all; clc;
Ns = [149,118,70];
Ncirc = [120,70,36];

for ii = 1:numel(Ns)
    for adivd = [1.8,2,2.1,2.2]
        close all;
        fprintf('\n++++++++ ii = %g (Ns = %g, Ncirc = %g) ++++++++',ii,Ns(ii),Ncirc(ii))
        fprintf('\n++++++++ adivd = %g ++++++++\n',adivd)
        
        runWilsonLoop('triangular',Ns(ii),Ncirc(ii),adivd,[],[],[],1)
    end
end