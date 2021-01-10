function [xxc, yyc] = LoadXY()
global params_
load xc.txt
load yc.txt

Nfe = params_.Nfe;
Nv = params_.Nv;

xxc = zeros(Nv, 4, Nfe);
yyc = zeros(Nv, 4, Nfe);

counter = 0;
for ii = 1 : Nv
    for jj = 1 : 4
        for kk = 1 : Nfe
            counter = counter + 1;
            xxc(ii, jj, kk) = xc(counter);
            yyc(ii, jj, kk) = yc(counter);
        end
    end
end
end