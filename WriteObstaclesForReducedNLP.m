function WriteObstaclesForReducedNLP(xc, yc, s_ub)
global params_
Nfe = params_.Nfe;
Nv = params_.Nv;
Nobs = params_.Nobs;
oc = params_.Obs;
% This work assumes that the radius of each obstacle is equal to the radius of each vehicle disk
radius_of_obs = params_.R;

Mv2v = zeros(Nv, Nv, 4, 4, Nfe);
Mv2o = zeros(Nv, Nobs, 4, Nfe);
Mv2s = zeros(Nv, 4, 4, Nfe);

for ind = 1 : Nfe
    for ii = 1 : (Nv - 1)
        for jj = (ii + 1) : Nv
            for i = 1 : 4
                for j = 1 : 4
                    if ((xc(ii, i, ind) - xc(jj, j, ind))^2 + (yc(ii, i, ind) - yc(jj, j, ind))^2 < (params_.R + radius_of_obs)^2 + 0.01 + s_ub)
                        Mv2v(ii, jj, i, j, ind) = 1;
                    end
                end
            end
        end
    end
end

for ind = 1 : Nfe
    for ii = 1 : Nv
        for jj = 1 : Nobs
            for k = 1 : 4
                if ((xc(ii, k, ind) - oc{jj}.x)^2 + (yc(ii, k, ind) - oc{jj}.y)^2 < (params_.R + radius_of_obs)^2 + 0.01 + s_ub)
                    Mv2o(ii, jj, k, ind) = 1;
                end
            end
        end
    end
end

for ind = 1 : Nfe
    for ii = 1 : Nv
        for a = 1 : 3
            for b = (a+1) : 4
                if ((xc(ii, a, ind) - xc(ii, b, ind))^2 + (yc(ii, a, ind) - yc(ii, b, ind))^2 < (params_.R + radius_of_obs)^2)
                    Mv2s(ii, a, b, ind) = 1;
                end
            end
        end
    end
end

delete('M_V2V');
fid = fopen('M_V2V', 'w');
for ii = 1 : Nv
    for jj = 1 : Nv
        for a = 1 : 4
            for b = 1 : 4
                for kk = 1 : Nfe
                    fprintf(fid, '%g  %g  %g  %g  %g  %g\r\n', ii, jj, a, b, kk, Mv2v(ii,jj,a,b,kk));
                end
            end
        end
    end
end
fclose(fid);

delete('M_V2O');
fid = fopen('M_V2O', 'w');
for ii = 1 : Nv
    for jj = 1 : Nobs
        for a = 1 : 4
            for kk = 1 : Nfe
                fprintf(fid, '%g  %g  %g  %g  %g\r\n', ii, jj, a, kk, Mv2o(ii,jj,a,kk));
            end
        end
    end
end
fclose(fid);

delete('M_V2S');
fid = fopen('M_V2S', 'w');
for ii = 1 : Nv
    for a = 1 : 4
        for b = 1 : 4
            for kk = 1 : Nfe
                fprintf(fid, '%g  %g  %g  %g  %g\r\n', ii, a, b, kk, Mv2s(ii,a,b,kk));
            end
        end
    end
end
fclose(fid);
end