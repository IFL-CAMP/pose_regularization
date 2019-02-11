
addpath('..', '../wrapper');

replicate_denoising('../../test/data/case5_preset_optical.mat');
replicate_denoising('../../test/data/case5_preset_em1.mat');
replicate_denoising('../../test/data/case5_preset_em2.mat');

function replicate_denoising(filename)

vars = load(filename);

result = regularizeMatrices4x4(vars.Poses_Noisy, vars.p, vars.q, vars.r, ...
    vars.alpha, vars.beta, vars.inner_factor, vars.steps, vars.inner_steps);

isequal(result, vars.Poses_Regularized)


end


