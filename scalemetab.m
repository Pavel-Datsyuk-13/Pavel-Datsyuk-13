function outspec = scalemetab(instore, met)
% Function that scales simulation data to the metabolite per govind paper/
% our in vivo data - Only used for RMS simulations in this case
% Justin - 5-20-2025

met = lower(met);
switch met
    case 'gln'
        scale = 4.0;
    case 'glu'
        scale = 8.0;
    case 'thr'
        scale = 0.33;
    case 'glc'
        scale = 1.16;
    case 'gly'
        scale = 0.9;
    case 'ins' 
        scale = 5.58;
    otherwise
        error('Unknown metabolite: %s, please edit the funciton to include it', met);
end
outspec = scale*instore;
end