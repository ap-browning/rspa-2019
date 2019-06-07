function cest = EstimateWavespeed(T,F)
%EstimateWavespeed Estimate the wavespeed from output of the continuous
%model [T,F]. The procedure is outlined in the supporting material
%document, section 3, page 4.

% Example usage:
%   [T,F]       = SpringsContinuous(pparams,nparams);
%   cest        = EstimateWavespeed(T,F);
    
% Region to lookat
tmax    = T(end);
lookat  = T > 0.8 * tmax & T < 0.9 * tmax;

Tfit    = T(lookat);
Ffit    = F(lookat);

pfit    = polyfit(Tfit,Ffit,1);
cest    = pfit(1);

end