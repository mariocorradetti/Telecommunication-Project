function [kH, kV, alphaH, alphaV] = attenuationCoeffs(f)
%frequency dependent coefficient ITU-R P838 
coeff_kH.a_vec = [-5.33980 -0.35351 -0.23789 -0.94158]';
coeff_kH.b_vec = [-0.10008 1.26970  0.86036  0.64552]';
coeff_kH.c_vec = [1.13098  0.45400  0.15354  0.16817]';
coeff_kH.m     = -0.18961;
coeff_kH.c     = 0.71147;

coeff_kV.a_vec = [-3.80595 -3.44965 -0.39902 0.50167]';
coeff_kV.b_vec = [0.56934 -0.22911 0.73042 1.07319]';
coeff_kV.c_vec = [0.81061 0.51059 0.11899 0.27195]';
coeff_kV.m     = -0.16398;
coeff_kV.c     = 0.63297;


coeff_alphaH.a_vec = [-0.14318 0.29591 0.32177 -5.37610 16.1721]';
coeff_alphaH.b_vec = [1.82442 0.77564 0.63773 -0.96230 -3.29980]';
coeff_alphaH.c_vec = [-0.55187 0.19822 0.13164 1.47828 3.43990]';
coeff_alphaH.m     = 0.67849;
coeff_alphaH.c     = -1.95537;


coeff_alphaV.a_vec = [-0.07771 0.56727 -0.20238 -48.2991 48.5833]';
coeff_alphaV.b_vec = [2.33840 0.95545 1.14520 0.791669 0.791459]';
coeff_alphaV.c_vec = [-0.76284 0.54039 0.26809 0.116226 0.116479]';
coeff_alphaV.m     = -0.053739;
coeff_alphaV.c     = 0.83433;


coeffs = coeff_kH;
kH = 10 ^ ( sum( coeffs.a_vec .* exp( -((log10(f)-coeffs.b_vec)./coeffs.c_vec).^2 ) ) + coeffs.m*log10(f) + coeffs.c);

coeffs = coeff_kV;
kV = 10 ^ (sum( coeffs.a_vec .* exp( -((log10(f)-coeffs.b_vec)./coeffs.c_vec).^2 ) ) + coeffs.m*log10(f) + coeffs.c);

coeffs = coeff_alphaH;
alphaH = sum( coeffs.a_vec .* exp( -((log10(f)-coeffs.b_vec)./coeffs.c_vec).^2 ) ) + coeffs.m*log10(f) + coeffs.c;

coeffs = coeff_alphaV;
alphaV = sum( coeffs.a_vec .* exp( -((log10(f)-coeffs.b_vec)./coeffs.c_vec).^2 ) ) + coeffs.m*log10(f) + coeffs.c;


end

