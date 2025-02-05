function mat_f = fourier_expand(pha, order)
mat_f  =  zeros(size(pha,1),order);
for i = 1:order
    mat_f(:,2*i-1) = cos(i*pha);
    mat_f(:,2*i)   = sin(i*pha);
end
end