function out = soft_threshold(z, threshold)
% Se abs(z) < threshold -> restituisce zero
% altrimenti ritorna z - threshold
    out = sign(z) .* max(abs(z) - threshold, 0);
end