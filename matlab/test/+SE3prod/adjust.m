function coeffs = adjust( coeffs, eigvals, t )

coeffs(1) = 0.5*coeffs(1);
coeffs(4) = 0.5*coeffs(4);
coeffs(5) = 0.5*coeffs(5);
coeffs(6) = 0.5*coeffs(6);

for i = 2:length(coeffs)-3;
    
    if t > 0.000001
        coeffs(i) = coeffs(i)*sin(sqrt(eigvals(i))*t/2)/(sin(sqrt(eigvals(i))*t) + 0.000000000001);
    else
        coeffs(i) = 0.5*coeffs(i);
    end

end

end