function f = pitchfork(z,mu,gamma)

f(1,1) = mu.*z(1)-z(1).^3+gamma;

end