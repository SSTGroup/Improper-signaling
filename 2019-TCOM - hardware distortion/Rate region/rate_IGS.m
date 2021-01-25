function r=rate_IGS(sigma2,p,q,a,b,f,g,f_til)

temp=((sigma2^2+p*a)^2-abs(q*f+p*f_til)^2) ...
    /((sigma2^2+p*b)^2-abs(q*g+p*f_til)^2);
r=0.5*log2(temp);