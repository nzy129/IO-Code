function f=inverse2(H,i,n)

fun=@(x)(x.^(i-1)).*((1-x).^(n-i));

Huge=@(z)(factorial(n)/(factorial(n-i))/(factorial(i-1))*integral(fun,0,z)-H);
f=fzero(Huge,[0,1]);
