function disp = CylinderAnalysis(point, a, b, p, GInf, G1, nu, t, tao)
x = point(1);
y = point(2);
r = sqrt(x^2+y^2);
K = (GInf+G1)*2*(1+nu)/(3*(1-2*nu));

% G = G_Inf + G1*exp(-t/tao)
% a<b
p1G = tao;
q0G = 2*GInf;
q1G = 2*(GInf+G1)*tao;
p0G = 1;
disp = zeros(length(t),2);
f = zeros(length(t),1);
for i=1:length(t)
    f(i)=a^2*p/(b^2-a^2)*(3*r*(1/(6*K + q0G)+ ( p1G/(q1G + 6*K*p1G)-1/(6*K + q0G) )*exp(-(6*K + q0G)/...
        (q1G + 6*K*p1G)*t(i))) + b^2/(q0G*r)*(1+(q0G/q1G*p1G-1)*exp(-q0G/q1G*t(i))));
    disp(i,1) = -f(i)*x/r;
    disp(i,2) = -f(i)*y/r;
end
end