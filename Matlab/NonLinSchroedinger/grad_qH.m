function [grad_qHqp, a_nonlin] = grad_qH(q,p,A, para)
% evaluation of the nonlinear right-hand side function
% grad_qH(q,p), depending on a scalar parameter
%
% A is the finite differences operator
a_nonlin  = (q.*q + p.*p).*q;
grad_qHqp = A*q + para*a_nonlin;

return;
end

