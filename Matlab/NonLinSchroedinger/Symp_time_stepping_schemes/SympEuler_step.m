
function [qn_p1, pn_p1] =SympEuler_step(dt,qn,pn,grad_qH, grad_pH)
%
% one time step step of size dt according to the 
% symplectic Euler scheme
% 
% 
% pn+1 = pn - dt*grad_qH(qn, pn+1) 
% qn+1 = qn + dt*grad_pH(qn, pn+1)
% Alternative ?
% pn+1 = pn - dt*grad_qH(qn+1, pn) 
% qn+1 = qn + dt*grad_pH(qn+1, pn)
%
% inputs:
% dt : time step
% qn : q-state vector at step n (may be a vector)
% pn : q-state vector at step n (may be a vector)
% right-hand side functions 
%  grad_qH(q,p), grad_pH(q,p)
%
% outputs: updated state
% qn_p1 % read q_{n+1}
% pn_p1 % read p_{n+1}


% suppress output onscreen when using fsolve
options = optimoptions('fsolve','Display','none');
%-----------------------------------------
% 1 step: compute implicitely defined qn+½ 
% solve associated 	NONLINEAR EQUATION
%-----------------------------------------

fun0 = @(s) pn   - dt*grad_qH(qn, s) - s;

%
% zero of fun0 is solution to 0 =  qn   + dt/2 grad_pH(s, pn) - s
% use qn as initial guess
pn_p1 = fsolve(fun0, pn, options);


%-----------------------------------------
% 2 step: compute explicitely defined qn+1 
%-----------------------------------------
qn_p1 = qn + dt*grad_pH(qn, pn_p1);

return;
end


