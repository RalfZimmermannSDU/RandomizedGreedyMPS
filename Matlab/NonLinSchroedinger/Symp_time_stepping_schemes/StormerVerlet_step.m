function [qn_p1, pn_p1] = StormerVerlet_step(dt,qn,pn,grad_qH, grad_pH)
%
% one time step step of size dt according to the St"ormer-Verlet scheme
% 
% qn+½ = qn   + dt/2 grad_pH(qn+½, pn)
% pn+1 = pn   - dt/2(grad_qH(qn+½, pn) + grad_qH(qn+½, pn+1)) 
% qn+1 = qn+½ + dt/2 grad_pH(qn+½, pn+1)
% Alternative ?
% pn+½ = pn   - dt/2 grad_qH(qn, pn+½)
% qn+1 = qn   + dt/2(grad_pH(qn, pn+½) + grad_pH(qn+1, pn+½)) 
% pn+1 = pn+½ - dt/2 grad_qH(qn+1, pn+½)
%
%
% inputs:
% dt : time step
% qn : q-state vector at step n (may be a vector)
% pn : q-state vector at step n (may be a vector)
% right-hand side functions 
%  grad_qH(q,p), grad_pH(q,p)
%
% outputs: updated states
% qn_p1      % read q_{n+1}
% pn_p1      % read p_{n+1}


% suppress output onscreen when using fsolve
options = optimoptions('fsolve','Display','none');
%-----------------------------------------
% 1 step: compute implicitely defined qn+½ 
% solve associated 	NONLINEAR EQUATION
%-----------------------------------------

fun0 = @(s) qn   + (dt/2)*grad_pH(s, pn) - s;

%
% zero of fun0 is solution to 0 =  qn   + dt/2 grad_pH(s, pn) - s
% use qn as initial guess
qn_half = fsolve(fun0, qn, options);


%-----------------------------------------
% 2 step: compute implicitely defined pn+1
% solve associated 	NONLINEAR EQUATION
%-----------------------------------------
% pn+1 = pn   - dt/2(grad_qH(qn+½, pn) + grad_qH(qn+½, pn+1))
grad_temp = grad_qH(qn_half, pn);

fun0 = @(s) pn - (dt/2)*(grad_temp + grad_qH(qn_half, s)) - s;

% zero of fun0 is solution to 
% 0 = pn   - dt/2(grad_qH(qn+½, pn) + grad_qH(qn+½, s)) - s
% use pn as initial guess
pn_p1 = fsolve(fun0, pn, options);

%-----------------------------------------
% 3 step: compute explicitely defined qn+1 
%-----------------------------------------
qn_p1 = qn_half + dt/2*grad_pH(qn_half, pn_p1);

return;
end


