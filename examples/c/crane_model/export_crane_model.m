clc;
clear all;
close all;

% addpath('../../external/casadi-octave-v3.2.2')
import casadi.*

%% Parameters (taken from Riens ACADO model)
tau1 = 0.012790605943772;   a1   = 0.047418203070092;
tau2 = 0.024695192379264;   a2   = 0.034087337273386;
g = 9.81;

%% Set up States & Controls
xC = SX.sym('xC');     %States  
vC = SX.sym('vC');
xL = SX.sym('xL');     
vL = SX.sym('vL');
theta = SX.sym('theta');
omega = SX.sym('omega'); 
uC = SX.sym('uC');  
uL = SX.sym('uL');
uCR = SX.sym('uCR');  %Controls
uLR = SX.sym('uLR');

x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega);

f_expl = vertcat(vC, ...
                  - 1/tau1 * (vC - a1 * uC), ...
                  vL,...
                  - 1/tau2 * (vL - a2 * uL), ...
                  uCR,...
                  uLR,...
                  omega, ...
                  - (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL);


nx = length(x);
u = vertcat(uCR, uLR);
nu = length(u);

odeFun = Function('odeFun',{x,u},{f_expl});

Sx = SX.sym('Sx',nx,nx);
Sp = SX.sym('Sp',nx,nu);
lambdaX = SX.sym('lambdaX',nx,1);

% Derive Variational Differential Equations
vdeX = SX.zeros(nx,nx);
vdeX = vdeX + jtimes(f_expl,x,Sx);
 
vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u);
vdeP = vdeP + jtimes(f_expl,x,Sp);

vdeFun = Function('vdeFun',{x,Sx,Sp,u},{f_expl,vdeX,vdeP});

jacX = SX.zeros(nx,nx) + jacobian(f_expl,x);
jacFun = Function('jacFun',{x,u},{f_expl,jacX});

% oj: The jtimes function optionally calculates the
% transposed-Jacobian-times-vector product, i.e.  reverse mode AD
%adj = jtimes(f_expl,[x;u],lambdaX,true);
adj = jtimes(f_expl,[u;x],lambdaX,true);
% 
adjFun = Function('adjFun',{x,lambdaX,u},{adj});
% 
S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
hess = S_forw.'*jtimes(adj,[x;u],S_forw);
hess2 = [];
for j = 1:nx+nu
     for i = j:nx+nu
         hess2 = [hess2; hess(i,j)];
     end
end
 
hessFun = Function('hessFun',{x,Sx,Sp,lambdaX,u},{adj,hess2});
% 
opts = struct('mex', false);
odeFun.generate(['ode_model'], opts);
vdeFun.generate(['vde_forw_model'], opts);
jacFun.generate(['jac_model'], opts);
adjFun.generate(['vde_adj_model'], opts);
hessFun.generate(['vde_hess_model'], opts);

% implicit fcn generation for impl integrators
x_dot = SX.sym('x_dot',nx,1);         
f_impl = SX.zeros(nx,1)+(x_dot - f_expl);

impl_odeFun = Function('impl_odeFun',{x,x_dot,u},{f_impl});
jac_x = SX.zeros(nx,nx) + jacobian(f_impl,x);
jac_xdot = SX.zeros(nx,nx) + jacobian(f_impl,x_dot);
jac_u = SX.zeros(nx,nu) + jacobian(f_impl,u);

impl_jacFun_x = Function('impl_jacFun_x',{x,x_dot,u},{jac_x});
impl_jacFun_xdot = Function('impl_jacFun_xdot',{x,x_dot,u},{jac_xdot});
impl_jacFun_u = Function('impl_jacFun_u',{x,x_dot,u},{jac_u});

opts = struct('mex', false);
% 
impl_odeFun.generate(['impl_ode'],opts);
impl_jacFun_x.generate(['impl_jac_x'],opts);
impl_jacFun_xdot.generate(['impl_jac_xdot'],opts);
impl_jacFun_u.generate(['impl_jac_u'],opts);

x0 = zeros(nx,1); x0(3)=0.8;
u0 = [40.108149413030752; -50.446662212534974];
k0 = 0*ones(nx,1);

% odeFun(x0,u0)

full(impl_odeFun(x0,k0,u0))
full(impl_jacFun_x(x0,k0,u0))
full(impl_jacFun_xdot(x0,k0,u0))
full(impl_jacFun_u(x0,k0,u0))