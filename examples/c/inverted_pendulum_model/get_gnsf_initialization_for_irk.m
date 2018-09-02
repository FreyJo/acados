import casadi.*

x0 = [  1.0000; -5.0000;  0.1000; -0.5000;  0.1000;  1.0000; ];
u0 = 1;


x1_0 = x0(1:gnsf.nx1);
xdot1_z1_0 = (gnsf.E\(gnsf.A*x1_0+ gnsf.B * u0 + gnsf.c));
xdot1_0 = xdot1_z1_0(1:gnsf.nx1);
z1_0 = xdot1_z1_0(1+gnsf.nx1:end);


x1 = gnsf.x(1:gnsf.nx1);
x1dot = gnsf.xdot(1:gnsf.nx1);
z1 = gnsf.z(1:gnsf.nz1);
u = gnsf.u;
f_lo_fun = Function(['f_lo_fun_jac_x1k1uz'], {x1, x1dot, z1, u}, ...
    {gnsf.f_lo_expr});

xdot2_z2_0 = gnsf.E_LO \ ( gnsf.A_LO * x0(gnsf.nx1+1:end) + ...
    f_lo_fun( x1_0, xdot1_0, z1_0, u0))