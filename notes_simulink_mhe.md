TODO MHE SIMULINK

obj.problem_class = "OCP" SIM & OCP
add MHE

when rendering s function template:
if MHE:
    Inputs:
        for arrival cost:
        - yref at stage 0
        - W at stage 0
    Outputs:
        - x[end] or x trajectory
        - u trajectory (process noise) - w

    Other setter:
        - remove yref for not initial stage


# Git:
- git add interfaces/acados_template
- git commit -m "christian template edits" *
- git checkout -b christian
- git fetch origin
- git checkout origin master
- git cherry-pick <commit-id>(*)

#
Test templated solver in this style (i.e. without Simulink)
https://github.com/acados/acados/blob/6d59a6f9692a3f531585d412f11af29980bca403/examples/acados_matlab_octave/test/test_template_pendulum_exact_hess.m#L215-L219
