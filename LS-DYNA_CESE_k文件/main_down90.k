$# LS-DYNA Keyword file created by LS-PrePost(R) 2024/R2(4.11.5)-19May2024
$# Created on Jun-5-2024 (17:38:17)
*KEYWORD
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
     0.005         0       0.0       0.01.000000E8         0
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid
5.00000E-5         0         0         0         0
*CESE_PART
$#     pid       mid     eosid
         1         0         3
*CESE_DATABASE_ELOUT
$#   outlv     dtout
         15.00000E-5
$#   elsid
         1
*CESE_CONTROL_TIMESTEP
$#    iddt       cfl     dtint
         0       0.51.30274E-7
*CESE_CONTROL_SOLVER
$#   icese     iflow     igeom    iframe
         0         1         2         0
*CESE_CONTROL_LIMITER
$#   idlmt      alfa     beta       epsr
         0       4.0       1.0       0.5
*CESE_BOUNDARY_SOLID_WALL_SET
$#    ssid      lcid        vx        vy        vz
         2         0       0.0       0.0       0.0
         3         0       0.0       0.0       0.0
*CESE_BOUNDARY_NON_REFLECTIVE_SET
$#    ssid
         4
*CESE_INITIAL
$#       u         v         w       rho         p         t
       0.0       0.0       0.0     1.225  101325.0       0.0
*CESE_BOUNDARY_PRESCRIBED_SET
$#    ssid    idcomp
         1         0
$#    lc_u     lc_v       lc_w    lc_rho     lc_p       lc_t
         0         0         0         0         0        -1
$#    sf_u     sf_v       sf_w    sf_rho     sf_p       sf_t
       0.0       0.0       0.0     1.225   6079500       0.0
*CESE_EOS_IDEAL_GAS
$#   eosid        cv        cp
         3     717.5    1004.5
*INCLUDE
mesh_down90.k
*END
