# AMPL program to optimize a quarter-plane symmetric APLC
# Created by ntz with scda.py on zimmer.stsci.edu at 2016-03-25 14:23
load amplgsl.dll;

#---------------------

param pi:= 4*atan(1);

#---------------------
param c := 10.00;

#---------------------
param Rmask := 4.000;
param rho0 := 3.50;
param rho1 := 10.00;

#---------------------
param N := 250;				# discretization parameter (pupil)
param M := 50;				# discretization parameter (mask)

param Nimg := 22;           # discretization parameter (image)
param rho2 := 10.00;        # NTZ: We paramaterize our image plane resolution by fpres = sampling rate at 
                            #      the shortest wavelength. Then Nimg is an integer function of fpres, oca,
#---------------------      #      and bw. This integer is not specified directly by the user, but computed "privately"
param bw := 0.10;           #      by the APLC class constructor.
param lam0 := 1.;
param dl := bw*lam0;
param Nlam := 3;

#---------------------

#---------------------
# steps in each plane
param dx := 1/(2*N);
param dy := dx;

param dmx := Rmask/M;
param dmy := dmx;

param dxi := rho2/Nimg;
param deta := dxi;

#---------------------
# coordinate vectors in each plane
set Xs := setof {i in 0.5..N-0.5 by 1} i*dx;
set Ys := setof {j in 0.5..N-0.5 by 1} j*dy;

set MXs := setof {i in 0.5..M-0.5 by 1} i*dmx;
set MYs := setof {j in 0.5..M-0.5 by 1} j*dmy;

set Xis := setof {i in 0..Nimg-1 by 1} i*dxi;
set Etas := setof {j in 0..Nimg-1 by 1} j*deta;
#---------------------
# Load telescope aperture
param TelAp {x in Xs, y in Ys};

read {y in Ys, x in Xs} TelAp[x,y] < "../InputMasks/TelAp/TelAp_quart_hex3X025c1_N0250.dat";
close "../InputMasks/TelAp/TelAp_quart_hex3X025c1_N0250.dat";

# Load FPM
param FPM {mx in MXs, my in MYs};

read {my in MYs, mx in MXs} FPM[mx,my] < "../InputMasks/FPM/FPM_quart_occspot_M050.dat"; 
close "../InputMasks/FPM/FPM_quart_occspot_M050.dat";

# Load Lyot stop
param LS {x in Xs, y in Ys};

read {y in Ys,x in Xs} LS[x,y] < "../InputMasks/LS/LS_quart_ann30D85_X025c1Pad05_N0250.dat";
close "../InputMasks/LS/LS_quart_ann30D85_X025c1Pad05_N0250.dat";

set Ls := setof {l in 1..Nlam} lam0*(1 - bw/2 + (l-1)*bw/(Nlam-1));

#---------------------

set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] > 0} (x,y);
set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
set Lyot := setof {x in Xs, y in Ys: LS[x,y] > 0} (x,y);

param TR := sum {(x,y) in Pupil} TelAp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.

var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;

#---------------------

set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);

#---------------------
var EBm_real_X {mx in MXs, y in Ys, lam in Ls};
var EBm_real {mx in MXs, my in MYs, lam in Ls};

subject to st_EBm_real_X {mx in MXs, y in Ys, lam in Ls}: EBm_real_X[mx,y,lam] = 2.*sum {x in Xs: (x,y) in Pupil} TelAp[x,y]*A[x,y]*cos(2.*pi*x*mx*(lam0/lam))*dx;
subject to st_EBm_real {(mx, my) in Mask, lam in Ls}: EBm_real[mx,my,lam] = 2.*(lam0/lam)*sum {y in Ys} EBm_real_X[mx,y,lam]*cos(2.*pi*y*my*(lam0/lam))*dy;

#---------------------
var ECm_real_X {x in Xs, my in MYs, lam in Ls};
var ECm_real {x in Xs, y in Ys, lam in Ls};

subject to st_ECm_real_X {x in Xs, my in MYs, lam in Ls}: ECm_real_X[x,my,lam] = 2.*sum {mx in MXs: (mx,my) in Mask} EBm_real[mx,my,lam]*cos(2.*pi*x*mx*(lam0/lam))*dmx;
subject to st_ECm_real {(x,y) in Lyot, lam in Ls}: ECm_real[x,y,lam] = 2.*(lam0/lam)*sum {my in MYs} ECm_real_X[x,my,lam]*cos(2.*pi*y*my*(lam0/lam))*dmy;

#---------------------
var ED_real_X {xi in Xis, y in Ys, lam in Ls};
var ED_real {xi in Xis, eta in Etas, lam in Ls};

subject to st_ED_real_X {xi in Xis, y in Ys, lam in Ls}: ED_real_X[xi,y,lam] = 2.*sum {x in Xs: (x,y) in Lyot} (TelAp[x,y]*A[x,y]-ECm_real[x,y,lam])*cos(2.*pi*x*xi*(lam0/lam))*dx;
subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] = 2.*(lam0/lam)*sum {y in Ys} ED_real_X[xi,y,lam]*cos(2.*pi*y*eta*(lam0/lam))*dy;

#---------------------

var ED00_real := 0.0;
subject to st_ED00_real: ED00_real = 4.*sum {x in Xs, y in Ys: (x,y) in Lyot} (A[x,y]*TelAp[x,y])*dx*dy;

#---------------------
maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;

subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/sqrt(2.); 
subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/sqrt(2.);

option solver gurobi;

option gurobi_options "outlev=1 lpmethod=2 crossover=0";

solve;

display solve_result_num, solve_result;
display ED00_real; 

#---------------------

param A_fin {y in Ys, x in Xs};
let {y in Ys, x in Xs} A_fin[x,y] := 0;
let {(x,y) in Pupil} A_fin[x,y] := A[x,y];

printf {y in Ys, x in Xs}: "%15g %15g %15g \n", x, y, A_fin[x,y] > "./nccs_survey_test/solutions/ApodSol_APLC_quart_hex3X025c1_N0250_FPM400M050_LSann30D85X025c1Pad05_Img100C_35WA100CA100_BW10Nlam03fpres2_linbarpre1.dat";
