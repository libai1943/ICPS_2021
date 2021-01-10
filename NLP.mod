param Nfe := 100;
param Nv;
var tf >= 0.1;
var dt = tf / (Nfe - 1);
param Nobs;
param OC{i in {1..Nobs}, k in {1..2}};
param TPBV{i in {1..Nv}, j in {1..8}};
param M_V2V{i in {1..Nv}, j in {1..Nv}, a in {1..4}, b in {1..4}, k in {1..Nfe}};
param M_V2O{i in {1..Nv}, j in {1..Nobs}, b in {1..4}, k in {1..Nfe}};
param M_V2S{i in {1..Nv}, a in {1..4}, b in {1..4}, k in {1..Nfe}};
param M{j in {1..3}};

param amax := 0.5;
param jmax := 0.1;
param vmax := 3.0;
param wmax := 0.3;
param phymax := 0.7;

param L_tractor_front_hang := 0.25;
param L_tractor_wheelbase := 1.5;
param L_tractor_rear_hang := 0.25;
param LHW := 1;
param L := 3.0;
param L_trailer_front_hang := 1;
param L_trailer_rear_hang := 1;

##### Decision variabes besides tf #####
var x{i in {1..Nv}, j in {1..4}, k in {1..Nfe}};
var y{i in {1..Nv}, j in {1..4}, k in {1..Nfe}};
var v{i in {1..Nv}, j in {1..4}, k in {1..Nfe}};
var theta{i in {1..Nv}, j in {1..4}, k in {1..Nfe}};
var xc{i in {1..Nv}, j in {1..4}, k in {1..Nfe}};
var yc{i in {1..Nv}, j in {1..4}, k in {1..Nfe}};
var phy{i in {1..Nv}, k in {1..Nfe}};


minimize objective_:
tf;

s.t. timer:
tf <= 0.50 * Nfe;

s.t. DIFF_dxdt {i in {1..Nv}, j in {2..Nfe}}:
x[i,1,j] = x[i,1,j-1] + dt * v[i,1,j-1] * cos(theta[i,1,j-1]);

s.t. DIFF_dydt {i in {1..Nv}, j in {2..Nfe}}:
y[i,1,j] = y[i,1,j-1] + dt * v[i,1,j-1] * sin(theta[i,1,j-1]);

s.t. DIFF_dtheta1dt {i in {1..Nv}, j in {2..Nfe}}:
theta[i,1,j] = theta[i,1,j-1] + dt * v[i,1,j-1] * tan(phy[i,j-1]) / L_tractor_wheelbase;

s.t. ALGE_x_2_to_Nv {i in {1..Nv}, j in {2..4}, k in {1..Nfe}}:
x[i,j,k] = x[i,j-1,k] - L * cos(theta[i,j,k]) - M[j-1] * cos(theta[i,j-1,k]);

s.t. ALGE_y_2_to_Nv {i in {1..Nv}, j in {2..4}, k in {1..Nfe}}:
y[i,j,k] = y[i,j-1,k] - L * sin(theta[i,j,k]) - M[j-1] * sin(theta[i,j-1,k]);

s.t. DIFF_theta_2_to_Nv {i in {2..Nfe}, j in {2..4}, kk in {1..Nv}}:
L * (theta[kk,j,i] - theta[kk,j,i-1]) = dt * (v[kk,j-1,i] * sin(theta[kk,j-1,i] - theta[kk,j,i])) - M[j-1] * cos(theta[kk,j-1,i] - theta[kk,j,i]) * (theta[kk,j-1,i] - theta[kk,j-1,i-1]);

s.t. DIFF_v_2_to_Nv {i in {2..Nfe}, j in {2..4}, kk in {1..Nv}}:
dt * v[kk,j,i] = dt * v[kk,j-1,i] * cos(theta[kk,j-1,i] - theta[kk,j,i]) + M[j-1] * sin(theta[kk,j-1,i] - theta[kk,j,i]) * (theta[kk,j-1,i] - theta[kk,j-1,i-1]);

s.t. SpecifyTractorXc1 {i in {1..Nv}, j in {1..Nfe}}:
xc[i,1,j] = x[i,1,j] + cos(theta[i,1,j]) * 0.75;

s.t. SpecifyTractorYc1 {i in {1..Nv}, j in {1..Nfe}}:
yc[i,1,j] = y[i,1,j] + sin(theta[i,1,j]) * 0.75;

s.t. SpecifyTractorXc2To4 {i in {1..Nv}, j in {1..Nfe}, kk in {2..4}}:
xc[i,kk,j] = x[i,kk,j];

s.t. SpecifyTractorYc2To4 {i in {1..Nv}, j in {1..Nfe}, kk in {2..4}}:
yc[i,kk,j] = y[i,kk,j];

s.t. EQ_init_x {i in {1..Nv}} :
x[i,1,1] = TPBV[i,1];

s.t. EQ_init_y {i in {1..Nv}} :
y[i,1,1] = TPBV[i,2];

s.t. EQ_init_theta1to4 {i in {1..Nv}, kk in {1..4}} :
theta[i,kk,1] = TPBV[i,kk+2];

s.t. EQ_bv_phy {i in {1..Nv}, index in {1,Nfe}} :
phy[i,index] = 0;

s.t. EQ_bv_v {i in {1..Nv}, index in {1,2,3,Nfe-2,Nfe-1,Nfe}} :
v[i,1,index] = 0;

s.t. EQ_terminal_x {i in {1..Nv}, j in {1..4}} :
TPBV[i,7] - 11 <= xc[i,j,Nfe] <= TPBV[i,7];

s.t. EQ_terminal_y {i in {1..Nv}, j in {1..4}} :
TPBV[i,8] - 0.1 <= yc[i,j,Nfe] <= TPBV[i,8] + 0.1;

s.t. Bonds_v1 {i in {1..Nfe}, kk in {1..Nv}}:
-vmax <= v[kk,1,i] <= vmax;

s.t. Bonds_a1 {i in {2..Nfe}, kk in {1..Nv}}:
v[kk,1,i] - v[kk,1,i-1] <= amax * dt;

s.t. Bonds_a2 {i in {2..Nfe}, kk in {1..Nv}}:
v[kk,1,i] - v[kk,1,i-1] >= -amax * dt;

s.t. Bonds_jerk1 {i in {3..Nfe}, kk in {1..Nv}}:
(v[kk,1,i] - 2 * v[kk,1,i-1] + v[kk,1,i-2]) >= -jmax * (dt^2);

s.t. Bonds_jerk2 {i in {3..Nfe}, kk in {1..Nv}}:
(v[kk,1,i] - 2 * v[kk,1,i-1] + v[kk,1,i-2]) <= jmax * (dt^2);

s.t. Bonds_phy {i in {1..Nfe}, kk in {1..Nv}}:
-phymax <= phy[kk,i] <= phymax;

s.t. Bonds_w1 {i in {2..Nfe}, kk in {1..Nv}}:
phy[kk,i] - phy[kk,i-1] <= wmax * dt;

s.t. Bonds_w2 {i in {2..Nfe}, kk in {1..Nv}}:
phy[kk,i] - phy[kk,i-1] >= -wmax * dt;

s.t. Bound_x {i in {1..Nv}, j in {1..Nobs}, k in {1..4}, m in {1..Nfe}}:
-20 + 1.414 <= xc[i,k,m] <= 20 - 1.414;

s.t. Bound_y {i in {1..Nv}, j in {1..Nobs}, k in {1..4}, m in {1..Nfe}}:
-20 + 1.414 <= yc[i,k,m] <= 20 - 1.414;

s.t. Bonds_theta {i in {1..Nfe}, kk in {1..Nv}, mm in {1..3}}:
-1.5708 <= theta[kk,mm+1,i] - theta[kk,mm,i] <= 1.5708;

s.t. VehicleToVehicle {i1 in {1..(Nv-1)}, i2 in {(i1+1)..Nv}, j1 in {1..4}, j2 in {1..4}, k in {1..Nfe}}:
if (M_V2V[i1,i2,j1,j2,k] == 1) then (-(xc[i1,j1,k] - xc[i2,j2,k])^2 - (yc[i1,j1,k] - yc[i2,j2,k])^2 + 8) else (0) <= 0;

s.t. VehicleToObstacle {i in {1..Nv}, j in {1..Nobs}, k in {1..4}, m in {1..Nfe}}:
if (M_V2O[i,j,k,m] == 1) then (-(xc[i,k,m] - OC[j,1])^2 - (yc[i,k,m] - OC[j,2])^2 + 8) else (0) <= 0;

s.t. VehicleToSelf {i in {1..Nv}, a in {1..3}, b in {(a+1)..4}, m in {1..Nfe}}:
if (M_V2S[i,a,b,m] == 1) then (-(xc[i,a,m] - xc[i,b,m])^2 - (yc[i,a,m] - yc[i,b,m])^2 + 8) else (0) <= 0;

data;
param Nv := include Nv;
param Nobs := include Nobs;
param: OC := include OC;
param: TPBV := include TPBV;
param: M := include M;
param: M_V2O := include M_V2O;
param: M_V2V := include M_V2V;
param: M_V2S := include M_V2S;