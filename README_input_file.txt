|||||||||||||||||||||||||||||||||||||||| NUMERICAL SETTINGS |||||||||||||||||||||||||||||||||||||||||||||||||||||||

test_case_ion_vacuum.inp

Cells Number
Left boundary position
Right boundary position
Grid type (uniform, Chebyshev)
Simulation type  (Euler -> single fluid, Euler_multi_fluid -> multi fluid)
t_start
t_end
timesteps
Time scheme (MidpointEuler,ForwardEuler)
limiter (0,1)
write every
EM fields activation (esc -> self-consistent electric field along x, force -> impose fixed external electric field along x axis)
Initial condition (Test_case) (plasma_expansion ->  Plasma expansion into vacuum ,discontinous-> Sod Shock problem,  continous-> uniform continous initial condition)
Collision simulation (collision_full -> elastic and inelastic collisional terms through Boltzmann operator ,collisionA -> approximation of elastic terms with thermal velocity (NOT YET IMPLMENTED FOR MULTIFLUID), collision_elastic -> only elastic term derivation through Boltzmann operator NOT IMPLEMENTED FOR ADIMENSIONAL, no_coll-> no collisional terms)
Non-dimensional bool (Dim,Nondim)

Note: BC here are ONLY for fluid variables within time_integrators

|||||||||||||||||||||||||||||||||||||||| INITIAL CONDITION |||||||||||||||||||||||||||||||||||||||||||||||||||||||

initial_conditions.inp

rho_init_electrons
T_init_electrons
u_init_electrons
rho_init_ions
T_init_ions
u_init_ions


|||||||||||||||||||||||||||||||||||||||| BOUNDARY CONDITIONS |||||||||||||||||||||||||||||||||||||||||||||||||||||||

boundary_conditions.inp


BC_Left_electrons	    (open, periodic, Neumann) For periodic both must be set, otherwise it will throw runtime error
rho_left_electrons
T_left_electrons
u_left_electrons
BC_Right_electrons     (open, periodic, Neumann) For periodic both must be set, otherwise it will throw runtime error
rho_right_electrons
T_right_electrons
u_right_electrons
BC_Left_ions	    (open, periodic, Neumann) For periodic both must be set, otherwise it will throw runtime error
rho_left_ions
T_left_ions
u_left_ions
BC_Right_ions		 (open, periodic, Neumann) For periodic both must be set, otherwise it will throw runtime error
rho_right_ions
T_right_ions
u_right_ions

|||||||||||||||||||||||||||||||||||||||| PHYSICAL CONSTANTS |||||||||||||||||||||||||||||||||||||||||||||||||||||||


physical_constants.inp

gamma_electrons
mass_electrons
charge_electrons
k_boltzmann
eps_0
gamma_ions
mass_ions
charge_ions
neutrals_density
neutrals_velocity
threshold energy for Xenon ionization
threshold energy for Xenon excitation
Still ions background density (for single fluid only)


|||||||||||||||||||||||||||||||||||||||| REFERENCE VARIABLES |||||||||||||||||||||||||||||||||||||||||||||||||||||||

FOR NONDIMENSIONAL ANALYSIS

characteristics_quantities.inp

rho_ref_el
rho_ref_ions
x_ref
time_ref
nu_coll_ref
phi_ref
B_mag
E_y
uy_el
uy_ion
P_ref_el
P_ref_ion

|||||||||||||||||||||||||||||||||||||||| ELECTRIC FIELD |||||||||||||||||||||||||||||||||||||||||||||||||||||||

POISSON CONDITIONS

EM_fields.dat

Ex_in
Ey_in
B_in
phi_0
phi_end
BC Left (Poisson) (Dirichlet, Neumann , periodic)
BC_Right (Poisson) (Dirichlet, Neumann , periodic) Periodic has to be there for boths
Neu_left (Neumann non-homogeneous BC)
Neu_right (Neumann non-homogeneous BC)
