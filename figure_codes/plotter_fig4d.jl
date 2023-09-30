#===========================
PLOTTER FILE FOR 
FIGURE 4D
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations, Loess
using HypothesisTests, LsqFit, EffectSizes
CairoMakie.activate!(type="svg")
size_cm = (5.0, 4.1);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

# Read the data file
df_data = DataFrame(CSV.File("../figure_data_files/master_data_RNA_DNA_p27_kE0.csv"))

# Read the parameter files
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df_param        = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_1.xlsx","Sheet1";infer_eltypes=true))
# Parameter order (Total 11): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, T0_exp
df_param_sans_S = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_1_sans_S.xlsx","Sheet1";infer_eltypes=true))

macaque_names = ["29915", "29925", "31041", "AV979",
"BA081", "BA209", "BB598", "BC094",
"BC179", "BC657", "BD536", "BD885",
"BG927", "BL669", "BO186", "BO413"]
# Progressors are placed in the top row, while controllers are placed in the remaining slots
ordering_idx = [3, 4, 1, 2]; # Order in which subplots must be placed

#=
Section 2: Define the required storage variables
=#

time_vec = 0:0.01:365; tspan = (0, 365);
N = 1e6; C = 1e6;
sol_V        = zeros(16,length(time_vec));
sol_V_sans_S = zeros(16,length(time_vec));
sol_E        = zeros(16,length(time_vec));
sol_E_sans_S = zeros(16,length(time_vec));
sol_ψ        = zeros(16,length(time_vec));
sol_ψ_sans_S = zeros(16,length(time_vec));
sol_S        = zeros(16,length(time_vec));
sol_S_sans_S = zeros(16,length(time_vec));

#=
Section 3: Solve the model for each macaque and store the values
=#

# Best-fit model

for i = 1:1:16
    # Section 3a: Get the parameter and initial condition vectors
    curr_p    = Vector{Float64}(df_param[i,2:13])
    init_vec  = [0.0,(10^(-2.76))/(curr_p[7]),0.0,0.0,0.0];
    curr_T0   = 10^(curr_p[12])
    curr_init = init_vec
    if i == 8 || i == 9 || i == 11 || i == 15
        curr_init[1] = curr_T0
    else
        curr_init[1] = curr_T0
        curr_init[2] = 10*curr_init[2]
    end

    # Section 3b: Solve the model
    prob = ODEProblem(model_1!,curr_init,tspan,curr_p)
    soln = solve(prob,saveat = 0.01,isoutofdomain=(u,p,t)->any(x->x<0,u))

    # Section 3c: Extract the required values and save them
    sol_V[i,:] .= log10.(2*curr_p[7].*max.(soln[2,:], vec((1e-8)*ones(1,length(time_vec)))))
    sol_E[i,:] .= soln[4,:]
    sol_ψ[i,:] .= (soln[4,:].*soln[5,:])
    sol_S[i,:] .= S_func!.(soln[4,:].*soln[5,:])
end

# Model sans S

for i = 1:1:16
    # Section 3a: Get the parameter and initial condition vectors
    curr_p    = Vector{Float64}(df_param_sans_S[i,2:12])
    init_vec  = [0.0,(10^(-2.76))/(curr_p[7]),0.0,0.0];
    curr_T0   = 10^(curr_p[11])
    curr_init = init_vec
    if i == 8 || i == 9 || i == 11 || i == 15
        curr_init[1] = curr_T0
    else
        curr_init[1] = curr_T0
        curr_init[2] = 10*curr_init[2]
    end

    # Section 3b: Solve the model
    prob = ODEProblem(model_1_sans_S!,curr_init,tspan,curr_p)
    soln = solve(prob,saveat = 0.01,isoutofdomain=(u,p,t)->any(x->x<0,u))

    # Section 3c: Extract the required values and save them
    sol_V_sans_S[i,:] .= log10.(2*curr_p[7].*max.(soln[2,:], vec((1e-8)*ones(1,length(time_vec)))))
    sol_E_sans_S[i,:] .= soln[4,:]
    sol_ψ_sans_S[i,:] .= (soln[4,:]) # Because for this model, k does not change with time
    sol_S_sans_S[i,:] .= S_func!.(soln[4,:])
end

# Defining the linear model

@. lin(x, p) = p[1] + p[2] * x;
p0 = [1.0, -1.0];

#=
Section 4: Plot
=#

#================================
FIGURE 4D: CORRELATION PLOT 
OF <AUC> OF S VS SPVL BY MODEL 1
================================#

# Estimate the required variables

AUC_m1_SIC = zeros(12, 1); AUC_m1_VIR = zeros(4, 1);
AUC_m7_SIC = zeros(12, 1); AUC_m7_VIR = zeros(4, 1);
VLf_m1_SIC = zeros(12, 1); VLf_m1_VIR = zeros(4, 1);
VLf_m7_SIC = zeros(12, 1); VLf_m7_VIR = zeros(4, 1);

ordering_idx = [3, 4, 7, 12, 1, 2, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16];
q = 1;
for i = 1:1:4
    idx = ordering_idx[q];
    AUC_m1_VIR[i] = mid_point_int!(time_vec[1:2801],sol_S[idx,1:2801])/28;
    AUC_m7_VIR[i] = mid_point_int!(time_vec[1:2801],sol_S_sans_S[idx,1:2801])/28;
    VLf_m1_VIR[i] = sol_V[idx,end];
    VLf_m7_VIR[i] = sol_V_sans_S[idx,end];
    q += 1;
end

for i = 1:1:12
    idx = ordering_idx[q];
    AUC_m1_SIC[i] = mid_point_int!(time_vec[1:2801],sol_S[idx,1:2801])/28;
    AUC_m7_SIC[i] = mid_point_int!(time_vec[1:2801],sol_S_sans_S[idx,1:2801])/28;
    VLf_m1_SIC[i] = sol_V[idx,end];
    VLf_m7_SIC[i] = sol_V_sans_S[idx,end];
    q += 1;
end

AUC_m1_SIC = vec(AUC_m1_SIC); AUC_m1_VIR = vec(AUC_m1_VIR);
AUC_m7_SIC = vec(AUC_m7_SIC); AUC_m7_VIR = vec(AUC_m7_VIR);
VLf_m1_SIC = vec(VLf_m1_SIC); VLf_m1_VIR = vec(VLf_m1_VIR);
VLf_m7_SIC = vec(VLf_m7_SIC); VLf_m7_VIR = vec(VLf_m7_VIR);
AUC_m1 = [AUC_m1_SIC; AUC_m1_VIR]; VLf_m1 = [VLf_m1_SIC; VLf_m1_VIR];
AUC_m7 = [AUC_m7_SIC; AUC_m7_VIR]; VLf_m7 = [VLf_m7_SIC; VLf_m7_VIR];

fit_4d = curve_fit(lin, AUC_m1, VLf_m1, p0); # Equation is: VLf = slp*AUC + int
int = coef(fit_4d)[1]; slp = coef(fit_4d)[2];
@. fit_model(x) = int + slp*x;

x_lim  = [-0.05, 1.0];
y_lim  = [-3.2, 6.0];
x_tick = [0.0, 0.25, 0.50, 0.75, 1.0];
y_tick = [-2, 0, 2, 4, 6];

fig = Figure(resolution = size_pt, fontsize=12);

# Setting the panel
ax = Axis(fig[1, 1], spinewidth = 1, xticks = x_tick, yticks = y_tick,  xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

# Obtaining the LOESS curve
model  = loess(AUC_m1, VLf_m1);
x_pred = range(extrema(AUC_m1)...; step = 0.5);
y_pred = predict(model, x_pred);

# Plotting commands
scatter!(ax, AUC_m1_SIC, VLf_m1_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(ax, AUC_m1_VIR, VLf_m1_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
# lines!(ax, sort(AUC_m1), fit_model.(sort(AUC_m1)), color = :black, linewidth = 1);
# lines!(ax, x_pred, y_pred, color = :black, linewidth = 1);

xlims!(ax, x_lim);
ylims!(ax, y_lim);
# order += 1;

fig
save("fig_4d.svg", fig, px_per_unit = 6)


