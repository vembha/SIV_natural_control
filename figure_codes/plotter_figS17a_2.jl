#===========================
PLOTTER FILE FOR 
FIGURE S17A_2
REVISION VERSION #3
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
using HypothesisTests, LsqFit, EffectSizes
CairoMakie.activate!(type="svg")
size_cm = (5.0, 4.1);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

# Read the data file
df_data = DataFrame(CSV.File("../../figure_data_files/v9_files/master_data_RNA_DNA_p27_kE0.csv"))

# Read the parameter files
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df_param        = DataFrame(XLSX.readtable("../../figure_data_files/v9_files/parameters_model_1.xlsx","Sheet1";infer_eltypes=true))
# Parameter order (Total 11): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, T0_exp
df_param_sans_S = DataFrame(XLSX.readtable("../../figure_data_files/v9_files/parameters_model_1_sans_S.xlsx","Sheet1";infer_eltypes=true))

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
sol_E        = zeros(16,length(time_vec));
sol_ψ        = zeros(16,length(time_vec));
sol_S        = zeros(16,length(time_vec));

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

#=
Section 4: Plot
=#

# Estimate the required variables

AUC_m1_SIC = zeros(12, 4); AUC_m1_VIR = zeros(4, 4);
VLf_m1_SIC = zeros(12, 1); VLf_m1_VIR = zeros(4, 1);

ordering_idx = [3, 4, 7, 12, 1, 2, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16];
t_end = [14, 28, 42, 90]; # Cutoff time points
idx_end = [1401, 2801, 4201, 9001]; # Cutoff indices
for k = 1:1:4 # Different cutoffs
    q = 1;

    # VIRs
    for i = 1:1:4 # Different macaques
        idx = ordering_idx[q];
        AUC_m1_VIR[i, k] = mid_point_int!(time_vec[1:idx_end[k]],sol_S[idx,1:idx_end[k]])/t_end[k];
        VLf_m1_VIR[i] = sol_V[idx,end];
        q += 1;
    end

    # SICs
    for i = 1:1:12
        idx = ordering_idx[q];
        AUC_m1_SIC[i, k] = mid_point_int!(time_vec[1:idx_end[k]],sol_S[idx,1:idx_end[k]])/t_end[k];
        VLf_m1_SIC[i] = sol_V[idx,end];
        q += 1;
    end
end


# AUC_m1_SIC = vec(AUC_m1_SIC); AUC_m1_VIR = vec(AUC_m1_VIR);
VLf_m1_SIC = vec(VLf_m1_SIC); VLf_m1_VIR = vec(VLf_m1_VIR);
AUC_m1 = [AUC_m1_SIC; AUC_m1_VIR]; VLf_m1 = [VLf_m1_SIC; VLf_m1_VIR];

# # Write the data to .xlsx file

XLSX.openxlsx("fig_S17a_1_data.xlsx", mode = "rw") do xf
    sheet = xf[1]
    XLSX.writetable!(sheet, [VLf_m1, AUC_m1[:,1], AUC_m1[:,2], AUC_m1[:,3], AUC_m1[:,4]], ["spVL", "tAUC_14", "tAUC_28", "tAUC_42", "tAUC_90"])
end

# Plot

x_lim  = [0.25, 4.75];
y_lim  = [0.3, 0.6];
x_tick = [1, 2, 3, 4];
y_tick = [0.3, 0.4, 0.5, 0.6];

fig = Figure(resolution = size_pt, fontsize=12);

x = [1, 2, 3, 4];
y = [0.39, 0.51, 0.50, 0.36];

# Setting the panel
ax = Axis(fig[1, 1], spinewidth = 1, xticks = (x_tick, ["14", "28", "42", "90"]), yticks = y_tick,  xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

# Plotting commands
barplot!(ax, x, y, color = :steelblue3, strokecolor = :black, strokewidth = 1);

xlims!(ax, x_lim);
ylims!(ax, y_lim);

fig
save("fig_S17a_2.svg", fig, px_per_unit = 6)

