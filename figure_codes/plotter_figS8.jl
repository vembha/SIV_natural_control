#===========================
PLOTTER FILE FOR 
FIGURES S8
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
using LsqFit, HypothesisTests, Distributions
using Colors, ColorSchemes
CairoMakie.activate!(type="svg")
size_cm = (18, 10.84);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

# Read the data file
df_data = DataFrame(CSV.File("../figure_data_files/master_data_RNA_DNA_p27_kE0.csv"))

# Read the parameter file
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df_param = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_1.xlsx","Sheet1";infer_eltypes=true))

macaque_names = ["29915", "29925", "31041", "AV979",
"BA081", "BA209", "BB598", "BC094",
"BC179", "BC657", "BD536", "BD885",
"BG927", "BL669", "BO186", "BO413"]
# Progressors are placed in the top row, while controllers are placed in the remaining slots
ordering_idx = [3, 4, 7, 12, 1, 2, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16]; # Order in which subplots must be placed

#=
Section 2: Read the required parameters
=#

# Read the parameter files
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df_param_SIC = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_1.xlsx","Controllers";infer_eltypes=true))
df_param_VIR = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_1.xlsx","Progressors";infer_eltypes=true))

λ_SIC  = df_param_SIC.lambda_mode;
λ_VIR  = df_param_VIR.lambda_mode;

dD_SIC = df_param_SIC.dD_mode;
dD_VIR = df_param_VIR.dD_mode;

r_SIC = df_param_SIC.r_mode;
r_VIR = df_param_VIR.r_mode;

fD_SIC = df_param_SIC.fD_mode;
fD_VIR = df_param_VIR.fD_mode;

ω_SIC  = df_param_SIC.omega_exp_mode;
ω_VIR  = df_param_VIR.omega_exp_mode;

#=
Section 3: Plot
=#

fig = Figure(resolution = size_pt, fontsize=12);

# Section 3a: λ

x1 = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 12));
x2 = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 4));
x1_mean = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
x2_mean = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
ax = Axis(fig[1, 1], xticks = ([0.75, 2.25], ["Controllers", "Progressors"]), spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

scatter!(x1_mean, fill(median(λ_SIC), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x2_mean, fill(median(λ_VIR), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x1, λ_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(x2, λ_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
xlims!(ax, [0, 3]);

MannWhitneyUTest(λ_SIC, λ_VIR)

# Section 3b: dD

x1 = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 12));
x2 = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 4));
x1_mean = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
x2_mean = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
ax = Axis(fig[1, 2], xticks = ([0.75, 2.25], ["Controllers", "Progressors"]), spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

scatter!(x1_mean, fill(median(dD_SIC), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x2_mean, fill(median(dD_VIR), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x1, dD_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(x2, dD_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
xlims!(ax, [0, 3]);

MannWhitneyUTest(dD_SIC, dD_VIR)

# Section 3c: r

x1 = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 12));
x2 = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 4));
x1_mean = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
x2_mean = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
ax = Axis(fig[1, 3], xticks = ([0.75, 2.25], ["Controllers", "Progressors"]), spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

scatter!(x1_mean, fill(median(log10.(r_SIC)), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x2_mean, fill(median(log10.(r_VIR)), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x1, log10.(r_SIC), strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(x2, log10.(r_VIR), strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
xlims!(ax, [0, 3]); ylims!(ax, [2, 4]);

MannWhitneyUTest(r_SIC, r_VIR)

# Section 3d: fD

x1 = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 12));
x2 = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 4));
x1_mean = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
x2_mean = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
ax = Axis(fig[2, 1], xticks = ([0.75, 2.25], ["Controllers", "Progressors"]), spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

scatter!(x1_mean, fill(median(fD_SIC), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x2_mean, fill(median(fD_VIR), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x1, fD_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(x2, fD_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
xlims!(ax, [0, 3]);

MannWhitneyUTest(fD_SIC, fD_VIR)

# Section 3e: ω

x1 = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 12));
x2 = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 4));
x1_mean = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
x2_mean = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
ax = Axis(fig[2, 2], xticks = ([0.75, 2.25], ["Controllers", "Progressors"]), spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

scatter!(x1_mean, fill(median(ω_SIC), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x2_mean, fill(median(ω_VIR), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x1, ω_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(x2, ω_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
xlims!(ax, [0, 3]);

MannWhitneyUTest(ω_SIC, ω_VIR)

fig
save("fig_S8.svg", fig, px_per_unit = 6)


