#===========================
PLOTTER FILE FOR 
FIGURE S17A_1
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
df = DataFrame(XLSX.readtable("fig_S17a_1_data.xlsx","Sheet1",infer_eltypes=true));

#=
Section 2: Plot accordingly
=#

# Get the required data
idx = 3; # 1 for S₁₄, 2 for S₄₂ and 3 for S₉₀

if idx == 1
    x = df.tAUC_14;
    x_lim  = [-0.05, 2.0];
    y_lim  = [-3.2, 6.0];
    x_tick = [0.0, 0.50, 1.0, 1.5, 2.0];
    y_tick = [-2, 0, 2, 4, 6];
elseif idx == 2
    x = df.tAUC_42;
    x_lim  = [-0.05, 2.0];
    y_lim  = [-3.2, 6.0];
    x_tick = [0.0, 0.50, 1.0, 1.5, 2.0];
    y_tick = [-2, 0, 2, 4, 6];
elseif idx == 3
    x = df.tAUC_90;
    x_lim  = [-0.05, 2.0];
    y_lim  = [-3.2, 6.0];
    x_tick = [0.0, 0.50, 1.0, 1.5, 2.0];
    y_tick = [-2, 0, 2, 4, 6];
end
y = df.spVL;
x_SIC = x[1:12]; x_VIR = x[13:16];
y_SIC = y[1:12]; y_VIR = y[13:16];

# Defining the linear model abd obtain the regression line
@. lin(x, p) = p[1] + p[2] * x;
p0 = [1.0, -1.0];
fit_S10c = curve_fit(lin, x, y, p0); # Equation is: VLf = slp*AUC + int
int = coef(fit_S10c)[1]; slp = coef(fit_S10c)[2];
@. fit_model(x) = int + slp*x;

# Setting the panel
fig = Figure(resolution = size_pt, fontsize=12);
ax = Axis(fig[1, 1], spinewidth = 1, xticks = x_tick, yticks = y_tick,  xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

# Plotting commands
scatter!(ax, x_SIC, y_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(ax, x_VIR, y_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
lines!(ax, sort(x), fit_model.(sort(x)), color = :black, linewidth = 1);

xlims!(ax, x_lim);
ylims!(ax, y_lim);

fig
if idx == 1
    save("fig_S17a_1_1.svg", fig, px_per_unit = 6)
elseif idx == 2
    save("fig_S17a_1_2.svg", fig, px_per_unit = 6)
elseif idx == 3
    save("fig_S17a_1_3.svg", fig, px_per_unit = 6)
end

