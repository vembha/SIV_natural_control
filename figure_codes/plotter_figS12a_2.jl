#===========================
PLOTTER FILE FOR 
FIGURE S12A_2
===========================#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
using HypothesisTests, LsqFit, EffectSizes
CairoMakie.activate!(type="svg")
size_cm = (5.0, 4.1);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

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
save("fig_S12a_2.svg", fig, px_per_unit = 6)

