#==========================
CODE FOR SENSITIVITY 
ANALYSIS OF THE MODEL #1 
USING PARAMETER ESTIMATES 
FROM THE BEST-FIT MODEL
==========================#

#=
Section 1: Loading packages
=#

using DiffEqSensitivity, GlobalSensitivity, DifferentialEquations
using CairoMakie, Statistics
CairoMakie.activate!(type="svg")
size_cm = (13.5, 4.5);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

#=
Section 2: Defining the required parameter values, 
the model system, and the objective function
=#

λ = 352.7; β_prime = 10^(-2.84); fD = 0.92; μI = 0.1;
μD = 0.069; r = 427.65; αE = 0.63; θE = 0.1;
μE = 1.0; ω = 10^(-2.50); T0 = 10^(4.21);
μT = λ/T0; u0 = [T0, (10^(-2.76))/r, 0.0, 0.0, 0.0]; λE⁺ = 10^(0.15);
tspan = (0.0, 365.0); param = [λ, β_prime, fD, μI, μD, r, λE⁺, αE, θE, μE, ω, μT];

function model_1!(du, u, p, t)

    T, I, D, E⁺, k⁺ = u; # Variables
    λ, β_prime, fD, μI, μD, r, λE⁺, αE, θE, μE, ω, μT = p; # Parameters

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - k⁺*E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + αE*E⁺*I/(θE + I) - μE*E⁺
    du[5] = dk⁺ = ω*(1 - k⁺)
end

prob = ODEProblem(model_1!, u0, tspan, param);
sol = solve(prob);

function spVL!(p)
    u0 = [T0, (10^(-2.76))/p[6], 0.0, 0.0]
    u0[1] = p[1]/p[12]
    prob1 = remake(prob; p=p)
    sol1  = solve(prob1, isoutofdomain = (u,p,t)->any(x->x<0,u))
    return log10(2*p[6]*max(sol1[2,end], 1e-8))
end

#=
Section 3: Perform Sobol sensitivity analysis
=#

# Sensitivity of spVL
res = gsa(spVL!, Sobol(), [[0.1*λ, 10*λ], [0.1*β_prime, 10*β_prime], [0.1*fD, 1.0], [0.1*μI, 10*μI], [0.1*μD, 10*μD], [0.1*r, 10*r], [0.1*λE⁺, 10*λE⁺], [0.1*αE, 10*αE], [0.1*θE, 10*θE], [0.1*μE, 10*μE], [0.1ω, 10*ω], [0.1*μT, 10*μT]]; samples = 10_000)

#=
Section 4: Plot
=#

fig = Figure(resolution = size_pt, fontsize=12);
ax = Axis(fig[1, 1], xticks = (1:12, ["λ", "β'", "fD", "μI", "μD", "r", "λE⁺", "αE", "θE", "μE", "ω", "μT"]), spinewidth = 1, xtickwidth = 1, ytickwidth = 1, yticks = [0.0, 0.25, 0.5, 0.75, 1.0], titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

barplot!([1:1:12; 1:1:12], [res.ST; res.S1], dodge = Int.([ones(12); 2*ones(12)]), color = Int.([ones(12); 2*ones(12)]), strokecolor = :black, strokewidth = 1);
ylims!(ax, [0.0, 1.0]);

fig
save("fig_S7.svg", fig, px_per_unit = 6)

