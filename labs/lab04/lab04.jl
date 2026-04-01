using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# начальные условия
x0_val = 0.1
y0_val = 1.1
u0 = [x0_val, y0_val]
tspan = (0.0, 63.0)
dt = 0.05

# случай 1. x'' + 3.7x = 0
function harmonical1!(du, u, p, t)
    g = 0.0 
    w2 = 3.7
    
    du[1] = u[2]
    du[2] = -g * u[2] - w2 * u[1]
end

prob1 = ODEProblem(harmonical1!, u0, tspan)
sol1 = solve(prob1, Tsit5(), saveat=dt)

# графики
plt1_sol = plot(sol1.t, [u[1] for u in sol1.u],
    title = "Случай 1: решение x(t)",
    xlabel = "t", ylabel = "x(t)",
    label = "x(t)", lw = 1.5, color = :blue)
    
plt1_phase = plot([u[1] for u in sol1.u], [u[2] for u in sol1.u],
    title = "Случай 1: фазовый портрет",
    xlabel = "x", ylabel = "x'",
    label = false, lw = 1.5, color = :blue, 
    aspect_ratio = :equal)
    
# случай 2. x'' + 3x' + 10x = 0
function harmonical2!(du, u, p, t)
   g = 3.0 
   w2 = 10.0
   
   du[1] = u[2]
   du[2] = -g * u[2] - w2 * u[1]
end

prob2 = ODEProblem(harmonical2!, u0, tspan)
sol2 = solve(prob2, Tsit5(), saveat=dt)
   
plt2_sol = plot(sol2.t, [u[1] for u in sol2.u],
    title = "Случай 2: решение x(t)",
    xlabel = "t", ylabel = "x(t)",
    label = "x(t)", lw = 1.5, color = :green)
    
plt2_phase = plot([u[1] for u in sol2.u], [u[2] for u in sol2.u],
    title = "Случай 2: фазовый портрет",
    xlabel = "x", ylabel = "x'",
    label = false, lw = 1.5, color = :green)
    
# случай. x'' + 3x' + 11x = 0.9sin(0.9t)
function harmonical3!(du, u, p, t)
   g = 3.0 
   w2 = 11.0
   f_t = 0.9 * sin(0.9 * t)
   
   du[1] = u[2]
   du[2] = -g * u[2] - w2 * u[1] + f_t
end

prob3 = ODEProblem(harmonical3!, u0, tspan)
sol3 = solve(prob3, Tsit5(), saveat=dt)

plt3_sol = plot(sol3.t, [u[1] for u in sol3.u],
    title = "Случай 3: решение x(t)",
    xlabel = "t", ylabel = "x(t)",
    label = "x(t)", lw = 1.5, color = :red)
    
plt3_phase = plot([u[1] for u in sol3.u], [u[2] for u in sol3.u],
    title = "Случай 3: фазовый портрет",
    xlabel = "x", ylabel = "x'",
    label = false, lw = 1.5, color = :red)
    
    
solutions = plot(plt1_sol, plt2_sol, plt3_sol, 
    layout = (1, 3), size = (1400, 420))
    
savefig(plotsdir(script_name, "lab04.png"))
println("Графики сохранены: lab04_var58.png")
 

