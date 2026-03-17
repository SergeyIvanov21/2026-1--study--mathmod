using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

x0 = 100_000.0
y0 = 10_000.0
 
# СЛУЧАЙ 1
function model1!(du, u, p, t)
    x, y = u[1], u[2]
    du[1] = -0.12*x - 0.9*y  + abs(sin(t))
    du[2] = -0.3*x  - 0.1*y  + abs(cos(t))
end
 
u0      = [x0, y0]
tspan1  = (0.0, 3.0)
 
prob1 = ODEProblem(model1!, u0, tspan1)
sol1  = solve(prob1, Tsit5(), dtmax = 0.005)
 
# Определяем победителя
x_final1 = sol1[1, end]
y_final1 = sol1[2, end]
if x_final1 > 0 && y_final1 <= 0
    winner1 = "Победа армии X"
elseif y_final1 > 0 && x_final1 <= 0
    winner1 = "Победа армии Y"
else
    winner1 = x_final1 > y_final1 ? "Преимущество у армии X" : "Преимущество у армии Y"
end
 
println("Случай 1 — Регулярные войска:")
println("  X(конец) ≈ $(round(x_final1, digits=1))")
println("  Y(конец) ≈ $(round(y_final1, digits=1))")
println("  $winner1\n")
 
p1 = plot(
    sol1.t, sol1[1, :],
    label  = "Армия X",
    color  = :steelblue,
    lw     = 2.5,
    title  = "Случай 1: регулярные войска",
    xlabel = "Время",
    ylabel = "Численность",
    legend = :topright
)
plot!(p1, sol1.t, sol1[2, :], label = "Армия Y", color = :firebrick, lw = 2.5)
hline!(p1, [0], color = :black, linestyle = :dash, label = "")
 

# СЛУЧАЙ 2
function model2!(du, u, p, t)
    x, y = u[1], u[2]
    du[1] = -0.25*x        - 0.96*y      + sin(2*t)  + 1.0
    du[2] = -0.25*x*y      - 0.3*y       + cos(20*t) + 1.0
end
 
tspan2 = (0.0, 0.8)   
 
prob2 = ODEProblem(model2!, u0, tspan2)
sol2  = solve(prob2, Tsit5(), dtmax = 0.001)
 
x_final2 = sol2[1, end]
y_final2 = sol2[2, end]
if x_final2 > 0 && y_final2 <= 0
    winner2 = "Победа армии X"
elseif y_final2 > 0 && x_final2 <= 0
    winner2 = "Победа армии Y"
else
    winner2 = x_final2 > y_final2 ? "Преимущество у армии X" : "Преимущество у армии Y"
end
 
println("Случай 2 — Регулярные войска + партизаны:")
println("  X(конец) ≈ $(round(x_final2, digits=1))")
println("  Y(конец) ≈ $(round(y_final2, digits=1))")
println("  $winner2\n")
 
p2 = plot(
    sol2.t, sol2[1, :],
    label  = "Армия X (рег.)",
    color  = :steelblue,
    lw     = 2.5,
    title  = "Случай 2: регулярные + партизаны",
    xlabel = "Время",
    ylabel = "Численность",
    legend = :topright
)
plot!(p2, sol2.t, sol2[2, :], label = "Армия Y (парт.)", color = :firebrick, lw = 2.5)
hline!(p2, [0], color = :black, linestyle = :dash, label = "")
 
# Сохранение
combined = plot(p1, p2, layout = (2, 1), size = (850, 900))
savefig(plotsdir(script_name, "lab03.png"))
println("Графики сохранены: lab03_var58.png")
 



