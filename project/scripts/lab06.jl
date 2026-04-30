using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# начальные условия
N = 17854.0
I0 = 199.0
R0 = 35.0
S0 = N - I0 - R0

alpha = 0.01
beta = 0.02

u0 = [S0, I0, R0]
tspan = (0.0, 200.0)

# система оду
function epidemic!(du, u, p, t)
    S = u[1]
    I = u[2]
    R = u[3]
    I_star = p[1]

    if I > I_star
       du[1] = -alpha * S
       du[2] = alpha * S - beta * I 
       du[3] = beta * I 
    else
       du[1] = 0.0 
       du[2] = -beta * I
       du[3] = beta * I
    end
end

# случай 1: I(0) <= I*
I_star_1 = 300.0
prob1 = ODEProblem(epidemic!, u0, tspan, [I_star_1])
sol1 = solve(prob1, Tsit5(), dtmax=0.1)

# случай 2: I(0) > I*
I_star_2 = 100.0
prob2 = ODEProblem(epidemic!, u0, tspan, [I_star_2])
sol2 = solve(prob2, Tsit5(), dtmax=0.1)

p1 = plot(sol1, label = ["S(t)" "I(t)" "R(t)"], lw=2, color = [:blue :red :green], 
    title = "Случай 1: I(0) <= I* (изоляция)", xlabel = "Время", ylabel = "Численность особей")
    
p2 = plot(sol2, label = ["S(t)" "I(t)" "R(t)"], lw=2, color = [:blue :red :green], 
    title = "Случай 2: I(0) > I* (эпидемия)", xlabel = "Время", ylabel = "Численность особей")
    
fig = plot(p1, p2, layout=(2,1), size=(900, 800))
savefig(plotsdir(script_name, "lab06.png"))
println("Графики сохранены: lab06_var58.png")
