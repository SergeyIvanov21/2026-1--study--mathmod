using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

N = 761.0
n0 = 6.0
u0 = [n0]

tspan_1 = (0.0, 12.0)
tspan_2 = (0.0, 0.02)
tspan_3 = (0.0, 0.02)

#случай 1
function ad_campaign1!(du, u, p, t)
    n = u[1]
    du[1] = (0.82 + 0.00003 * n) * (N - n)
end

#случай 2
function ad_campaign2!(du, u, p, t)
    n = u[1]
    du[1] = (0.00003 + 0.82 * n) * (N - n)
end

#случай 3
function ad_campaign3!(du, u, p, t)
    n = u[1]
    du[1] = (0.2 * sin(t) + 0.8 * cos(t) * n) * (N - n)
end

prob1 = ODEProblem(ad_campaign1!, u0, tspan_1)
sol1 = solve(prob1, Tsit5(), dtmax=0.1)

prob2 = ODEProblem(ad_campaign2!, u0, tspan_2)
sol2 = solve(prob2, Tsit5(), dtmax=0.1)

prob3 = ODEProblem(ad_campaign3!, u0, tspan_3)
sol3 = solve(prob3, Tsit5(), dtmax=0.1)

#поиск момента максимальной скорости для случая 2
function speed_case2(n)
   return(0.00003 + 0.82 * n) * (N - n)
end

t_dense2 = range(tspan_2[1], tspan_2[2], length=10000)
speeds2 = [speed_case2(sol2(t)[1]) for t in t_dense2]

max_speed_value, max_speed_index = findmax(speeds2)
max_speed_time = t_dense2[max_speed_index]
max_speed_n = sol2(max_speed_time)[1]

println("Случай 2")
println("Максимальная скорость распространения рекламы: ")
println("t = ", round(max_speed_time, digits=6))
println("n(t) = ", round(max_speed_n, digits=3), " человек за единицу времени")

#построение графика
t1 = range(tspan_1[1], tspan_1[2], length=1000)
n1 = [sol1(t)[1] for t in t1]

t2 = range(tspan_2[1], tspan_2[2], length=1000)
n2 = [sol2(t)[1] for t in t2]

t3 = range(tspan_3[1], tspan_3[2], length=1000)
n3 = [sol3(t)[1] for t in t3]

p1 = plot(t1, n1, title="Случай 1", label="n(t)", xlabel="t", ylabel="n(t)", 
    lw=2, color=:blue, legend=:bottomright)
    
p2 = plot(t2, n2, title="Случай 2", label="n(t)", xlabel="t", ylabel="n(t)", 
    lw=2, color=:red, legend=:bottomright)
    
scatter!(p2, [max_speed_time], [max_speed_n], label="Макс скорость", 
    color=:black, markersize=5)
   
p3 = plot(t3, n3, title="Случай 3", label="n(t)", xlabel="t", ylabel="n(t)", 
    lw=2, color=:green, legend=:bottomright)
    
fig = plot(p1, p2, p3, layout=(3,1), size=(900, 1100))
savefig(plotsdir(script_name, "lab07.png"))
println("Графики сохранены: lab07_var58.png")
