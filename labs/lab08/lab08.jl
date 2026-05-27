using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

p_cr = 42.0
N = 45.0
q = 1.0

tau1 = 28.0
tau2 = 22.0

p1 = 8.1
p2 = 10.5

M1_0 = 7.2
M2_0 = 9.1

u0 = [M1_0, M2_0]
tspan = (0.0, 30.0)

# Коэффициенты
a1 = p_cr / (tau1^2 * p1^2 * N * q)
a2 = p_cr / (tau2^2 * p2^2 * N * q)

b = p_cr / (tau1^2 * p1^2 * tau2^2 * p2^2 * N * q)

c1 = (p_cr - p1) / (tau1 * p1)
c2 = (p_cr - p2) / (tau2 * p2)

println("Коэффициенты:")
println("a1 = %.10f\n", a1)
println("a2 = %.10f\n", a2)
println("b  = %.10f\n", b)
println("c1 = %.10f\n", c1)
println("c2 = %.10f\n", c2)

# Случай 1

function case1!(du, u, p, t)
    M1 = u[1]
    M2 = u[2]

    du[1] = M1 - (b / c1) * M1 * M2 - (a1 / c1) * M1^2
    du[2] = (c2 / c1) * M2 - (b / c1) * M1 * M2 - (a2 / c1) * M2^2
end

# Случай 2

function case2!(du, u, p, t)
    M1 = u[1]
    M2 = u[2]

    du[1] = M1 - (b / c1 + 0.00048) * M1 * M2 - (a1 / c1) * M1^2
    du[2] = (c2 / c1) * M2 - (b / c1) * M1 * M2 - (a2 / c1) * M2^2
end

# Решение систем

prob1 = ODEProblem(case1!, u0, tspan)
sol1 = solve(prob1, Tsit5(), dtmax=0.1)

prob2 = ODEProblem(case2!, u0, tspan)
sol2 = solve(prob2, Tsit5(), dtmax=0.1)

theta = 0:0.01:30

y1 = reduce(hcat, sol1.(theta))
y2 = reduce(hcat, sol2.(theta))

# Стационарное состояние случая 1

D = a1 * a2 - b^2

M1_st = (c1 * a2 - b * c2) / D
M2_st = (a1 * c2 - b * c1) / D

println()
println("Стационарное состояние случая 1:")
println("M1* = %.4f\n", M1_st)
println("M2* = %.4f\n", M2_st)

# График случая 1

plot1 = plot(
    theta, y1[1, :],
    label="Фирма 1",
    linewidth=2,
    xlabel="θ = t / c1",
    ylabel="M",
    title="Случай 1"
)

plot!(
    theta, y1[2, :],
    label="Фирма 2",
    linewidth=2
)

# График случая 2

plot2 = plot(
    theta, y2[1, :],
    label="Фирма 1",
    linewidth=2,
    xlabel="θ = t / c1",
    ylabel="M",
    title="Случай 2"
)

plot!(
    theta, y2[2, :],
    label="Фирма 2",
    linewidth=2
)

# Вывод двух графиков

plot(plot1, plot2, layout=(2, 1), size=(900, 800))
savefig(plotsdir(script_name, "lab08.png"))
println("Графики сохранены: lab08_var58.png")
