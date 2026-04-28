using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# начальные условия
x0 = 6.0 
y0 = 23.0
u0 = [x0, y0]
tspan = (0.0, 200.0)

# коэффициенты
a = 0.38
b = 0.043
c = 0.36
d = 0.052

function predator_prey!(du, u, p, t)
   x = u[1] 
   y = u[2]
   
   du[1] = -a * x + b * x * y
   du[2] = c* y - d * x * y
end

prob = ODEProblem(predator_prey!, u0, tspan)
sol = solve(prob, Tsit5(), dtmax = 0.08)

x_val = sol[1, :]
y_val = sol[2, :]
time_val = sol.t

# 1 график
p1 = plot(time_val, x_val, label="Хищники (х)", color=:red)
plot!(p1, time_val, y_val, label="Жертвы (у)", color=:blue)
title!(p1, "Изменение численности хищников и жертв во времени")
xlabel!(p1, "Время t")
ylabel!(p1, "Численность популяции")

# 2 график
p2 = plot(y_val, x_val, label="Фазовая траетория", color=:green)
scatter!(p2, [y0], [x0], label="Начальная точка", color=:black)

# стационарное состояние
x_eq = c / d
y_eq = a / b
scatter!(p2, [y_eq], [x_eq], label = "Начальная точка", color = :black)

title!(p2, "Фазовый портрет: численность хищников от численности жертв")
xlabel!(p2, "Жертвы (у)")
ylabel!(p2, "Хищники (х)")

solutions = plot(p1, p2, layout=(2,1), size=(900, 800))

savefig(plotsdir(script_name, "lab05.png"))
println("Графики сохранены: lab05_var58.png")

