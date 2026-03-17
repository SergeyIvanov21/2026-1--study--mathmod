using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# вариант 58
s = 20.2          # Начальное расстояние 
n = 5.1           # Отношение скоростей
fi = 3*pi/4       # Направление движения лодки (можно выбрать любое)

# Начальные расстояния для двух случаев
r0_case1 = s / (n + 1)  # Случай 1: катер между полюсом и лодкой
r0_case2 = s / (n - 1)  # Случай 2: лодка между полюсом и катером

println("Начальные расстояния:")
println("Случай 1: r0 = ", round(r0_case1, digits=3))
println("Случай 2: r0 = ", round(r0_case2, digits=3))

# Константа из уравнения dr/dθ = r/√(n² - 1)
const_factor = 1 / sqrt(n^2 - 1)
println("Коэффициент: 1/√(n²-1) = ", round(const_factor, digits=3))

# Дифференциальное уравнение движения катера
# θ - угол, r - радиус
function patrol_boat_ode(r, p, theta)
    return r / sqrt(n^2 - 1)
end

# Решение для случая 1
theta_span1 = (0.0, 2*pi)      # Интервал углов
prob1 = ODEProblem(patrol_boat_ode, r0_case1, theta_span1)
sol1 = solve(prob1, abstol=1e-8, reltol=1e-8)

# Решение для случая 2 (начинаем с противоположной стороны)
theta_span2 = (-pi, pi)         # Интервал углов
prob2 = ODEProblem(patrol_boat_ode, r0_case2, theta_span2)
sol2 = solve(prob2, abstol=1e-8, reltol=1e-8)

# Движение лодки (прямолинейно под углом fi)
boat_x(t) = t * tan(fi)         # Параметрическое уравнение лодки
# В полярных координатах лодка движется по лучу r = v·t, θ = fi
# Для графика нам нужна прямая линия под углом fi

# Построение графиков
# Случай 1
p1 = plot(sol1, proj=:polar, 
          label="Катер (случай 1)", 
          title="Случай 1: r0 = s/(n+1) = $(round(r0_case1, digits=2))",
          c=:red, lw=2, 
          xlabel="Угол θ", ylabel="Радиус r")

# Добавляем траекторию лодки (луч под углом fi)
# Для отображения лодки в полярных координатах построим луч
theta_boat = range(0, 2*pi, length=100)
r_boat = range(0, maximum(sol1.u)*1.2, length=100)

# Создаем массив точек для луча
boat_r = [r for r in r_boat]
boat_theta = [fi for _ in r_boat]
plot!(p1, boat_theta, boat_r, proj=:polar,
      label="Лодка (θ = $(round(fi, digits=2)))", 
      c=:green, lw=2, linestyle=:dash)

# Случай 2
p2 = plot(sol2, proj=:polar, 
          label="Катер (случай 2)", 
          title="Случай 2: r0 = s/(n-1) = $(round(r0_case2, digits=2))",
          c=:blue, lw=2,
          xlabel="Угол θ", ylabel="Радиус r")

# Добавляем траекторию лодки
boat_r2 = [r for r in range(0, maximum(sol2.u)*1.2, length=100)]
boat_theta2 = [fi for _ in boat_r2]
plot!(p2, boat_theta2, boat_r2, proj=:polar,
      label="Лодка (θ = $(round(fi, digits=2)))", 
      c=:green, lw=2, linestyle=:dash)

plot(p1, p2, layout=(2, 1), size=(700, 900))
savefig(plotsdir(script_name, "lab02.png"))

# Находим точку пересечения
# Для случая 1: катер догонит лодку, когда θ = fi и r совпадут
# В этот момент r_boat = v·t, а r_catcher из решения уравнения

println("Катер догонит лодку, когда его траектория пересечет луч лодки (θ = fi)")
println("Для случая 1: точка пересечения при θ = $fi")
println("Для случая 2: точка пересечения также при θ = $fi")
println("Радиус в точке встречи можно найти из уравнения r(θ) = r0·exp(θ/√(n²-1))")
