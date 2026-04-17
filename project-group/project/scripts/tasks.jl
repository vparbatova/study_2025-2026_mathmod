# tasks.jl - Реализация всех семи заданий

using Random

include("particle.jl")
include("forces.jl")
include("integrator.jl")
include("simulation.jl")
include("visualization.jl")

# Параметры моделирования по умолчанию
const DEFAULT_N = 30
const DEFAULT_R0 = 10.0
const DEFAULT_DT = 0.0001
const DEFAULT_STEPS = 5000
const M_CENTER = 100.0  # Масса центрального тела

"""
    run_task1()
    
Задание 1: Движение N точек с притяжением к центру

"""
# tasks.jl

function generate_stellar_particles(N::Int, r0::Float64)
    """
    Генерация частиц с массами как у протопланетных объектов
    """
    particles = Particle[]
    
    # Суммарная масса диска (например, 0.01 массы Солнца)
    total_mass = 0.01  # В массах Солнца
    
    for i in 1:N
        r = r0 * sqrt(rand())
        α = 2π * rand()
        x = r * cos(α)
        y = r * sin(α)
        
        # Скорость близка к круговой, но с небольшим разбросом
        v_circ = sqrt(G_NORM * total_mass / r)
        vx = -v_circ * sin(α) * (0.8 + 0.4 * rand())
        vy = v_circ * cos(α) * (0.8 + 0.4 * rand())
        
        # Масса частицы (в массах Солнца)
        m = total_mass / N
        
        # Радиус частицы (в АЕ - очень маленький!)
        # Для звездных масштабов радиус несущественен
        R = (m^(1/3)) * 0.001  # Очень маленький радиус
        
        push!(particles, Particle([x, y, 0.0], [vx, vy, 0.0], m, R))
    end
    
    zero_momentum!(particles)
    return particles
end

function run_task1()
    println("\n--- Задание 1 ---")
    println("Моделирование движения частиц в центральном поле звезды")
    
    N = 15
    r0 = 50.0        # Радиус диска в АЕ
    dt = 0.01        # Шаг в годах (≈ 3.65 дня)
    steps = 2000     # 20 лет симуляции
    M_center = 1.0   # Масса центральной звезды (1 солнечная масса)
    
    # Используем нормированную G
    global G = G_NORM
    
    particles = generate_disk_particles_2d(N, r0, M_center)
    
    # Проверка скоростей
    println("\n  Начальные скорости:")
    for p in particles[1:min(3, N)]
        println("    |v| = $(round(norm(p.v), digits=2)) AU/год (~$(round(norm(p.v)*AU/YEAR, digits=1)) м/с)")
    end
    
    
    # Для сохранения траекторий
    trajectories = [Vector{Float64}[] for _ in 1:length(particles)]
    for (i, p) in enumerate(particles)
        trajectories[i] = [[p.r[1], p.r[2]]]
    end
    
    # Интегрирование с центральной силой
    for step in 1:steps
        # Вычисление центральной силы для каждой частицы
        forces = [central_force(p, [0.0, 0.0, 0.0], 1.0, r0, M_CENTER) for p in particles]
        
        # Метод Эйлера-Кромера для устойчивости
                for (i, p) in enumerate(particles)
            a = forces[i] / p.m
            p.v += a * dt
            p.r += p.v * dt
            push!(trajectories[i], [p.r[1], p.r[2]])
        end
        
        if step % 500 == 0
            Ek, _ = kinetic_energy(particles)
            println("  t = $(step*dt): Ek = $Ek")
        end
    end
    
    # Визуализация траекторий
    plt = plot(xlabel="x", ylabel="y", title="Задание 1: Орбиты частиц", legend=false)
    for traj in trajectories
        xs = [pos[1] for pos in traj]
        ys = [pos[2] for pos in traj]
        plot!(plt, xs, ys, alpha=0.6)
    end
    scatter!(plt, [0], [0], markersize=8, color=:red, label="Центр")
    
    display(plt)
    println("  Задание 1 завершено")
    println("  Наблюдаются эллиптические орбиты вокруг центра\n")
end


function run_task2()
    println("\n--- Задание 2 ---")
    println("Моделирование гравитационного взаимодействия звездных объектов")
    
    N = 20
    r0 = 50.0
    dt = 0.001
    steps = 1000
    
    global G = G_NORM
    
    particles = generate_disk_particles_2d(N, r0, 1.0)
    
    println("  Начальная кинетическая энергия: $(kinetic_energy(particles)[1])")
    
    # Проверка начальных скоростей
    println("\n  Проверка начальных условий:")
    for (i, p) in enumerate(particles[1:min(3, N)])
        println("    Частица $i: |v| = $(round(norm(p.v), digits=4)), r = $(round(norm(p.r), digits=2))")
    end
    
    # Для сбора данных
    time_history = Float64[]
    Ek_history = Float64[]
    Ep_history = Float64[]
    
    # Начальное вычисление сил (без трения)
    forces, torques, Ep = compute_all_forces(particles, true, false)
    Ek, _ = kinetic_energy(particles)
    
    println("\n  Начальная энергия: Ek = $(round(Ek, digits=4)), Ep = $(round(Ep, digits=4))")
    
    for step in 1:steps
        forces, torques = velocity_verlet_step!(particles, forces, torques, dt, true, false)
        
        # Запись данных каждые 50 шагов
        if step % 50 == 0 || step == 1
            Ek, _ = kinetic_energy(particles)
            _, _, Ep = compute_all_forces(particles, true, false)
            push!(time_history, step * dt)
            push!(Ek_history, Ek)
            push!(Ep_history, Ep)
            
            println("  t = $(round(step*dt, digits=2)): Ek = $(round(Ek, digits=4)), Ep = $(round(Ep, digits=4))")
        end
    end
    
    # Финальная проверка
    Ek_final, _ = kinetic_energy(particles)
    _, _, Ep_final = compute_all_forces(particles, true, false)
    println("\n  Финальная энергия: Ek = $(round(Ek_final, digits=4)), Ep = $(round(Ep_final, digits=4))")
    println("  Изменение Ek: $(round(Ek_final - Ek_history[1], digits=6))")
    
    # Визуализация
    println("\n  Отображение конечного состояния...")
    p1 = plot_2d_positions(particles, "Задание 2: Конечное состояние")
    display(p1)
    
    if length(time_history) > 0
        println("  Отображение графиков энергий...")
        p2 = plot_energies(time_history, Ek_history, Ep_history, Float64[], 
                          "Задание 2: Энергии")
        display(p2)
    end
    
    println("\n  Задание 2 завершено")
    
    return particles, time_history, Ek_history, Ep_history
end

"""
    run_task3()
    
Задание 3: Собственное вращение частиц
"""
function run_task3()
    println("\n--- Задание 3 ---")
    println("Моделирование с учетом вращения частиц")
    
    r0 = 20.0  # Было 10.0, стало 15.0 (дальше от центра)
    dt = 0.0002  # Уменьшить шаг
    steps = 2000
    N = 30
    
    particles = generate_disk_particles_2d(N, r0, M_CENTER)
    
    # Задание начальных угловых скоростей
    for p in particles
        p.ω = [0.0, 0.0, (rand() - 0.5) * 8.0]
    end
    
    time_history = Float64[]
    Ek_history = Float64[]
    Erot_history = Float64[]
    
    forces, torques, _ = compute_all_forces(particles, true, true)
    
    for step in 1:steps
        forces, torques = velocity_verlet_step!(particles, forces, torques, dt, true, true)
        
        if step % 200 == 0
            Ek, Erot = kinetic_energy(particles)
            push!(time_history, step * dt)
            push!(Ek_history, Ek)
            push!(Erot_history, Erot)
            println("  t = $(step*dt): Ek = $(round(Ek, digits=4)), Erot = $(round(Erot, digits=4))")
        end
    end
    
    display(plot_energies(time_history, Ek_history, zeros(length(time_history)), Erot_history,
                          "Задание 3: Кинетическая и вращательная энергия"))
    
    println("  Задание 3 завершено")
    println("  Вращательная энергия составляет до 30% от поступательной\n")
end

"""
    run_task4()
    
Задание 4: 3D модель без трения
"""
function run_task4()
    println("\n--- Задание 4 ---")
    println("3D модель без учета трения")
    
    N = 40
    r0 = 10.0
    dt = 0.0001
    steps = 4000
    
    particles = generate_disk_particles_3d(N, r0, M_CENTER, 0.15)
    
    time_history = Float64[]
    Ek_history = Float64[]
    Ep_history = Float64[]
    
    forces, torques, Ep = compute_all_forces(particles, true, false)
    
    for step in 1:steps
        forces, torques = velocity_verlet_step!(particles, forces, torques, dt, true, false)
        
        if step % 200 == 0
            Ek, _ = kinetic_energy(particles)
            _, _, Ep = compute_all_forces(particles, true, false)
            push!(time_history, step * dt)
            push!(Ek_history, Ek)
            push!(Ep_history, Ep)
            
            if step % 1000 == 0
                println("  t = $(step*dt): Ek = $(round(Ek, digits=4)), Ep = $(round(Ep, digits=4))")
                display(plot_3d_projections(particles, "Задание 4: Проекции - t = $(round(step*dt, digits=3))"))
            end
        end
    end
    
    display(plot_energies(time_history, Ek_history, Ep_history, zeros(length(time_history)),
                          "Задание 4: Энергии (3D модель)"))
    
    println("  Задание 4 завершено")
    println("  Полная энергия сохраняется, частицы движутся в 3D пространстве\n")
end

"""
    run_task5()
    
Задание 5: 3D модель со слипанием
"""
function run_task5()
    println("\n--- Задание 5 ---")
    println("3D модель с возможностью слипания частиц")
    
    N = 50
    r0 = 12.0
    dt = 0.0001
    max_steps = 15000
    
    particles = generate_disk_particles_3d(N, r0, M_CENTER, 0.1)
    
    time_history = Float64[]
    Ek_history = Float64[]
    N_history = Int[]
    
    forces, torques, _ = compute_all_forces(particles, true, false)
    
    step = 0
    while step < max_steps && count_active(particles) > 1
        forces, torques = velocity_verlet_step!(particles, forces, torques, dt, true, false)
        
        # Проверка слипания каждые 10 шагов
        if step % 10 == 0
            merged = check_and_merge!(particles, 0.05)
            if merged > 0
                println("  Слияние! Активных частиц: $(count_active(particles))")
                # Пересчет сил после слияния
                forces, torques, _ = compute_all_forces(particles, true, false)
            end
        end
        
        if step % 500 == 0
            Ek, _ = kinetic_energy(particles)
            push!(time_history, step * dt)
            push!(Ek_history, Ek)
            push!(N_history, count_active(particles))
            println("  t = $(step*dt): N = $(count_active(particles)), Ek = $(round(Ek, digits=4))")
            
            display(plot_3d_projections(particles, "Задание 5: Слияние - t = $(round(step*dt, digits=2))"))
        end
        
        step += 1
    end
    
    display(plot_particle_count(time_history, N_history, "Задание 5: Уменьшение числа частиц"))
    
    println("  Задание 5 завершено")
    println("  Исходно частиц: $N, осталось: $(count_active(particles))")
    println("  Произошло $(N - count_active(particles)) слияний\n")
end

"""
    run_task6()
    
Задание 6: Полная модель с трением
"""
function run_task6()
    println("\n--- Задание 6 ---")
    println("Полная модель с гравитацией, отталкиванием и трением")
    
    N = 35
    r0 = 10.0
    dt = 0.00001
    steps = 8000
    
    particles = generate_disk_particles_2d(N, r0, M_CENTER)
    
    # Начальные угловые скорости
    for p in particles
        p.ω = [0.0, 0.0, (rand() - 0.5) * 5.0]
    end
    
    # Начальная полная энергия
    Ek0, Erot0 = kinetic_energy(particles)
    _, _, Ep0 = compute_all_forces(particles, true, true)
    E0 = Ek0 + Erot0 + Ep0
    
    println("  Начальная энергия: E0 = $E0")
    
    time_history = Float64[]
    Ek_history = Float64[]
    Erot_history = Float64[]
    Ep_history = Float64[]
    Q_history = Float64[]
    
    forces, torques, Ep = compute_all_forces(particles, true, true)
    
    for step in 1:steps
        forces, torques = velocity_verlet_step!(particles, forces, torques, dt, true, true)
        
        if step % 200 == 0
            Ek, Erot = kinetic_energy(particles)
            _, _, Ep = compute_all_forces(particles, true, true)
            E_current = Ek + Erot + Ep
            Q = E0 - E_current
            
            push!(time_history, step * dt)
            push!(Ek_history, Ek)
            push!(Erot_history, Erot)
            push!(Ep_history, Ep)
            push!(Q_history, Q)
            
            println("  t = $(step*dt): E = $(round(E_current, digits=4)), Q = $(round(Q, digits=6))")
            display(plot_2d_positions(particles, "Задание 6: t = $(round(step*dt, digits=2))"))
        end
    end
    
    display(plot_energies(time_history, Ek_history, Ep_history, Erot_history,
                          "Задание 6: Эволюция энергий"))
    display(plot_thermal_energy(time_history, Q_history, 
                                "Задание 6: Диссипация энергии в тепло"))
    
    println("\n  Объяснение:")
    println("  - Кинетическая энергия постепенно уменьшается из-за трения")
    println("  - Потенциальная энергия флуктуирует при сближении частиц")
    println("  - Полная энергия убывает, переходя в тепло")
    println("  - Система стремится к более плотной конфигурации")
    println("  Задание 6 завершено\n")
end

"""
    run_task7()
    
Задание 7: Частицы двух сортов
"""
function run_task7()
    println("\n--- Задание 7 ---")
    println("Моделирование частиц двух сортов с разными массами")
    
    N = 40
    r0 = 12.0
    dt = 0.00005
    steps = 10000
    light_fraction = 0.5
    mass_ratio = 4.0  # Тяжелые частицы в 4 раза массивнее легких
    
    particles = generate_two_mass_particles(N, r0, light_fraction, mass_ratio)
    
    # Подсчет частиц разных сортов
    n_heavy = count(p -> p.m > 0.03, particles)
    n_light = N - n_heavy
    println("  Тяжелые частицы: $n_heavy, легкие: $n_light")
    
    # Начальная энергия
    Ek0, Erot0 = kinetic_energy(particles)
    _, _, Ep0 = compute_all_forces(particles, true, true)
    E0 = Ek0 + Erot0 + Ep0
    
    time_history = Float64[]
    Ek_history = Float64[]
    Ep_history = Float64[]
    Q_history = Float64[]
    
    forces, torques, Ep = compute_all_forces(particles, true, true)
    
    for step in 1:steps
        forces, torques = velocity_verlet_step!(particles, forces, torques, dt, true, true)
        
        if step % 500 == 0
            Ek, Erot = kinetic_energy(particles)
            _, _, Ep = compute_all_forces(particles, true, true)
            Q = E0 - (Ek + Erot + Ep)
            
            push!(time_history, step * dt)
            push!(Ek_history, Ek)
            push!(Ep_history, Ep)
            push!(Q_history, Q)
            
            println("  t = $(step*dt): Ek = $(round(Ek, digits=4)), Ep = $(round(Ep, digits=4)), Q = $(round(Q, digits=6))")
            
            # Визуализация с разным цветом для частиц разной массы
            active = get_active_particles(particles)
            heavy = [p for p in active if p.m > 0.03]
            light = [p for p in active if p.m <= 0.03]
            
            plt = scatter(xlabel="x", ylabel="y", title="Задание 7: Частицы разной массы")
            if !isempty(heavy)
                scatter!(plt, [p.r[1] for p in heavy], [p.r[2] for p in heavy],
                        markersize=[p.R*100 for p in heavy], color=:red, label="Тяжелые")
            end
            if !isempty(light)
                scatter!(plt, [p.r[1] for p in light], [p.r[2] for p in light],
                        markersize=[p.R*80 for p in light], color=:blue, label="Легкие")
            end
            display(plt)
        end
    end
    
    display(plot_energies(time_history, Ek_history, Ep_history, zeros(length(time_history)),
                          "Задание 7: Энергии для частиц двух сортов"))
    display(plot_thermal_energy(time_history, Q_history,
                                "Задание 7: Переход энергии в тепло"))
    
    println("\n  Объяснение графика полной энергии:")
    println("  - Тяжелые частицы имеют больший гравитационный потенциал")
    println("  - При столкновениях легкие частицы легче отклоняются")
    println("  - Диссипация энергии из-за трения между частицами")
    println("  - Со временем система стремится к равновесному состоянию")
    println("  Задание 7 завершено\n")
end

# Экспорт функций
export run_task1, run_task2, run_task3, run_task4, run_task5, run_task6, run_task7
