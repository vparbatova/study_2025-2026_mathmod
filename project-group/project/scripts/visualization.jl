# visualization.jl - Функции для визуализации результатов

using Plots

# Настройка стиля Plots
Plots.default(linewidth=2, framestyle=:box, grid=true, gridalpha=0.3)

"""
    plot_2d_positions(particles, title)

Отображение положений частиц в 2D
"""
function plot_2d_positions(particles::Vector{Particle}, title::String="Положения частиц")
    active = get_active_particles(particles)
    
    xs = [p.r[1] for p in active]
    ys = [p.r[2] for p in active]
    sizes = [p.R * 100 for p in active]
    masses = [p.m for p in active]
    
    scatter(xs, ys, markersize=sizes, marker_z=masses, 
            xlabel="x", ylabel="y", title=title,
            color=:viridis, legend=false, size=(600, 600))
end

"""
    plot_3d_projections(particles, title)

Отображение проекций 3D системы на плоскости XY, XZ, YZ
"""
function plot_3d_projections(particles::Vector{Particle}, title::String="3D проекции")
    active = get_active_particles(particles)
    
    xs = [p.r[1] for p in active]
    ys = [p.r[2] for p in active]
    zs = [p.r[3] for p in active]
    sizes = [p.R * 80 for p in active]
    
    p1 = scatter(xs, ys, markersize=sizes, xlabel="x", ylabel="y", title="XY проекция", legend=false)
    p2 = scatter(xs, zs, markersize=sizes, xlabel="x", ylabel="z", title="XZ проекция", legend=false)
    p3 = scatter(ys, zs, markersize=sizes, xlabel="y", ylabel="z", title="YZ проекция", legend=false)
    
    plot(p1, p2, p3, layout=(1,3), size=(900, 400), title=title)
end

"""
    plot_energies(time_history, Ek_history, Ep_history, Erot_history, title)

Построение графиков энергий
"""
function plot_energies(time_history::Vector{Float64}, Ek_history::Vector{Float64}, 
                      Ep_history::Vector{Float64}, Erot_history::Vector{Float64}, 
                      title::String="Зависимость энергии от времени")
    
    if length(time_history) == 0
        println("Предупреждение: нет данных для построения графика энергий")
        return plot(title=title)
    end
    
    # Полная энергия
    Etot = [Ek_history[i] + Ep_history[i] + (i <= length(Erot_history) ? Erot_history[i] : 0.0) 
            for i in 1:length(time_history)]
    
    # Создание графика
    p = plot(xlabel="Время", ylabel="Энергия", title=title, legend=:topright)
    
    # Добавление линий
    plot!(p, time_history, Ek_history, label="Кинетическая", linewidth=2)
    plot!(p, time_history, Ep_history, label="Потенциальная", linewidth=2)
    
    if length(Erot_history) == length(time_history) && sum(abs.(Erot_history)) > 0
        plot!(p, time_history, Erot_history, label="Вращательная", linewidth=2, linestyle=:dash)
    end
    
    plot!(p, time_history, Etot, label="Полная", linewidth=3, linestyle=:dot, color=:black)
    
    return p
end

"""
    plot_thermal_energy(time_history, Q_history, title)

Построение графика тепловой энергии
"""
function plot_thermal_energy(time_history::Vector{Float64}, Q_history::Vector{Float64}, 
                             title::String="Энергия, перешедшая в тепло")
    
    if length(time_history) == 0
        println("Предупреждение: нет данных для построения графика тепловой энергии")
        return plot(title=title)
    end
    
    p = plot(time_history, Q_history, xlabel="Время", ylabel="Q (тепло)", 
             title=title, linewidth=2, color=:red, legend=false)
    
    # Добавление заливки под кривой
    plot!(p, time_history, Q_history, fillrange=0, fillalpha=0.3, color=:red)
    
    return p
end

"""
    plot_particle_count(time_history, N_history, title)

Построение графика изменения числа частиц
"""
function plot_particle_count(time_history::Vector{Float64}, N_history::Vector{Int}, 
                             title::String="Эволюция числа частиц")
    
    if length(time_history) == 0
        println("Предупреждение: нет данных для построения графика числа частиц")
        return plot(title=title)
    end
    
    p = plot(time_history, N_history, xlabel="Время", ylabel="Число частиц",
             title=title, linewidth=2, color=:blue, marker=:circle, 
             markersize=4, legend=false)
    
    return p
end

"""
    plot_orbits(trajectories, title)

Построение траекторий частиц (для задания 1)
"""
function plot_orbits(trajectories::Vector{Vector{Vector{Float64}}}, 
                     title::String="Траектории частиц")
    
    p = plot(xlabel="x", ylabel="y", title=title, legend=false, size=(800, 600))
    
    for (i, traj) in enumerate(trajectories)
        if length(traj) > 0
            xs = [pos[1] for pos in traj]
            ys = [pos[2] for pos in traj]
            plot!(p, xs, ys, alpha=0.5 + i/length(trajectories)*0.3, linewidth=1)
        end
    end
    
    # Отметка центра
    scatter!(p, [0], [0], markersize=10, color=:red, label="Центр")
    
    return p
end

"""
    plot_particles_colored_by_mass(particles, title)

Визуализация частиц с разным цветом для разных масс
"""
function plot_particles_colored_by_mass(particles::Vector{Particle}, 
                                        title::String="Частицы разной массы")
    active = get_active_particles(particles)
    
    if isempty(active)
        println("Нет активных частиц")
        return plot(title=title)
    end
    
    # Разделение на тяжелые и легкие
    avg_mass = sum(p.m for p in active) / length(active)
    heavy = [p for p in active if p.m > avg_mass]
    light = [p for p in active if p.m <= avg_mass]
    
    p = plot(xlabel="x", ylabel="y", title=title, size=(700, 600), legend=:topright)
    
    if !isempty(heavy)
        scatter!(p, [p.r[1] for p in heavy], [p.r[2] for p in heavy],
                markersize=[p.R*100 for p in heavy], color=:red, 
                label="Тяжелые частицы", alpha=0.8)
    end
    
    if !isempty(light)
        scatter!(p, [p.r[1] for p in light], [p.r[2] for p in light],
                markersize=[p.R*80 for p in light], color=:blue, 
                label="Легкие частицы", alpha=0.6)
    end
    
    return p
end

"""
    print_energy_info(Ek, Ep, Erot, step, dt)

Вывод информации об энергиях в консоль
"""
function print_energy_info(Ek::Float64, Ep::Float64, Erot::Float64, 
                           step::Int, dt::Float64, prefix::String="")
    t = step * dt
    Etot = Ek + Ep + Erot
    println("$prefix t = $(round(t, digits=4)): Ek = $(round(Ek, digits=6)), " *
            "Ep = $(round(Ep, digits=6)), Erot = $(round(Erot, digits=6)), " *
            "E = $(round(Etot, digits=6))")
end
# visualization.jl - добавить масштабирование

function plot_2d_positions(particles::Vector{Particle}, title::String="Положения частиц")
    active = get_active_particles(particles)
    
    xs = [p.r[1] for p in active]
    ys = [p.r[2] for p in active]
    masses = [p.m / maximum([p.m for p in active]) * 100 for p in active]
    
    scatter(xs, ys, markersize=sqrt.(masses) * 10, 
            xlabel="x (AU)", ylabel="y (AU)", title=title,
            legend=false, size=(700, 700))
end
