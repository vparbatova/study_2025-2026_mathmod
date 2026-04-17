# integrator.jl - Метод Верле для интегрирования уравнений движения

using LinearAlgebra

"""
    velocity_verlet_step!(particles, forces, torques, dt, use_gravity, use_friction)
    
Один шаг метода Верле по скорости
"""
function velocity_verlet_step!(particles::Vector{Particle}, forces::Vector{Vector{Float64}}, 
                               torques::Vector{Vector{Float64}}, dt::Float64, 
                               use_gravity::Bool, use_friction::Bool)
    
    # Сохраняем силы для второго полушага
    forces_old = [copy(f) for f in forces]
    torques_old = [copy(τ) for τ in torques]
    
    # Первый полушаг: обновление скоростей (половина шага)
    for i in eachindex(particles)
        if !particles[i].active
            continue
        end
        
        a = forces[i] / particles[i].m
        particles[i].v += a * (dt / 2.0)
        
        # Обновление угловой скорости
        ε = torques[i] / particles[i].I
        particles[i].ω += ε * (dt / 2.0)
    end
    
    # Обновление положений (полный шаг)
    for i in eachindex(particles)
        if particles[i].active
            particles[i].r += particles[i].v * dt
        end
    end
    
    # Пересчет сил на новых положениях
    forces_new, torques_new, _ = compute_all_forces(particles, use_gravity, use_friction)
    
    # Второй полушаг: обновление скоростей
    for i in eachindex(particles)
        if !particles[i].active
            continue
        end
        
        a_avg = (forces[i] + forces_new[i]) / (2.0 * particles[i].m)
        particles[i].v += a_avg * dt
        
        ε_avg = (torques[i] + torques_new[i]) / (2.0 * particles[i].I)
        particles[i].ω += ε_avg * dt
    end
    
    return (forces_new, torques_new)
end

"""
    velocity_verlet!(particles, dt, steps, use_gravity, use_friction, callback)
    
Выполнение серии шагов интегрирования
"""
function velocity_verlet!(particles::Vector{Particle}, dt::Float64, steps::Int, 
                          use_gravity::Bool, use_friction::Bool, callback=nothing)
    
    # Начальное вычисление сил
    forces, torques, _ = compute_all_forces(particles, use_gravity, use_friction)
    
    history = []
    
    for step in 1:steps
        forces, torques = velocity_verlet_step!(particles, forces, torques, dt, use_gravity, use_friction)
        
        # Вызов callback функции для мониторинга
        if callback !== nothing
            if step % 100 == 0 || step == steps
                callback(particles, step, dt)
            end
        end
    end
    
    return history
end
