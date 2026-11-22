using Plots

struct Caixa
    largura::Float64
    altura::Float64
end

mutable struct Particula2D
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
    m::Float64
    r::Float64
end

function mover!(p::Particula2D, dt::Real)
    p.x += p.vx * dt
    p.y += p.vy * dt
end

function colidir_paredes!(p::Particula2D, caixa::Caixa)
    if p.x - p.r <= 0.0 && p.vx < 0
        p.vx = -p.vx
        p.x = p.r
    elseif p.x + p.r >= caixa.largura && p.vx > 0
        p.vx = -p.vx
        p.x = caixa.largura - p.r
    end

    if p.y - p.r <= 0.0 && p.vy < 0
        p.vy = -p.vy
        p.y = p.r
    elseif p.y + p.r >= caixa.altura && p.vy > 0
        p.vy = -p.vy
        p.y = caixa.altura - p.r
    end
end

function colidir!(p1::Particula2D, p2::Particula2D)
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    dist2 = dx^2 + dy^2
    soma_r = p1.r + p2.r

    if dist2 <= soma_r^2 && dist2 > 0
        dvx = p2.vx - p1.vx
        dvy = p2.vy - p1.vy

        if dx * dvx + dy * dvy < 0
            dist = sqrt(dist2)
            nx = dx / dist
            ny = dy / dist

            vn = dvx * nx + dvy * ny
            m1, m2 = p1.m, p2.m
            impulso = (2 * vn) / (m1 + m2)

            p1.vx += impulso * m2 * nx
            p1.vy += impulso * m2 * ny
            p2.vx -= impulso * m1 * nx
            p2.vy -= impulso * m1 * ny
        end
    end
end

function verificar_colisoes!(particulas::Vector{Particula2D})
    n = length(particulas)
    for i in 1:n
        for j in (i + 1):n
            colidir!(particulas[i], particulas[j])
        end
    end
end

energia_cinetica(ps) = sum(0.5 * p.m * (p.vx^2 + p.vy^2) for p in ps)

function criar_gas(n::Int, caixa::Caixa; r=0.3, m=1.0, v_max=3.0)
    particulas = Particula2D[]

    while length(particulas) < n
        # Posição aleatória respeitando a margem do raio nas paredes
        x = r + rand() * (caixa.largura - 2r)
        y = r + rand() * (caixa.altura - 2r)

        # Rejeita se encostar em alguma partícula já criada
        sobreposta = any((p.x - x)^2 + (p.y - y)^2 < (p.r + r)^2 for p in particulas)
        if sobreposta
            continue
        end

        # Velocidade aleatória entre -v_max e +v_max
        vx = (2 * rand() - 1) * v_max
        vy = (2 * rand() - 1) * v_max

        push!(particulas, Particula2D(x, y, vx, vy, m, r))
    end

    return particulas
end

# 1. Configuração do mini gás
caixa = Caixa(10.0, 10.0)
particulas = criar_gas(20, caixa, r=0.3, v_max=3.0)

# Deixa 19 azuis e pinta a primeira de vermelho para você conseguir acompanhá-la
cores = fill(:dodgerblue, length(particulas))
cores[1] = :crimson

dt = 0.02
passos = 400   # 8 segundos simulados

E_inicial = energia_cinetica(particulas)

# 2. Simulação
n = length(particulas)
tempos = Float64[]
pos_x = Matrix{Float64}(undef, passos + 1, n)
pos_y = Matrix{Float64}(undef, passos + 1, n)

for passo in 0:passos
    push!(tempos, passo * dt)
    for (i, p) in enumerate(particulas)
        pos_x[passo + 1, i] = p.x
        pos_y[passo + 1, i] = p.y
    end

    for p in particulas
        mover!(p, dt)
        colidir_paredes!(p, caixa)
    end
    verificar_colisoes!(particulas)
end

println("Energia -> antes: $(round(E_inicial, digits=4)) | depois: $(round(energia_cinetica(particulas), digits=4))")

# 3. Gráfico do rastro das trajetórias
rastro = plot(
    pos_x[:, 1], pos_y[:, 1],
    xlims = (0, caixa.largura),
    ylims = (0, caixa.altura),
    aspect_ratio = :equal,
    framestyle = :box,
    label = "Partícula rastreada",
    color = :crimson,
    lw = 1.5,
    title = "Trajetória de 1 molécula no meio do gás"
)
display(rastro)

# 4. Animação 2D
anim = @animate for k in 1:2:length(tempos)
    scatter(
        pos_x[k, :], pos_y[k, :],
        xlims = (0, caixa.largura),
        ylims = (0, caixa.altura),
        aspect_ratio = :equal,
        framestyle = :box,
        markersize = 10,
        color = cores,
        label = "",
        title = "t = $(round(tempos[k], digits=2)) s",
        size = (450, 450)
    )
end

gif(anim, "caixa_2d.gif", fps=25)