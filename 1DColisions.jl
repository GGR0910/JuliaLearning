using Plots

struct Trilho
    inicio::Float64
    fim::Float64
end

mutable struct Particula
    x::Float64
    v::Float64
    m::Float64
    r::Float64
    cnt::Int64
end

function mover!(p::Particula, dt::Real)
    p.x += p.v * dt
end
function colidir!(p1::Particula, p2::Particula)
    dx = p2.x - p1.x
    dv = p2.v - p1.v

    encostadas = abs(dx) <= (p1.r + p2.r)
    aproximando = (dx * dv) < 0

    if encostadas && aproximando
        m1, m2 = p1.m, p2.m
        v1, v2 = p1.v, p2.v

        p1.v = ((m1 - m2) * v1 + 2 * m2 * v2) / (m1 + m2)
        p2.v = ((m2 - m1) * v2 + 2 * m1 * v1) / (m1 + m2)

        p1.cnt += 1
        p2.cnt += 1
    end
end

# Verifica todos os pares possíveis de partículas (funciona para 2 ou N)
function verificar_colisoes!(particulas::Vector{Particula})
    n = length(particulas)
    for i in 1:n
        for j in (i + 1):n
            colidir!(particulas[i], particulas[j])
        end
    end
end

function colidir_parede!(p::Particula, trilho::Trilho)
    # Parede da esquerda
    if p.x - p.r <= trilho.inicio && p.v < 0
        p.v = -p.v
        p.x = trilho.inicio + p.r   # encosta exatamente na parede sem vazar para fora
        # (se quiser contar batidas na parede no seu contador, incremente aqui)
    end

    # Parede da direita
    if p.x + p.r >= trilho.fim && p.v > 0
        p.v = -p.v
        p.x = trilho.fim - p.r
    end
end

# Auditoria física: embalo (P) e energia cinética (Ec)
momento_total(ps) = sum(p.m * p.v for p in ps)
energia_cinetica(ps) = sum(0.5 * p.m * p.v^2 for p in ps)

# 1. Configuração: P1 andando, P2 parada (massas iguais)
trilho = Trilho(0.0, 10.0)
particulas = [
    Particula(1.0, 1.0, 80, 0.3, 0),   # P1: vem da esquerda a +2.0
    Particula(5.0, 0.0, 1.0, 0.3, 0),   # P2: parada no meio
    Particula(8.0, -3.0, 4.0, 0.3, 0)   # P3: vem da direita a -1.0
]
cores = [:dodgerblue, :crimson, :orange]

# dt menor deixa o momento exato do toque mais preciso
dt = 0.05
passos = 800   # 80 * 0.05 = 4.0 segundos

# Auditoria ANTES da simulação
P_inicial = momento_total(particulas)
E_inicial = energia_cinetica(particulas)

# 2. Simulação
n_particulas = length(particulas)
tempos = Float64[]
posicoes = Matrix{Float64}(undef, passos + 1, n_particulas)

for passo in 0:passos
    push!(tempos, passo * dt)
    for (i, p) in enumerate(particulas)
        posicoes[passo + 1, i] = p.x
    end

    for p in particulas
        mover!(p, dt)
        verificar_colisoes!(particulas)
        colidir_parede!(p, trilho)
    end
    verificar_colisoes!(particulas)
end

# Auditoria DEPOIS da simulação
println("Momento  -> antes: $(round(P_inicial, digits=4)) | depois: $(round(momento_total(particulas), digits=4))")
println("Energia  -> antes: $(round(E_inicial, digits=4)) | depois: $(round(energia_cinetica(particulas), digits=4))")
println("Velocidades finais -> P1: $(round(particulas[1].v, digits=2)) | P2: $(round(particulas[2].v, digits=2)) | P3: $(round(particulas[3].v, digits=2))")
println("Número de colisões: $(particulas[1].cnt) | $(particulas[2].cnt) | $(particulas[3].cnt)")

# 3. Gráfico Posição × Tempo
grafico = plot(
    tempos, posicoes,
    xlabel = "Tempo (s)",
    ylabel = "Posição (x)",
    title = "Colisão Elástica 1D",
    label = ["Partícula 1" "Partícula 2"],
    color = [:dodgerblue :crimson],
    lw = 2,
    ylims = (trilho.inicio, trilho.fim)
)
display(grafico)

# 4. Animação (pegando 1 a cada 2 quadros com 1:2:end para renderizar mais rápido)
anim = @animate for k in 1:2:length(tempos)
    t = tempos[k]
    xs = posicoes[k, :]
    ys = zeros(n_particulas)

    plot(
        [trilho.inicio, trilho.fim], [0.0, 0.0],
        color = :gray,
        lw = 3,
        label = "",
        xlims = (trilho.inicio - 0.5, trilho.fim + 0.5),
        ylims = (-1.0, 1.0),
        yaxis = false,
        xlabel = "Posição no trilho (x)",
        title = "t = $(round(t, digits=2)) s",
        size = (600, 200)
    )
    vline!([trilho.inicio, trilho.fim], color = :black, lw = 4, label = "")
    scatter!(
        xs, ys,
        markersize = 12,
        color = cores,
        label = ""
    )
end

gif(anim, "colisao_1d.gif", fps=15)