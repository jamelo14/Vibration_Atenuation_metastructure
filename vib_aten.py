import numpy as np
import matplotlib.pyplot as plt

Ndiv = 100
Nper = 100

# Valores iniciais para theta e dtheta
theta0 = np.pi
dtheta0 = 0

# const
k = 1 
m = 1
csi = 0.3
F = 10
omega = .8

# Frequência natural amortecida
wd = np.sqrt(1 - csi**2) * np.sqrt(k/m)

# Arrays de variáveis
c = np.array([k, m, csi, F, omega])
r = np.array([theta0, dtheta0])

# Período
periodo = 2*np.pi/wd

# Tempo final
tf = Nper*periodo

# Número da divisão de onde Poincaré será tirado
n_da_divisao = 4

# Passo de tempo
h = 2*np.pi / (omega * Ndiv)#tf/(Nper*Ndiv)

# cor do gráfico
cor = "red"

# Pode-se plotar a mesma equação variando uma das constantes (omega, omegaN, csi ou F) num dado
# intervalo. Nas especificações abaixo, é possível configurar tal plotagem
# Plot único ou plotar para um intervalo de uma variável
plot_uni = False

# A constante a ser variada.
# Inserir o índice equivalente da lista c
variavel = 4

# Valores que a constante vai assumir (max. de 7 valores)
valores_da_var = [2 ,4, 6, 8, 10, 12, 14]

# Título do gráfico
title = "θ vs t"
# Subtítulo das abscisas
labelx = "tempo [s]"

# Subtítulo das ordenadas
labely = "θ [rad]"

simbolo_da_var = "l"

################################## Inicio do programa ##################################
#    Não fuce em nada que você não sabe o que significa! Criar isso foi trabalhoso!    #
    
def somarNaLista(lista, valor):
	return list(map(lambda x: x + valor, lista))
	
def restoNaLista(lista, valor):
    return lista
    return list(map(lambda x: (abs(x)%valor)*x/abs(x), lista))

tempos = np.zeros(Ndiv*Nper)
thetas = np.zeros(Ndiv*Nper)
dthetas = np.zeros(Ndiv*Nper)

# Matriz da seção de Poincaré
poincare = [[], []]

# Matrizpara o espaço de fase
esp_fas = [[], []]

print(f"h - {h}\nwd - {wd}\ntf - {tf}\ntamanho - {int(tf/h)}")

tempos[0] = 0
thetas[0] = theta0
dthetas[0] = dtheta0

def ddt(t, r, c):
    
    return c[3]*np.cos(c[4]*t) - c[0]/c[1] * np.sin(r[0]) - 2*c[2]*np.sqrt(c[0]/c[1])*r[1] 
    
def dt(t, r, c):
    return r[1]
    
def DELTA_rk4(f, t, r, c):
    global h
    
    k1 = f(t, r, c)
    k2 = f(t + h/2, somarNaLista(r, h/2*k1), c)
    k3 = f(t + h/2, somarNaLista(r, h/2*k2), c)
    k4 = f(t + h    , somarNaLista(r, h*k3)    , c)
    
    return h/6*(k1 + 2*k2 + 2*k3 + k4)
    
def main(c, cor):
    global theta0, dtheta0, labelx, labely, titulo, r, Nper, n_da_divisao, poincare, esp_fas
    
    t = 0
    i = 0
    
    var = r.copy()
     
    var_old =[0,0]     
     
    for periodo in range(0, Nper):
        print(f"Calculando o {periodo + 1}º período")
        for divisao in range(0, Ndiv):
            tempos[i] = t
            
            thetas[i] = var[0]
            dthetas[i] = var[1]
            
            var_old[0] = var[0]
            var_old[1] = var[1]
          
            var[0] += DELTA_rk4(dt, t, var_old, c)
            var[1] += DELTA_rk4(ddt, t, var_old, c)     	
            
            if periodo > (0.8*Nper):
                esp_fas[0].append(thetas[i])
                esp_fas[1].append(dthetas[i])
                
                if divisao == n_da_divisao:
                    poincare[0].append(thetas[i])
                    poincare[1].append(dthetas[i])
                
            t += h
            i+=1
            
    plt.scatter(tempos, restoNaLista(thetas, np.pi), color=cor, s=2)
    
    return tempos, thetas, dthetas
    
def itera(c):
    global variavel, tf, h, tempos, valores_da_var
    
    minimos_theta = []
    maximos_theta = []
    minimos_dtheta = []
    maximos_dtheta = []
    
    numero_de_iteracoes = len(valores_da_var)
    
    i = 0
    
    graficos = []
    
    while i < numero_de_iteracoes:
        print(f"{i}ª iteração da variável {variavel}: valor de {variavel} = {valores_da_var[i]}")
        c[variavel] = valores_da_var[i]
        
        tempos, thetas, dthetas = main(c, None)

        minimos_theta.append(min(thetas))
        minimos_dtheta.append(min(dthetas))
        
        maximos_theta.append(max(thetas))
        maximos_dtheta.append(max(dthetas))
        
        i+=1
    
    legendas = list(map(lambda y: simbolo_da_var + ' = ' + str(y), valores_da_var))  
    plt.legend(legendas)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(title)
    plt.xticks(tempos[::int(1/h)])
    plt.yticks([min(minimos_theta), max(maximos_theta), 0, min(minimos_theta) - 1, max(maximos_theta) + 1])
    plt.axis([0, None, min(minimos_theta) - 1 , max(maximos_theta) + 1])
    plt.show()
    
        
def plot_unico(c, cor):
    t, thetas, dthetas= main(c, cor)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(title)
    plt.xticks(tempos[0::int(1/h)])
    plt.yticks([min(thetas), max(thetas), 0, min(thetas) - 1, max(thetas) + 1])
    plt.axis([0, None, min(thetas) - 1 , max(thetas) + 1])

    plt.show()    
    
    return

if plot_uni == True:
    plot_unico(c, cor)
    
if plot_uni == False:
    itera(c)

plt.scatter(esp_fas[0], esp_fas[1], s=2, color="black")
plt.scatter(poincare[0], poincare[1], color="red", marker=",")

plt.show()
