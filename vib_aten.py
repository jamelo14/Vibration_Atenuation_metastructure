import numpy as np
import matplotlib.pyplot as plt

Ndiv = 100
Nper = 100

# Valores iniciais para theta e dtheta
theta0 = 0
dtheta0 = 0

# const
k = 9
m = 2
csi = 0.5
F = 12
omega = 12

h = 0.01

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
valores_da_var = np.linspace(2, 4, 100) #[2 ,4, 6, 8, 10, 12, 14]

# Título do gráfico
title = "θ vs t"
# Subtítulo das abscisas
labelx = "tempo [s]"

# Subtítulo das ordenadas
labely = "θ [rad]"

simbolo_da_var = "l"

################################## Inicio do programa ##################################
#    Não fuce em nada que você não sabe o que significa! Criar isso foi trabalhoso!    #

minimos_theta = []
maximos_theta = []
minimos_dtheta = []
maximos_dtheta = []

def somarNaLista(lista, valor):
	return list(map(lambda x: x + valor, lista))
	
def restoNaLista(lista, valor):
    return lista
    #return list(map(lambda x: (abs(x)%valor)*x/abs(x), lista))

tempos = np.zeros(Ndiv*Nper)
thetas = np.zeros(Ndiv*Nper)
dthetas = np.zeros(Ndiv*Nper)

# Matriz da seção de Poincaré
poincare = [[], []]

# Matriz para o espaço de fase
esp_fas = [[], []]

tempos[0] = 0
thetas[0] = theta0
dthetas[0] = dtheta0

def ddt(t, r, c):
    return c[3]*np.cos(c[4]*t) - c[0]/c[1] * r[0] - 2*c[2]*np.sqrt(c[0]/c[1])*r[1] 
    
def dt(t, r, c):
    return r[1]
    
def DELTA_rk4(f, t, r, c):
    global h
    
    k1 = f(t, r, c)
    k2 = f(t + h/2, somarNaLista(r, h/2*k1), c)
    k3 = f(t + h/2, somarNaLista(r, h/2*k2), c)
    k4 = f(t + h    , somarNaLista(r, h*k3)    , c)
    
    return (h/6*(k1 + 2*k2 + 2*k3 + k4))
    
def plotPoincare():
    plt.scatter(esp_fas[0], esp_fas[1], s=2, color="black")
    plt.scatter(poincare[0], poincare[1], color="red", marker=",")
    plt.xlabel("theta [rad]")
    plt.ylabel("dtheta [rad/s]")
    plt.xticks(  list(map( lambda x: round(x, 1), [min(esp_fas[0]), max(esp_fas[0]), ( min(esp_fas[0]) + max(esp_fas[0]) )/2, min(esp_fas[0]) - 2, max(esp_fas[0])+2 ] )) )
    plt.yticks(  list(map( lambda x: round(x, 1), [min(esp_fas[1]), max(esp_fas[1]), ( min(esp_fas[1]) + max(esp_fas[1]) )/2, min(esp_fas[1]) - 2, max(esp_fas[1])+2 ] )) )
    plt.axis( list(map( lambda x: round(x, 1),  [min(esp_fas[0]) - 2, max(esp_fas[0])+2, min(esp_fas[1]) - 2, max(esp_fas[1])+2 ] )) )
    plt.title("Espaço de fases e seção de Poincaré")
    plt.show()    

def plotVib():
    plt.plot(valores_da_var, maximos_theta) # Plota omega x theta max
    plt.show()

def main(c, cor):
    global theta0, dtheta0, labelx, labely, titulo, r, Nper, n_da_divisao, poincare, esp_fas
    
    temposNT= []
    thetasNotTrans = []
    dthetasNotTrans = []
    
    t = 0
    i = 0
        
    var = [r[0], r[1]]
     
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
 
          
            if periodo > (0.9*Nper):
                esp_fas[0].append(thetas[i])
                esp_fas[1].append(dthetas[i])
                
                temposNT.append(t)
                thetasNotTrans.append(thetas[i])
                dthetasNotTrans.append(dthetas[i])
                
                if divisao == n_da_divisao:
                    poincare[0].append(thetas[i])
                    poincare[1].append(dthetas[i])
                
            t += h
            i+=1
            
    print(len(temposNT))
    print(len(thetasNotTrans))
    
    return temposNT, thetasNotTrans, dthetasNotTrans
    
def itera(c):
    global variavel, tf, h, tempos, valores_da_var
    
    times =[]
    t = []
    dt = []
    
    numero_de_iteracoes = len(valores_da_var)
    
    i = 0
    
    graficos = []
    
    while i < numero_de_iteracoes:
        
        print(f"{i}ª iteração da variável {variavel}: valor de {variavel} = {valores_da_var[i]}")
        
            
        # Passo de tempo
        # OBS.: Entender melhor o motivo de ser assim
        h = 2*np.pi / ( valores_da_var[i] * Ndiv)
        
        c[variavel] = valores_da_var[i]
        
        times, t, dt = main(c, None)

        minimos_theta.append(min(t))
        minimos_dtheta.append(min(dt))
        
        maximos_theta.append(max(t))
        maximos_dtheta.append(max(dt))
        
        #plt.plot(times, t)
        #plt.show()
        
        i+=1
    
        
def plot_unico(c, cor):
    t, thetasNT, dthetasNT= main(c, cor)
    
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(title)
    plt.xticks( list(map( lambda x: round(x, 2), tempos[0:-1:int( (Ndiv*Nper)/6 )] )) )
    plt.yticks([min(thetas), max(thetas), 0, min(thetas) - 1, max(thetas) + 1])
    plt.axis([0, None, min(thetas) - 1 , max(thetas) + 1])
    
    plt.plot(tempos, thetas, color=cor, linewidth=2)
    plt.xlabel("s")
    
    plt.show()    
    
    return

if plot_uni == True:
    plot_unico(c, cor)
    
if plot_uni == False:
    itera(c)
    
    plotVib()

#plotPoincare()
