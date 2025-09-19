import numpy as np
import matplotlib.pyplot as plt

def cot(x):
    return 1 / np.tan(x)

def air_velocity(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    Vr = -np.sin(theta)
    Vt = (1 / (2 * np.pi * r)) + np.cos(theta)

    Vx = Vr * np.cos(theta) - Vt * np.sin(theta)
    Vy = Vr * np.sin(theta) + Vt * np.cos(theta)
    return Vx, Vy

def generar_curva(C, signo=1, n_puntos=5000):
    t = np.linspace(np.deg2rad(15.5568), np.deg2rad(176.891), n_puntos)
    if signo == 1:
        x = (5 * (C - t) * cot(t)) / (2*np.pi)
        y = 5 * (C - t) / (2*np.pi)
    else:
        x = (5 * (C - t) * cot(t)) / (2*np.pi)
        y = -5 * (C - t) / (2*np.pi)
    return np.column_stack((x, y))

# Curvas
C1 = 3.7699111843077517
C2 = 4.39822971502571
curva1 = generar_curva(C1, signo=1)
curva2 = generar_curva(C2, signo=-1)
curva3 = generar_curva(np.pi, signo=-1)
curva4 = generar_curva(np.pi, signo=1)
paredes = [curva1, curva2, curva3, curva4]

# Parámetros de la simulación
dt = 0.1
steps = 30000
enormal = 0.95
etangencial = 0.70

# Posición inicial
x, y = -8.0, 0.0
vx, vy = 1.0, 0.0  # velocidad inicial

trayectoria_x = [x]
trayectoria_y = [y]

for _ in range(steps):
    # Velocidad del aire en la posición actual
    Vx_air, Vy_air = air_velocity(x, y)

    # Aquí decides:
    # O bien haces que la partícula siga exactamente al aire:
    # vx, vy = Vx_air, Vy_air
    #
    # O bien haces una mezcla: la partícula se ajusta al aire poco a poco
    # (más realista, si quieres simular inercia):
    alpha = 0.1  # qué tan rápido se ajusta al aire
    vx += alpha * (Vx_air - vx)
    vy += alpha * (Vy_air - vy)

    # Nueva posición
    x_new = x + vx * dt
    y_new = y + vy * dt

    # Verificar colisiones con paredes
    for curva in paredes:
        distancias = np.sqrt((curva[:, 0] - x_new) ** 2 + (curva[:, 1] - y_new) ** 2)
        idx_min = np.argmin(distancias)
        if distancias[idx_min] < 0.1:  # umbral de choque
            # Punto más cercano
            px, py = curva[idx_min]

            # Vector tangente aproximado (derivada local)
            if 1 < idx_min < len(curva) - 2:
                tx = curva[idx_min + 1, 0] - curva[idx_min - 1, 0]
                ty = curva[idx_min + 1, 1] - curva[idx_min - 1, 1]
            else:
                tx, ty = 1, 0
            tangente = np.array([tx, ty])
            tangente = tangente / np.linalg.norm(tangente)

            # Normal = perpendicular a la tangente
            normal = np.array([-tangente[1], tangente[0]])

            # Descomponer velocidad
            v = np.array([vx, vy])
            v_t = np.dot(v, tangente) * tangente
            v_n = np.dot(v, normal) * normal

            # Aplicar coeficientes
            v_rebote = etangencial * v_t - enormal * v_n
            vx, vy = v_rebote

            # Empujar ligeramente fuera de la pared
            x_new = x + vx * dt + normal[0] * 1e-4
            y_new = y + vy * dt + normal[1] * 1e-4

            break

    # Guardar trayectoria
    trayectoria_x.append(x_new)
    trayectoria_y.append(y_new)

    x, y = x_new, y_new

# ---------------------------
# Graficar
# ---------------------------
plt.figure(figsize=(8,6))
plt.plot(curva1[:,0], curva1[:,1], 'k')
plt.plot(curva2[:,0], curva2[:,1], 'k')
plt.plot(curva3[:,0], curva3[:,1], 'k')
plt.plot(curva4[:,0], curva4[:,1], 'k')
plt.scatter(trayectoria_x, trayectoria_y, c='r', s=10, label="Trayectoria")

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Trayectoria de una partícula con rebotes")
plt.legend()
plt.axis("equal")
plt.show()



