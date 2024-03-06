import numpy as np
import matplotlib.pyplot as plt

# Параметры задачи
L = 1.0  # Длина области
T = 0.5  # Время интегрирования
Nx = 100  # Количество узлов по координате
Nt = 200  # Количество временных шагов
c=0.9
a = 1.0  # Параметр a
b = 1.0  # Параметр b

# Инициализация сетки
x = np.linspace(0, L, Nx + 1)
V = np.zeros(Nx + 1)
B = np.zeros(Nx + 1)
V_new = np.zeros(Nx + 1)
B_new = np.zeros(Nx + 1)

# Начальные условия (прямоугольные импульсы)
V[x >= 0.4] = 1.0
B[x <= 0.6] = 0.5

# Решение методом Лакса-Вендроффа
for n in range(Nt):
    for i in range(1, Nx):
        V_new[i] = 0.5 * (V[i + 1] + V[i - 1]) - 0.5 * a * c * (B[i + 1] - B[i - 1])
        B_new[i] = 0.5 * (B[i + 1] + B[i - 1]) - 0.5 * b * c * (V[i + 1] - V[i - 1])
    V[:] = V_new
    B[:] = B_new

# Графики решений
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(x, V, label="V")
plt.xlabel("x")
plt.ylabel("V")
plt.title("Решение dV/dt = a*dB/dx")
plt.grid()
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(x, B, label="B", color="orange")
plt.xlabel("x")
plt.ylabel("B")
plt.title("Решение dB/dt = b*dV/dx")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()