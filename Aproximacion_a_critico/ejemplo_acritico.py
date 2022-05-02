#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

"""
RECORDAR: actualizar con la últmia versión del repositorio Cinepy
"""
# import sys
# sys.path.append("/home/pablo/CinePy")
# from modules.point_kinetics.io_acritico import read_acritico

from io_acritico import read_acritico


# Archivos que se van a leer
archivos = [
            "tasa_tp3.D1.bin",
            "tasa_tp3.D2.bin",
            ]

fig, ax = plt.subplots()
for archivo in archivos:
    # Leo el archivo
    data, t, dt, header = read_acritico(archivo)
    # Convierto a cps
    R = data / dt
    # Grafico
    ax.plot(t, R, label= "Detector " + archivo.split('.')[-2])


ax.set_yscale('log')
ax.set_xlabel('Tiempo [s]')
ax.set_ylabel("Cuentas [cps]")
ax.legend()

plt.tight_layout()
plt.show()
