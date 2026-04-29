#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from acritico import funcion_fx


# El nombre al pdf sólo se asigna con la variable REACTOR, ojo
# simetrica = {'H':54, 'H0':54, 'Hr':54, 'd':0}
# REACTOR = "simetrica"
REACTOR = "RA3"

fx = lambda x: funcion_fx(x, REACTOR)
#fx = lambda x: funcion_fx(x, datos=simetrica)

# --- Configuración para A4 Apaisado (11.69 x 8.27 pulgadas) ---
fig, ax = plt.subplots(figsize=(11.69, 8.27))

# --- Eje X: Extracción de la barra [%] ---
xticks_labels = [0, 6, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 
                 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68,
                 70, 72, 74, 78, 82, 88, 100]
xticks_positions = fx(np.array(xticks_labels))

# Ticks menores cada 1 unidade para las líneas verticales de referencia
minor_xticks_t = np.arange(0, 101, 1)
minor_xticks_pos = fx(minor_xticks_t)

ax.set_xticks(xticks_positions)
ax.set_xticklabels(xticks_labels, fontsize=7) # Tamaño de fuente 
ax.set_xticks(minor_xticks_pos, minor=True)

# --- Eje Y: Inversa del contaje normalizado ---
ax.set_yticks(np.arange(0, 1.01, 0.04))
ax.set_yticks(np.arange(0, 1.01, 0.01), minor=True)

# --- Estilo de la Grilla (Papel milimetrado) ---
ax.grid(True, which='major', color='black', linestyle='-', linewidth=0.7)
ax.grid(True, which='minor', color='black', linestyle='-', linewidth=0.3, alpha=0.4)

# --- Refuerzo de los bordes del gráfico (Ejes externos) ---
for spine in ax.spines.values():
    spine.set_linewidth(0.7)
    spine.set_color('black')
    spine.set_visible(True)

ax.set_facecolor('white')

# Etiquetas de ejes
ax.set_xlabel('Extracción de la barra [%]', fontsize=11, labelpad=10)
ax.set_ylabel('Inversa del contaje normalizado', fontsize=11, labelpad=10)
ax.set_xlim(fx(np.array([0])), fx(np.array([100])))
ax.set_ylim(0, 1)

# Encabezados superiores (posicionamiento relativo para A4)
plt.text(0.0, 1.04, 'Detectores graficados:', transform=ax.transAxes, fontsize=10, fontweight='bold')
plt.text(0.40, 1.04, 'Barra extraída:', transform=ax.transAxes, fontsize=10, fontweight='bold')
plt.text(0.90, 1.04, f'Reactor: {REACTOR}', transform=ax.transAxes, fontsize=10, fontweight='bold')

# --- Ajustes de márgenes para maximizar espacio ---
plt.subplots_adjust(left=0.07, right=0.97, top=0.90, bottom=0.12)

# Guardado en PDF
plt.savefig(f"grafico_acritico_{REACTOR}.pdf", orientation='landscape')
plt.show()
