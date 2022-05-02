#! /usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def funcion_fx(z, reactor=None, datos=None):
    """
    Curva de calibración en reactividad de una barra de control normalizada.

    Se obtiene usando teoría de perturbacioens a primer orden a un grupo de
    energía. Se utiliza una reactor homogéneo desnudo con una barra central.

    Referencia: ITE-06REC-203

    Parámetros
    ----------
        z: float or numpy array
            Valor/es en donde se quiere evaluar a la f(x)
        reactor: string ('RA1', 'RA3')
            Se cargan los datos geométricos de los reactores
        datos: dict
            Se especifican los datos manualmente. Estos son:
                datos['H']: Altura del núcleo [cm]
                datos['H0']: Longitud de la barra [cm]
                datos['Hr']: Recorrido de la barra [cm]
                datos['d']: Ahorro por reflector [cm]
            Si se especifican 'datos' no se tienen en cuenta los datos de
            'reactor'.

    Resultados
    ----------
        fx: float or numpy array (igual que z)
            Valores de f(x) evaluados en z
    """
    if (reactor is None) and (datos is None):
       raise ValueError("No se epsecifican datos del problema")
    try:
        len(z)
    except TypeError:
        z = [z]
    z = np.asarray(z)

    if datos is None:
        # Por si se escribe al nombre del reactor con guión
        reactor = reactor.replace('-', '')
        if reactor == 'RA1':
            datos = {
                     'H': 54,
                     'H0': 54,
                     'Hr': 63,
                     'd': 8,
                     }
        elif reactor == 'RA3':
            datos = {
                     'H': 61.5,
                     'H0': 61.5,
                     'Hr': 66,
                     'd': 8,
                     }
        else:
            raise ValueError(f"No se tienen datos del reactor {reactor}")

    H = datos['H']
    H0 = datos['H0']
    Hr = datos['Hr']
    d = datos['d']

    # Normalización
    z = z * Hr / 100

    fx = np.zeros(len(z))

    A1 = np.zeros(len(z))
    A2 = np.zeros(len(z))
    A3 = np.sin(np.pi * (2 * z - H0) / (H + 2 * d))
    A4 = np.zeros(len(z))
    A4[:] = 2 * np.pi * H0 / (H + 2 * d) + 2 * np.sin(np.pi * H0 / (H + 2 * d))

    _tmp = (H - H0) / 2 + d

    ind1 = z <= _tmp
    A1[ind1] = 2 * np.pi * H0 / (H + 2 * d)
    A2[ind1] = np.sin(np.pi * (H0 + 2 * z[ind1]) / (H + 2 * d))
    fx[ind1] = (A4[ind1] + A3[ind1] - A2[ind1] - A1[ind1]) / A4[ind1]

    ind2 = (z > _tmp) & (z < (H0 + _tmp))
    A1[ind2] = np.pi * (H0 + H + 2 * d - 2 * z[ind2]) / (H + 2 * d)
    A2[ind2] = 0.0
    fx[ind2] = (A4[ind2] + A3[ind2] - A2[ind2] - A1[ind2]) / A4[ind2]

    ind3 = z >= (H0 + _tmp)
    fx[ind3] = 1

    return fx


def inversa_fx(fx_val, fx):
    """
    Invierte a la función fx en los puntos fx_val
            x0 = inversa(fx)(fx_val)

    Parámetros
    ----------
        fx_val: float or list of floats
            Valores en los cuales se quiere invertir a la fx
        fx: function
            Función fx que se quiere invertir

    Resultado
    ---------
        x0: float or numpy array (lo mismo que fx_val)
            Valores de la inversa de fx en fx_val
    """

    from scipy.optimize import brentq
    # Para aceptar tanto escalares como arrays
    try:
        len(fx_val)
    except TypeError:
        fx_val = [fx_val]

    x0 = np.zeros(len(fx_val))
    for i, y in enumerate(fx_val):
        fun = lambda x: fx(x) - y
        x0[i] = brentq(fun, 0, 100)

    return x0


def grafica_acritico(x, y_data, fx, y_err=None):
    """
    Grafica los datos de la aproximación a crítico.

    Si se especifica y_err se hace el gráfico con barras de error.

    Parámetros
    ----------
        x: list of floats
            Porcentajes de extracción de barras de control
        y: list of list of floats
            Cada elemento de la lista es una lista con los valores de la
            inversa del contaje normalizado.
            Cada elemento de la lista corresponde a un detector.
        fx: function
            Función f(x) utilizada
        y_err: list of list of floats (opcional)
            Cada elemento de la lista es una lista con las incertezas asociadas
            a la inversa del contaje normalizado.
            Cada elemento de la lista corresponde a un detector.

    Resultados
    ----------
        None

    """
    x = np.asarray(x)
    if y_err is not None:
        y_err = np.asarray(y_err)

    fig, ax = plt.subplots(figsize=(15, 7))

    if y_err is None:
        # Gráfico sin incertezas
        for j, y in enumerate(y_data):
            y = np.asarray(y)
            ax.plot(fx(x), y, 's', lw=1, ms=7, label=f"Detector {j+1}")
    else:
        # Gráfico con incertezas
        for j, y in enumerate(y_data):
            y = np.asarray(y)
            ax.errorbar(fx(x), y, yerr=y_err[j], label=f"Detector {j+1}",
                        fmt='s', elinewidth=1, capsize=3)

    ax.legend(loc='upper right')
    # Creo los ticks para el eje x
    xticks = np.asarray([0, 8, 12])
    xticks = np.append(xticks, np.arange(12, 78, 2))
    xticks = np.append(xticks, np.asarray([78, 82, 100]))
    ax.set_xticks(fx(xticks))
    ax.set_xticklabels(xticks)
    ax.set_xlim([0, 1])
    # Creo los ticks para el eje y
    n_yticks = 20
    yticks = np.arange(0, n_yticks + 1) / n_yticks
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{s:.2f}" for s in yticks])
    ax.set_ylim([0, 1])

    ax.set_xlabel("Procentaje de extracción [%]")
    ax.set_ylabel("R$_o$ / R")

    fig.tight_layout()
    plt.show()
    return None


if __name__ == "__main__":
    ###########################################################################
    # Función f(x) 
    ###########################################################################
    # Posiciones de barras
    test = [0, 22, 41, 50, 100]
    print(f"Porcentajes en donde se determina la f(x): {test}")
    fx_test = funcion_fx(test, "RA1")
    print(f"Valores de f(test): {fx_test}")
    my_fun = lambda x: funcion_fx(x, "RA1")
    test_rec = inversa_fx(fx_test, my_fun)
    print(f"Valores invirtiendo la f(x): {test_rec}")
    # Gráfico fx
    fig, ax = plt.subplots()
    x = np.linspace(0, 100, 200)
    ax.plot(x, funcion_fx(x, 'RA1'), label="f(x) del RA1")
    datos_simetrica = {'H':54, 'H0':54, 'Hr':54, 'd':0}
    ax.plot(x, funcion_fx(x, datos=datos_simetrica), label="f(x) simétrica")
    ax.set_xlabel("Porcentaje de extracción [%]")
    ax.set_ylabel("f(x)")
    ax.legend(loc='lower right')
    fig.tight_layout()

    ###########################################################################
    # Graficación aproximación a crítico
    ###########################################################################
    # Posiciones de la barra de control
    bcs = [0, 10, 30]
    # Tasa de cuentas normalizada del detector 1
    tasa1 = [1, 0.5, 0.4]
    # Incerteza del detector 1
    tasa1_err = [0.01, 0.01, 0.02]
    # Tasa de cuentas normalizada del detector 2
    tasa2 = [1, 0.6, 0.2]
    # Incerteza del detector 2
    tasa2_err = [0.01, 0.01, 0.02]
    grafica_acritico(bcs, [tasa1, tasa2], my_fun, y_err=[tasa1_err, tasa2_err])

