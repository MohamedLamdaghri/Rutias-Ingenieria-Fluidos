#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cálculo del caudal Q conocida la pérdida de carga Hf

"""

import math
import tkinter as tk
from tkinter import messagebox

# -------------------- CONSTANTES --------------------
G = 9.81
PI = math.pi


# -------------------- FUNCIONES DE CÁLCULO --------------------

def area(D):
    return PI * (D ** 2) / 4


def reynolds(Q, D, nu):
    return (Q / area(D)) * D / nu


def colebrook(Re, k, D, f0=0.02, tol=1e-8, max_iter=100):
    """Ecuación de Colebrook-White."""
    if Re < 2000:
        return 64 / Re, 0  # flujo laminar
    
    f = f0
    for i in range(max_iter):
        f_new = (1 / (-2 * math.log10((k / (3.7 * D)) + (2.51 / (Re * math.sqrt(f)))))) ** 2
        if abs(f_new - f) / f_new < tol:
            return f_new, i + 1
        f = f_new
    return f, max_iter


def swamee_jain_f(Re, k, D):
    """Estimación explícita de f por Swamee–Jain (1976)."""
    if Re < 2000:
        return 64 / Re
    return 0.25 / (math.log10((k / (3.7 * D)) + (5.74 / (Re ** 0.9)))) ** 2


def q_from_hf(f, D, L, Hf, g=G):
    """Despeje de Darcy-Weisbach para Q."""
    num = (PI ** 2) * g * Hf * (D ** 5)
    den = 8 * f * L
    return math.sqrt(abs(num) / den)


def resolver_Q(D, L, k, nu, Hf, tol=1e-6, max_iter=200):
    """Itera Q -> Re -> f -> Q hasta converger."""
    
    # Paso 1: estimación inicial de Q con f≈0.02 para calcular Re
    f = 0.02
    Q = q_from_hf(f, D, L, Hf)
    Re_est = reynolds(Q, D, nu)
    
    # Paso 2: inicializar f con Swamee–Jain (más eficiente)
    f = swamee_jain_f(Re_est, k, D)
    
    # Paso 3: iterar hasta convergencia
    for n in range(max_iter):
        V = Q / area(D)
        Re = reynolds(Q, D, nu)
        f, _ = colebrook(Re, k, D, f)
        
        # Ajuste por pérdidas menores
        
        Q_new = q_from_hf(f, D, L, Hf, G)
        err = abs(Q_new - Q) / Q_new
        if err < tol:
            Q = Q_new
            break
        Q = Q_new
    
    # Paso 4: resultados finales
    V = Q / area(D)
    Re = reynolds(Q, D, nu)
    
    if Re < 2000:
        regimen = "Laminar"
    elif Re < 4000:
        regimen = "Transitorio"
    else:
        regimen = "Turbulento"
    
    return {
        "Q": Q,
        "V": V,
        "f": f,
        "Re": Re,
        "regimen": regimen,
        "iteraciones": n + 1,
        "converge": err < tol
    }


# -------------------- INTERFAZ GRÁFICA --------------------

def calcular():
    try:
        D = float(entry_D.get())
        L = float(entry_L.get())
        k = float(entry_k.get())
        nu = float(entry_nu.get())
        Hf = float(entry_Hf.get())

        res = resolver_Q(D, L, k, nu, Hf)

        texto = f"""
===== RESULTADOS =====
Caudal Q        = {res['Q']:.8f} m³/s
Velocidad V     = {res['V']:.6f} m/s
Número de Re    = {res['Re']:.2f}
Factor f        = {res['f']:.6f}
Régimen         = {res['regimen']}
Iteraciones     = {res['iteraciones']}
Convergencia    = {'Sí' if res['converge'] else 'No'}
=======================
"""
        text_result.delete("1.0", tk.END)
        text_result.insert(tk.END, texto.strip())

    except ValueError as e:
        messagebox.showerror("Error de entrada", str(e))
    except Exception as e:
        messagebox.showerror("Error inesperado", str(e))


# Crear ventana principal
root = tk.Tk()
root.title("Cálculo de caudal Q conocida Hf (Colebrook–White mejorado)")
root.geometry("520x520")
root.resizable(False, False)

# ----- ENTRADAS -----
tk.Label(root, text="Datos de entrada (unidades SI)", font=("Arial", 12, "bold")).pack(pady=5)

frame_inputs = tk.Frame(root)
frame_inputs.pack(pady=5)

labels = ["Diámetro D [m]", "Longitud L [m]", "Rugosidad k [m]", "Viscosidad ν [m²/s]", "Pérdida de carga Hf [m]"]
defaults = ["0.03", "1000", "0.001", "0.001", "-41.743"]
entries = []

for i, (lab, val) in enumerate(zip(labels, defaults)):
    tk.Label(frame_inputs, text=lab, anchor="w", width=25).grid(row=i, column=0, pady=3, sticky="w")
    e = tk.Entry(frame_inputs, width=15)
    e.insert(0, val)
    e.grid(row=i, column=1, pady=3)
    entries.append(e)

entry_D, entry_L, entry_k, entry_nu, entry_Hf = entries

# ----- BOTÓN DE CÁLCULO -----
tk.Button(root, text="Calcular Q", command=calcular, bg="#4CAF50", fg="white", font=("Arial", 12, "bold"), width=20).pack(pady=10)

# ----- RESULTADOS -----
tk.Label(root, text="Resultados", font=("Arial", 12, "bold")).pack(pady=5)
text_result = tk.Text(root, height=12, width=55, font=("Courier New", 10))
text_result.pack(pady=5)

# ----- EJECUCIÓN -----
root.mainloop()
