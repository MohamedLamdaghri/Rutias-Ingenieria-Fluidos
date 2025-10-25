#!/usr/bin/env python3
import math
import tkinter as tk
from tkinter import ttk, messagebox

def darcy_friction_factor(Re: float, D: float, k: float, error: float = 1e-8):
    if Re <= 0:
        raise ValueError("Re debe ser > 0")
    if D <= 0:
        raise ValueError("D debe ser > 0")
    if k < 0:
        raise ValueError("k debe ser ≥ 0")
    if error <= 0:
        raise ValueError("error (tolerancia) debe ser > 0")

    if Re < 2300.0:
        f = 64.0 / Re
        return f, "Laminar (f = 64/Re)"

    rel_roughness = k / D
    Re_limit = float("inf")
    if rel_roughness > 0.0:
        Re_limit = 3500.0 / rel_roughness

    def swamee_jain(Re_, D_, k_):
        term = (k_/(3.7*D_)) + (5.74/(Re_**0.9))
        return 0.25 / (math.log10(term)**2)

    if rel_roughness == 0.0:
        f = 0.02 if Re <= 1e6 else 0.018
        while True:
            denom = 2.0*math.log10(Re*math.sqrt(f)) - 0.8
            f_new = 1.0/(denom*denom)
            er = abs((f_new - f) / f_new)
            f = f_new
            if er <= error:
                break
        method = "Turbulento, liso (Prandtl–Kármán)"
    elif Re < Re_limit:
        f = swamee_jain(Re, D, k)
        while True:
            denom = -2.0*math.log10((rel_roughness/3.71) + (2.51/(Re*math.sqrt(f))))
            f_new = 1.0/(denom*denom)
            er = abs((f_new - f) / f_new)
            f = f_new
            if er <= error:
                break
        method = "Turbulento, transición (Colebrook–White)"
    else:
        if k == 0.0:
            return darcy_friction_factor(Re, D, 0.0, error)
        f = 1.0 / (2.0*math.log10(D/k) + 1.14)**2
        method = "Turbulento, completamente rugoso (Prandtl–Kármán)"
    return f, method

class DarcyApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Darcy–Weisbach | Factor de fricción")
        self.geometry("560x420")
        self.minsize(520, 380)
        self._build_ui()

    def _build_ui(self):
        pad = {"padx": 10, "pady": 8}

        frm = ttk.Frame(self)
        frm.pack(fill="both", expand=True, **pad)

        # Inputs
        self.Re_var = tk.StringVar(value="120000")
        self.D_var  = tk.StringVar(value="0.1")
        self.k_var  = tk.StringVar(value="0.0001")
        self.err_var= tk.StringVar(value="1e-8")

        row = 0
        ttk.Label(frm, text="Reynolds (Re):").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.Re_var, width=20).grid(row=row, column=1, sticky="w")
        ttk.Label(frm, text="(adimensional)").grid(row=row, column=2, sticky="w"); row+=1

        ttk.Label(frm, text="Diámetro (D):").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.D_var, width=20).grid(row=row, column=1, sticky="w")
        ttk.Label(frm, text="m (misma unidad que k)").grid(row=row, column=2, sticky="w"); row+=1

        ttk.Label(frm, text="Rugosidad (k):").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.k_var, width=20).grid(row=row, column=1, sticky="w")
        ttk.Label(frm, text="m (misma unidad que D)").grid(row=row, column=2, sticky="w"); row+=1

        ttk.Label(frm, text="Tolerancia (error):").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.err_var, width=20).grid(row=row, column=1, sticky="w")
        ttk.Label(frm, text="p. ej. 1e-8").grid(row=row, column=2, sticky="w"); row+=1

        # Buttons
        btns = ttk.Frame(frm)
        btns.grid(row=row, column=0, columnspan=3, sticky="w", pady=(12, 4))
        calc_btn = ttk.Button(btns, text="Calcular", command=self.on_calculate)
        calc_btn.grid(row=0, column=0, padx=(0,6))
        clear_btn = ttk.Button(btns, text="Limpiar", command=self.on_clear)
        clear_btn.grid(row=0, column=1, padx=(0,6))
        quit_btn = ttk.Button(btns, text="Salir", command=self.destroy)
        quit_btn.grid(row=0, column=2)

        # Results
        sep = ttk.Separator(frm, orient="horizontal")
        sep.grid(row=row+1, column=0, columnspan=3, sticky="ew", pady=(8,6))

        res = ttk.LabelFrame(frm, text="Resultados")
        res.grid(row=row+2, column=0, columnspan=3, sticky="nsew")
        frm.rowconfigure(row+2, weight=1)
        frm.columnconfigure(2, weight=1)

        self.f_var      = tk.StringVar(value="—")
        self.method_var = tk.StringVar(value="—")
        self.rr_var     = tk.StringVar(value="—")
        self.regime_var = tk.StringVar(value="—")

        grid_r = 0
        ttk.Label(res, text="f (Darcy):").grid(row=grid_r, column=0, sticky="w", padx=8, pady=6)
        ttk.Label(res, textvariable=self.f_var, font=("TkDefaultFont", 11, "bold")).grid(row=grid_r, column=1, sticky="w"); grid_r+=1

        ttk.Label(res, text="Método:").grid(row=grid_r, column=0, sticky="w", padx=8, pady=6)
        ttk.Label(res, textvariable=self.method_var).grid(row=grid_r, column=1, sticky="w"); grid_r+=1

        ttk.Label(res, text="Régimen:").grid(row=grid_r, column=0, sticky="w", padx=8, pady=6)
        ttk.Label(res, textvariable=self.regime_var).grid(row=grid_r, column=1, sticky="w"); grid_r+=1

        ttk.Label(res, text="Rugosidad relativa k/D:").grid(row=grid_r, column=0, sticky="w", padx=8, pady=6)
        ttk.Label(res, textvariable=self.rr_var).grid(row=grid_r, column=1, sticky="w"); grid_r+=1

        # Log / notes
        self.notes = tk.Text(res, height=6, wrap="word")
        self.notes.grid(row=grid_r, column=0, columnspan=2, sticky="nsew", padx=8, pady=(6,8))
        res.rowconfigure(grid_r, weight=1)
        res.columnconfigure(1, weight=1)

        # Bind Enter to calculate
        self.bind("<Return>", lambda e: self.on_calculate())

    def parse_float(self, s: str, name: str) -> float:
        try:
            return float(s.replace(",", "."))
        except Exception:
            raise ValueError(f"'{name}' no es un número válido.")

    def on_clear(self):
        self.f_var.set("—")
        self.method_var.set("—")
        self.rr_var.set("—")
        self.regime_var.set("—")
        self.notes.delete("1.0", "end")

    def on_calculate(self):
        try:
            Re = self.parse_float(self.Re_var.get(), "Re")
            D  = self.parse_float(self.D_var.get(), "D")
            k  = self.parse_float(self.k_var.get(), "k")
            err= self.parse_float(self.err_var.get(), "error")

            f, method = darcy_friction_factor(Re, D, k, err)
            rr = k / D
            regime = "Laminar" if Re < 2300 else "Turbulento"

            self.f_var.set(f"{f:.8f}")
            self.method_var.set(method)
            self.regime_var.set(regime)
            self.rr_var.set(f"{rr:.6g}")

            # Notes
            self.notes.delete("1.0", "end")
            if Re < 2300:
                self.notes.insert("end", "Flujo laminar; f = 64/Re es exacto en esta región.\n")
            else:
                if k == 0.0:
                    self.notes.insert("end", "Se ha asumido tubería lisa (k = 0).\n")
                elif rr < 1e-6:
                    self.notes.insert("end", "Rugosidad relativa muy pequeña; el resultado se aproxima al caso liso.\n")
                if Re < 4000:
                    self.notes.insert("end", "Reynolds cerca de la transición; verifica la validez del modelo.\n")

        except Exception as e:
            messagebox.showerror("Error de entrada", str(e))

if __name__ == "__main__":
    app = DarcyApp()
    app.mainloop()
