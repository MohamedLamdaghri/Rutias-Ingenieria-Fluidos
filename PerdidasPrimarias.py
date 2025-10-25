import tkinter as tk
from tkinter import ttk, messagebox
import math

def darcy_friction_factor(Re, D, k, error):
    """
    Calculates the Darcy-Weisbach friction factor based on:
    - Laminar flow (f = 64/Re)
    - Turbulent smooth pipe (Von Karman–Prandtl)
    - Turbulent rough pipe: transition (Colebrook–White)
    - Turbulent rough pipe: fully rough (Von Karman–Prandtl 2)
    """

    if Re < 2300:
        f = 64.0 / Re
        method = "Laminar flow"
    else:
        rel_roughness = k / D
        Re_limit = 3500.0 / rel_roughness if rel_roughness > 0 else float("inf")

        if rel_roughness == 0:
            # Smooth pipe: Von Karman–Prandtl iterative method
            f = 25 / (math.log10((k/3.7*D)+(5.74/Re**0.9)))**2 # initial guess
            er = 1
            while er > error:
                f_new = 1.0 / (2 * math.log10(Re * math.sqrt(f)) - 0.8)**2
                er = abs((f_new - f) / f_new)
                f = f_new
            method = "Turbulent, smooth pipe (Von Karman–Prandtl)"
        elif Re < Re_limit:
            # Colebrook–White equation, iterative solution
            f = 25 / (math.log10((k/3.7*D)+(5.74/Re**0.9)))**2 # initial guess
            er = 1
            while er > error:
                f_new = 1.0 / (-2.0 * math.log10((rel_roughness / 3.71) + (2.51 / (Re * math.sqrt(f)))))**2
                er = abs((f_new - f) / f_new)
                f = f_new
            method = "Turbulent, transition (Colebrook–White)"
        else:
            f = 1.0 / (2 * math.log10(D / k) + 1.14)**2
            method = "Turbulent, fully rough (Von Karman–Prandtl 2)"
    return f, method

def process_darcy():
    try:
        # Get input values
        Re = float(entry1.get() or 0)
        D = float(entry2.get() or 0)
        k = float(entry3.get() or 0)
        error = float(entry4.get() or 0.0001)
        L = float(entry5.get() or 0)
        Q = float(entry6.get() or 0)
        g = 9.81  # gravitational acceleration [m/s²]

        # Calculate Darcy friction factor
        f, method = darcy_friction_factor(Re, D, k, error)

        # Calculate head loss (Darcy–Weisbach using Q)
        if D > 0:
            hf = (8 * f * L * Q**2) / (math.pi**2 * g * D**5)
        else:
            hf = float('nan')

        # Display results
        output_var.set(f"f = {f:.6f} ({method})    \nHead loss = {hf:.4f} m")

    except ValueError:
        messagebox.showerror("Error", "Please introduce only numbers.")

# --- GUI SETUP ---
root = tk.Tk()
root.title("Darcy-Weisbach Friction Factor and Head Loss (with Flow Q)")
root.geometry("700x500")
root.resizable(True, True)

frame = ttk.Frame(root, padding=10)
frame.pack(fill="both", expand=True)

# Inputs
ttk.Label(frame, text="Re : Reynolds number").pack(pady=2, anchor="w")
entry1 = ttk.Entry(frame)
entry1.pack(pady=2, fill="x")

ttk.Label(frame, text="D : Pipe diameter [m]").pack(pady=2, anchor="w")
entry2 = ttk.Entry(frame)
entry2.pack(pady=2, fill="x")

ttk.Label(frame, text="k : Roughness height [m]").pack(pady=2, anchor="w")
entry3 = ttk.Entry(frame)
entry3.pack(pady=2, fill="x")

ttk.Label(frame, text="Error (e.g. 0.0001)").pack(pady=2, anchor="w")
entry4 = ttk.Entry(frame)
entry4.pack(pady=2, fill="x")

ttk.Label(frame, text="L : Pipe length [m]").pack(pady=2, anchor="w")
entry5 = ttk.Entry(frame)
entry5.pack(pady=2, fill="x")

ttk.Label(frame, text="Q : Volumetric flow rate [m³/s]").pack(pady=2, anchor="w")
entry6 = ttk.Entry(frame)
entry6.pack(pady=2, fill="x")

# Calculate button
ttk.Button(frame, text="Calculate f and Head Loss", command=process_darcy).pack(pady=15)

# Output
ttk.Label(frame, text="Output:").pack(pady=2, anchor="w")
output_var = tk.StringVar()
output = ttk.Entry(frame, textvariable=output_var, state="readonly", justify="center")
output.pack(pady=5, fill="x")

root.mainloop()
