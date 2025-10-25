import math
import tkinter as tk
from tkinter import messagebox

g = 9.81
PI = math.pi

def area(D):
    return PI * (D ** 2) / 4


def reynolds(Q, D, v):
    return (Q / area(D)) * D / v


def colebrook(Re, k, D, f0=0.02, tol=1e-8):
    if Re < 2000:
        return 64 / Re, 0  # flujo laminar
    
    f = f0
    # enough big to enter the loop
    er = 100
    while(er > tol):
        f_new = (1 / (-2 * math.log10((k / (3.7 * D)) + (2.51 / (Re * math.sqrt(f)))))) ** 2
        er = abs((f_new - f) / f_new)

        # update f
        f = f_new
    return f

def swameeJainF(Re, k, D):
    if Re < 2000:
        return 64 / Re
    return 0.25 / (math.log10((k / (3.7 * D)) + (5.74 / (Re ** 0.9)))) ** 2

def getQfromHf(f, D, L, Hf):
    num = (PI ** 2) * g * Hf * (D ** 5)
    den = 8 * f * L
    return math.sqrt(abs(num) / den)

def initialGuessQ(D , L, k, v, Hf):
    aFactor = -0.965*math.sqrt(g*(D ** 5)*Hf/(L))
    bFactor = math.log((k/(3.7*D)) + math.sqrt((3.17*(v**2)*L)/(g*(D**3)*Hf)))
    return aFactor * bFactor

def solveQ(D, L, k, v, Hf, tol=1e-6):
    Q = initialGuessQ(D , L, k, v, Hf)
    Re_est = reynolds(Q, D, v)
    f = swameeJainF(Re_est, k, D)

    # enough big error value to enter the loop
    err = 100
    while(err > tol):
        V = Q / area(D)
        Re = reynolds(Q, D, v)

        f = colebrook(Re, k, D, f)
        Q_new = getQfromHf(f, D, L, Hf)
        err = abs(Q_new - Q) / Q_new

        #Update Q
        Q = Q_new

    V = Q / area(D)
    Re = reynolds(Q, D, v)
    
    if Re < 2300:
        regimen = "Laminar"
    else:
        regimen = "Turbulent"
    
    return Q, Re, f, regimen

def process():
    try:
        D = float(entry_D.get())
        L = float(entry_L.get())
        k = float(entry_k.get())
        v = float(entry_nu.get())
        Hf = float(entry_Hf.get())
        er = float(entry_Er.get())

        Q, Re, f, regime = solveQ(D, L, k, v, abs(Hf), (er/100))

        output_var.set(f"Q (volumetric flow m3/s): {Q:.6f}     \nf (friction factor): {f:.6f}     \nRe (reynolds)  {Re:.4f}    \nregime: {regime}")

    except ValueError as e:
        messagebox.showerror("Input error", str(e))
    except Exception as e:
        messagebox.showerror("Unexpected error", str(e))


root = tk.Tk()
root.title("Calculo de Caudal : Mohamed Lamdaghri")
root.geometry("800x400")
root.resizable(True, True)

tk.Label(root, text="Input data (SI)", font=("Arial", 12, "bold")).pack(pady=5)

frame_inputs = tk.Frame(root)
frame_inputs.pack(pady=5)

labels = ["D [m]", "L [m]", "k [m]", "ν [m²/s]", "Hf [m]", "Er [%]"]
defaults = ["0.03", "1000", "0.001", "0.000001", "-41.743", "2"]
entries = []

for i, (lab, val) in enumerate(zip(labels, defaults)):
    tk.Label(frame_inputs, text=lab, anchor="w", width=25).grid(row=i, column=0, pady=3, sticky="w")
    e = tk.Entry(frame_inputs, width=15)
    e.insert(0, val)
    e.grid(row=i, column=1, pady=3)
    entries.append(e)

entry_D, entry_L, entry_k, entry_nu, entry_Hf, entry_Er= entries

tk.Button(root, text="Get Q", command=process, bg="#4CAF50", fg="white", font=("Arial", 12, "bold"), width=20).pack(pady=10)


# Output
tk.Label(root, text="Output:").pack(pady=2, anchor="w")
output_var = tk.StringVar()
output = tk.Entry(root, textvariable=output_var, state="readonly", justify="center")
output.pack(pady=5, fill="x")

root.mainloop()