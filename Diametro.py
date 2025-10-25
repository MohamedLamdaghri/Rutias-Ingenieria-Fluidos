import tkinter as tk
from tkinter import ttk, messagebox
import math

class DiameterIterationCalculator:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Iteración del Diámetro del Conducto")
        self.root.geometry("800x600")

        # Propiedades del fluido (agua aprox. 20°C)
        self.ro = 1000.0      # densidad [kg/m3]
        self.mu = 1e-3        # viscosidad dinámica [Pa·s]
        self.g  = 9.81        # gravedad [m/s2]

        # Variables de entrada
        self.roughness_var = tk.DoubleVar(value=0.0)  # k (ε) [m]
        self.length_var    = tk.DoubleVar()
        self.flow_var      = tk.DoubleVar()
        self.hf_var        = tk.DoubleVar()

        self.current_frame = None
        self.setup_main_frame()

    # ---------- UI ----------
    def setup_main_frame(self):
        if self.current_frame:
            self.current_frame.destroy()

        self.current_frame = ttk.Frame(self.root, padding="10")
        self.current_frame.pack(fill='both', expand=True)

        ttk.Label(self.current_frame, text="Cálculo Iterativo del Diámetro (D)", 
                  font=("Arial", 14, "bold")).grid(row=0, column=0, columnspan=2, pady=10)

        ttk.Label(self.current_frame, text="Rugosidad (ε) [m]:").grid(row=1, column=0, sticky='w', pady=5)
        ttk.Entry(self.current_frame, textvariable=self.roughness_var).grid(row=1, column=1, pady=5)

        ttk.Label(self.current_frame, text="Longitud (L) [m]:").grid(row=2, column=0, sticky='w', pady=5)
        ttk.Entry(self.current_frame, textvariable=self.length_var).grid(row=2, column=1, pady=5)

        ttk.Label(self.current_frame, text="Caudal (Q) [m³/s]:").grid(row=3, column=0, sticky='w', pady=5)
        ttk.Entry(self.current_frame, textvariable=self.flow_var).grid(row=3, column=1, pady=5)

        ttk.Label(self.current_frame, text="Pérdida de carga (hf) [m]:").grid(row=4, column=0, sticky='w', pady=5)
        ttk.Entry(self.current_frame, textvariable=self.hf_var).grid(row=4, column=1, pady=5)

        ttk.Button(self.current_frame, text="Iniciar Iteración", 
                   command=self.iterate_diameter).grid(row=6, column=0, columnspan=2, pady=15)

    # ---------- Núcleo de cálculo ----------
    @staticmethod
    def _log10(x):
        return math.log(x, 10.0)

    def D0_inicial(self, Q, L, Hf, k, nu):
        """
        D0 = 0.66 * [ k^1.25 * (Q^2 L / (g Hf))^4.75 + nu * Q^9.4 * (L/(g Hf))^5.2 ]^0.04
        """
        if Hf <= 0 or L <= 0 or Q <= 0:
            raise ValueError("Q, L y hf deben ser positivos.")
        term1 = (k**1.25) * ((Q*Q * L) / (self.g * Hf))**4.75
        term2 = (nu) * (Q**9.4) * ((L) / (self.g * Hf))**5.2
        return 0.66 * (term1 + term2)**0.04

    def iterate_diameter(self):
        try:
            k  = float(self.roughness_var.get())
            L  = float(self.length_var.get())
            Q  = float(self.flow_var.get())
            hf = float(self.hf_var.get())

            if any(v <= 0 for v in [L, Q, hf]):
                messagebox.showerror("Error", "L, Q y hf deben ser mayores que cero.")
                return
            if k < 0:
                messagebox.showerror("Error", "La rugosidad no puede ser negativa.")
                return

            nu = self.mu / self.ro  # viscosidad cinemática [m2/s]

            # D inicial desde la fórmula solicitada
            D = self.D0_inicial(Q, L, hf, k, nu)

            error = 100.0
            iteration = 0
            max_iterations = 200
            Re = 0.0
            f = 0.02
            regimen = ""

            iteration_data = []

            while error > 0.1 and iteration < max_iterations:
                # Velocidad y Reynolds con el D actual
                A = math.pi * D**2 / 4.0
                V = Q / A
                Re = abs(self.ro * V * D / self.mu)

                # Factor de fricción según reglas
                f, regimen = self.calculate_f_with_rules(Re, D, k)

                # Actualización de D (Darcy–Weisbach)
                D_new = ((8.0 * f * L * Q**2) / (math.pi**2 * self.g * hf)) ** 0.2

                error = abs((D_new - D) / D_new) * 100.0

                iteration_data.append({
                    'iter': iteration + 1,
                    'D': D,
                    'Re': Re,
                    'f': f,
                    'error': error,
                    'V': V,
                    'regimen': regimen
                })

                D = D_new
                iteration += 1

            # último registro
            iteration_data.append({
                'iter': iteration,
                'D': D,
                'Re': Re,
                'f': f,
                'error': error,
                'V': Q / (math.pi * D**2 / 4.0),
                'regimen': regimen
            })

            self.show_result(iteration_data, D, Re, f, error, iteration, regimen)

        except Exception as e:
            messagebox.showerror("Error", f"Error en el cálculo: {str(e)}")

    # ---- Selección de ecuación de fricción según tus reglas ----
    def calculate_f_with_rules(self, Re, D, k):
        """Devuelve (f, regimen_str)."""
        if Re <= 2300.0:
            return 64.0 / Re, "laminar"

        if k == 0.0:
            # 1ª von Kármán–Prandtl (tubo liso)
            return self.f_von_karman_suave(Re), "turbulento (liso, VK1)"

        # Umbral de transición: 3500/(k/D) = 3500*D/k
        threshold = 3500.0 * D / k

        if Re < threshold:
            return self.f_colebrook(Re, k, D), "turbulento (transición, Colebrook)"
        else:
            return self.f_von_karman_rugoso(D, k), "turbulento (rugoso, VK2)"

    # ---- Ecuaciones de fricción ----
    def _newton_solve(self, h, dh, x0, tol=1e-12, max_iter=50):
        x = x0
        for _ in range(max_iter):
            fx = h(x)
            dfx = dh(x)
            if dfx == 0 or not math.isfinite(dfx):
                break
            x_new = x - fx/dfx
            if not math.isfinite(x_new) or x_new <= 0:
                x_new = 0.5*x
            if abs(x_new - x) < tol*max(1.0, x_new):
                return max(x_new, 1e-12)
            x = x_new
        return max(x, 1e-12)

    def f_von_karman_suave(self, Re):
        # 1/sqrt(f) = 2*log10(Re*sqrt(f)) - 0.8
        LOG10 = math.log(10.0)
        def h(x):  # x = sqrt(f)
            return 1.0/x - 2.0*math.log(Re*x)/LOG10 + 0.8
        def dh(x):
            return -1.0/(x*x) - 2.0/(LOG10*x)
        x0 = 1.0 / math.sqrt(0.02)
        x  = self._newton_solve(h, dh, x0)
        return min(1.0, max(1e-6, x*x))

    def f_colebrook(self, Re, k, D):
        # 1/sqrt(f) = -2*log10( k/(3.71D) + 2.51/(Re*sqrt(f)) )
        LOG10 = math.log(10.0)
        A = (k/(3.71*D)) if D > 0 else 0.0
        B = 2.51/ Re
        def h(x):   # x = sqrt(f)
            return 1.0/x + 2.0*math.log(A + B/x)/LOG10
        def dh(x):
            denom = A + B/x
            return -(1.0/(x*x)) + 2.0*((-B/(x*x)) / (denom*LOG10))
        x0 = 1.0 / math.sqrt(0.02)
        x  = self._newton_solve(h, dh, x0)
        return min(1.0, max(1e-6, x*x))

    def f_von_karman_rugoso(self, D, k):
        # 1/sqrt(f) = 2*log10(D/k) + 1.14  -> f explícita
        denom = 2.0*self._log10(D/k) + 1.14
        return min(1.0, max(1e-6, 1.0/(denom*denom)))

    # ---------- Resultados ----------
    def show_result(self, iteration_data, final_D, final_Re, final_f, final_error, iterations, regimen_final):
        result_window = tk.Toplevel(self.root)
        result_window.title("Resultado de la Iteración")
        result_window.geometry("980x620")

        main_frame = ttk.Frame(result_window, padding="10")
        main_frame.pack(fill='both', expand=True)

        result_text = (
            f"Diámetro final: {final_D:.5f} m\n"
            f"Número de Reynolds: {final_Re:.2f}\n"
            f"Factor f final: {final_f:.6f}\n"
            f"Régimen/ecuación: {regimen_final}\n"
            f"Error final: {final_error:.4f}%\n"
            f"Iteraciones: {iterations}"
        )

        ttk.Label(main_frame, text="Resultado Final", font=("Arial", 12, "bold")).pack(pady=5)
        ttk.Label(main_frame, text=result_text, font=("Arial", 10)).pack(pady=5)

        ttk.Separator(main_frame, orient='horizontal').pack(fill='x', pady=10)

        ttk.Label(main_frame, text="Tabla de Iteraciones", font=("Arial", 12, "bold")).pack(pady=5)

        columns = ('Iter', 'D (m)', 'Re', 'f', 'Error (%)', 'V (m/s)', 'Régimen')
        tree = ttk.Treeview(main_frame, columns=columns, show='headings', height=12)

        tree.heading('Iter', text='Iteración')
        tree.heading('D (m)', text='Diámetro (m)')
        tree.heading('Re', text='Reynolds')
        tree.heading('f', text='Factor f')
        tree.heading('Error (%)', text='Error (%)')
        tree.heading('V (m/s)', text='Velocidad (m/s)')
        tree.heading('Régimen', text='Régimen/Ecuación')

        tree.column('Iter', width=70)
        tree.column('D (m)', width=120)
        tree.column('Re', width=120)
        tree.column('f', width=100)
        tree.column('Error (%)', width=100)
        tree.column('V (m/s)', width=120)
        tree.column('Régimen', width=240)

        for data in iteration_data:
            tree.insert('', 'end', values=(
                data['iter'],
                f"{data['D']:.6f}",
                f"{data['Re']:.2f}",
                f"{data['f']:.6f}",
                f"{data['error']:.4f}",
                f"{data['V']:.4f}",
                data['regimen']
            ))

        scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscroll=scrollbar.set)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        ttk.Button(main_frame, text="Cerrar", command=result_window.destroy).pack(pady=10)

    def run(self):
        self.root.mainloop()


if __name__ == "__main__":
    app = DiameterIterationCalculator()
    app.run()
