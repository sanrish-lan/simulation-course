import customtkinter as ctk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from numba import njit
import time


@njit(fastmath=True)
def solve_heat_equation(L, T_end, dx, dt, lam, rho, c, T_left, T_right, T_init):
    nx = int(L / dx) + 1
    nt = int(T_end / dt)

    x = np.linspace(0, L, nx)
    T = np.full(nx, T_init, dtype=np.float64)

    T[0] = T_left
    T[-1] = T_right

    h2 = dx * dx
    coeff_AC = lam / h2
    coeff_B = (2 * lam / h2) + (rho * c / dt)
    factor_F = (rho * c / dt)

    alpha = np.zeros(nx)
    beta = np.zeros(nx)

    for _ in range(nt):
        T_prev = T.copy()
        alpha[0] = 0.0
        beta[0] = T_left

        for i in range(1, nx - 1):
            F_i = -factor_F * T_prev[i]

            denom = coeff_B - coeff_AC * alpha[i - 1]
            alpha[i] = coeff_AC / denom
            beta[i] = (coeff_AC * beta[i - 1] - F_i) / denom

        T[-1] = T_right
        for i in range(nx - 2, -1, -1):
            T[i] = alpha[i] * T[i + 1] + beta[i]

    return x, T


ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")


class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Моделирование теплопроводности (Серебро)")
        self.geometry("1100x700")

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.sidebar_frame = ctk.CTkFrame(self, width=250, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")

        self.logo_label = ctk.CTkLabel(self.sidebar_frame, text="Параметры", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))

        self.entries = {}

        self.create_group("Свойства материала (Ag)", [
            ("L", "Длина (м)", "0.1"),
            ("lam", "Теплопроводность λ", "429.0"),
            ("rho", "Плотность ρ", "10490.0"),
            ("c", "Теплоемкость c", "235.0")
        ], start_row=1)

        self.create_group("Граничные условия", [
            ("T_left", "T слева (°C)", "300.0"),
            ("T_right", "T справа (°C)", "90.0"),
            ("T_init", "T начальная (°C)", "20.0")
        ], start_row=6)

        self.create_group("Параметры расчета", [
            ("t_end", "Время (с)", "2.0"),
            ("dx", "Шаг dx (м)", "0.1"),
            ("dt", "Шаг dt (с)", "0.1")
        ], start_row=10)

        self.calc_button = ctk.CTkButton(self.sidebar_frame, text="РАССЧИТАТЬ", command=self.run_simulation)
        self.calc_button.grid(row=14, column=0, padx=20, pady=20, sticky="ew")

        self.result_label = ctk.CTkLabel(self.sidebar_frame, text="Ожидание...", text_color="gray")
        self.result_label.grid(row=15, column=0, padx=20, pady=10)

        self.plot_frame = ctk.CTkFrame(self, corner_radius=10)
        self.plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

        plt.style.use('dark_background')
        self.fig, self.ax = plt.subplots()
        self.fig.set_facecolor('#2b2b2b')
        self.ax.set_facecolor('#1a1a1a')

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=5, pady=5)

        toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        toolbar.update()
        toolbar.pack(side="bottom", fill="x")

    def create_group(self, title, items, start_row):
        lbl = ctk.CTkLabel(self.sidebar_frame, text=title, anchor="w", text_color="#3B8ED0")
        lbl.grid(row=start_row, column=0, padx=20, pady=(10, 0), sticky="w")

        current_row = start_row + 1
        for key, label_text, default_val in items:
            frame = ctk.CTkFrame(self.sidebar_frame, fg_color="transparent")
            frame.grid(row=current_row, column=0, padx=20, pady=2, sticky="ew")

            lbl_w = ctk.CTkLabel(frame, text=label_text, width=120, anchor="w")
            lbl_w.pack(side="left")

            entry = ctk.CTkEntry(frame, width=80)
            entry.insert(0, default_val)
            entry.pack(side="right")

            self.entries[key] = entry
            current_row += 1

    def get_float(self, key):
        return float(self.entries[key].get())

    def run_simulation(self):
        try:
            L = self.get_float("L")
            lam = self.get_float("lam")
            rho = self.get_float("rho")
            c = self.get_float("c")

            T_left = self.get_float("T_left")
            T_right = self.get_float("T_right")
            T_init = self.get_float("T_init")

            t_end = self.get_float("t_end")
            dx = self.get_float("dx")
            dt = self.get_float("dt")

            if dx <= 0 or dt <= 0 or L <= 0:
                self.result_label.configure(text="Ошибка: Шаги и длина > 0", text_color="red")
                return

            start_time = time.time()

            x, T = solve_heat_equation(L, t_end, dx, dt, lam, rho, c, T_left, T_right, T_init)

            exec_time = time.time() - start_time

            self.ax.clear()
            self.ax.plot(x, T, color='#1f6aa5', linewidth=2, label=f't = {t_end} с')
            self.ax.axhline(y=T_init, color='gray', linestyle='--', alpha=0.5, label='Нач. темп.')

            self.ax.set_title(f"Распределение температуры (Сетка: {len(x)} узлов)", color="white")
            self.ax.set_xlabel("Координата x (м)", color="white")
            self.ax.set_ylabel("Температура T (°C)", color="white")
            self.ax.grid(True, color='#444444', linestyle='--')
            self.ax.legend(facecolor='#333333', labelcolor='white')

            self.ax.tick_params(axis='x', colors='white')
            self.ax.tick_params(axis='y', colors='white')
            for spine in self.ax.spines.values():
                spine.set_color('#555555')

            self.canvas.draw()

            center_idx = len(T) // 2
            res_text = f"Расчет выполнен за {exec_time:.4f} с\nT в центре: {T[center_idx]:.2f} °C"
            self.result_label.configure(text=res_text, text_color="green")

        except ValueError:
            self.result_label.configure(text="Ошибка: Введите числа", text_color="red")
        except Exception as e:
            self.result_label.configure(text=f"Ошибка: {str(e)}", text_color="red")


if __name__ == "__main__":
    app = App()
    app.mainloop()