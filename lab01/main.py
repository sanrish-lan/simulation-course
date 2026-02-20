import customtkinter as ctk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from numba import njit




@njit
def physics_step(x, y, vx, vy, dt, k, m, g):
    v_abs = np.sqrt(vx ** 2 + vy ** 2) # Считаем полную скорость
    ax = - (k / m) * v_abs * vx # Находим мгновенное ускорение X
    ay = - g - (k / m) * v_abs * vy # Находим мгновенное ускорение Y

    vx_new = vx + ax * dt # Обновляем скорость
    vy_new = vy + ay * dt
    x_new = x + vx_new * dt # Находим новую координату
    y_new = y + vy_new * dt
    return x_new, y_new, vx_new, vy_new


@njit
def solve_full_flight(v0, angle, m, rho, Cd, S, dt):
    g = 9.81
    k = 0.5 * Cd * rho * S # Считаем коэффициент сопротивления воздуха
    rad = np.radians(angle) # Переводим угол из градусов в радианы для математики
    vx, vy = v0 * np.cos(rad), v0 * np.sin(rad) # Разлагаем скорость на горизонтальную и вертикальную
    x, y = 0.0, 0.0


    path_x = [0.0]
    path_y = [0.0]
    max_h = 0.0

    while y >= 0: # Пока тело не коснулось земли (высота не стала отрицательной)
        x, y, vx, vy = physics_step(x, y, vx, vy, dt, k, m, g) # Считаем следующую точку
        if y >= 0:
            path_x.append(x)
            path_y.append(y)
            if y > max_h: max_h = y # Если мы сейчас выше, чем были раньше, обновляем рекорд высоты
        if x > 10000: break

    return np.array(path_x), np.array(path_y), x, max_h, vx, vy


@njit
def solve_batch(x, y, vx, vy, dt, k, m, g, num_steps):
    batch_x = []
    batch_y = []
    max_h_local = y

    for _ in range(num_steps):
        x, y, vx, vy = physics_step(x, y, vx, vy, dt, k, m, g)
        if y < 0: break
        batch_x.append(x)
        batch_y.append(y)
        if y > max_h_local: max_h_local = y

    return x, y, vx, vy, batch_x, batch_y, max_h_local


ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Flight Sim + Numba JIT Acceleration")
        self.geometry("1280x850")

        self.run_index = 0
        self.colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6']
        self.is_running = False

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)


        self.sidebar = ctk.CTkFrame(self, width=320)
        self.sidebar.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        ctk.CTkLabel(self.sidebar, text="Параметры", font=("Arial", 20, "bold")).pack(pady=15)

        self.inputs = {}
        fields = [
            ("Нач. скорость (м/с)", "v0", "100.0"),
            ("Угол запуска (°)", "angle", "45.0"),
            ("Масса тела (кг)", "mass", "1.0"),
            ("Плотность воздуха", "rho", "1.29"),
            ("Cd (сопротивление)", "cd", "0.15"),
            ("Площадь (м²)", "area", "0.01"),
            ("Шаг dt (с)", "dt", "0.001")
        ]

        for label_text, key, default in fields:
            f = ctk.CTkFrame(self.sidebar, fg_color="transparent")
            f.pack(fill="x", padx=20, pady=2)
            ctk.CTkLabel(f, text=label_text).pack(side="left")
            entry = ctk.CTkEntry(f, width=80)
            entry.insert(0, default)
            entry.pack(side="right")
            self.inputs[key] = entry

        ctk.CTkLabel(self.sidebar, text="Настройки анимации", font=("Arial", 14, "bold")).pack(pady=(20, 5))

        self.anim_var = ctk.BooleanVar(value=True)
        self.cb_anim = ctk.CTkCheckBox(self.sidebar, text="Включить анимацию", variable=self.anim_var)
        self.cb_anim.pack(pady=5, padx=20, anchor="w")

        ctk.CTkLabel(self.sidebar, text="Скорость отрисовки:").pack(padx=20, anchor="w")
        self.speed_slider = ctk.CTkSlider(self.sidebar, from_=1, to=500, number_of_steps=100)
        self.speed_slider.set(50)
        self.speed_slider.pack(pady=5, padx=20, fill="x")

        self.btn_run = ctk.CTkButton(self.sidebar, text="ЗАПУСТИТЬ", font=("Arial", 14, "bold"),
                                     command=self.start_simulation)
        self.btn_run.pack(pady=(25, 10), padx=20, fill="x")

        self.btn_clear = ctk.CTkButton(self.sidebar, text="Очистить графики", fg_color="#c0392b",
                                       command=self.clear_all)
        self.btn_clear.pack(pady=10, padx=20, fill="x")

        self.right_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.right_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        self.right_frame.rowconfigure(0, weight=3)
        self.right_frame.rowconfigure(1, weight=1)
        self.right_frame.columnconfigure(0, weight=1)

        self.plot_frame = ctk.CTkFrame(self.right_frame)
        self.plot_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 10))
        self.setup_plot()

        self.table_frame = ctk.CTkFrame(self.right_frame)
        self.table_frame.grid(row=1, column=0, sticky="nsew")
        self.setup_table()

    def setup_plot(self):
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.fig.patch.set_facecolor('#2b2b2b')
        self.ax.set_facecolor('#1e1e1e')
        self.ax.tick_params(colors='white')
        self.ax.grid(True, color='#333333', linestyle='--')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=5, pady=5)

    def setup_table(self):
        headers = ["Шаг (с)", "Дальность (м)", "Высота (м)", "Скорость (м/с)"]
        h_frame = ctk.CTkFrame(self.table_frame, fg_color="#34495e", height=35)
        h_frame.pack(fill="x", padx=2, pady=2)
        for i, h in enumerate(headers):
            lbl = ctk.CTkLabel(h_frame, text=h, font=("Arial", 12, "bold"))
            lbl.place(relx=i * 0.25, rely=0.5, anchor="w", x=15)
        self.scroll_table = ctk.CTkScrollableFrame(self.table_frame, fg_color="transparent")
        self.scroll_table.pack(fill="both", expand=True)

    def start_simulation(self):
        if self.is_running: return
        try:
            v0 = float(self.inputs['v0'].get())
            angle = float(self.inputs['angle'].get())
            m = float(self.inputs['mass'].get())
            rho = float(self.inputs['rho'].get())
            Cd = float(self.inputs['cd'].get())
            S = float(self.inputs['area'].get())
            self.dt = float(self.inputs['dt'].get())
        except:
            return

        self.is_running = True
        self.btn_run.configure(state="disabled")
        self.k_val = 0.5 * Cd * rho * S
        self.m_val = m

        color = self.colors[self.run_index % len(self.colors)]
        self.line, = self.ax.plot([], [], color=color, lw=2, label=f"dt={self.dt}")
        self.point, = self.ax.plot([], [], color=color, marker='o', ms=6)
        self.ax.legend(facecolor='#2b2b2b', labelcolor='white')

        if self.anim_var.get():
            rad = np.radians(angle)
            self.st = {'x': 0.0, 'y': 0.0, 'vx': v0 * np.cos(rad), 'vy': v0 * np.sin(rad)}
            self.full_x, self.full_y = [0.0], [0.0]
            self.max_h = 0.0
            self.animate_loop()
        else:
            px, py, dist, mh, vx_f, vy_f = solve_full_flight(v0, angle, m, rho, Cd, S, self.dt)
            self.line.set_data(px, py)
            self.ax.relim();
            self.ax.autoscale_view()
            self.update_table(self.dt, dist, mh, np.sqrt(vx_f ** 2 + vy_f ** 2))
            self.finish_sim()

    def animate_loop(self):
        batch_size = int(self.speed_slider.get())
        nx, ny, nvx, nvy, bx, by, mh_local = solve_batch(
            self.st['x'], self.st['y'], self.st['vx'], self.st['vy'],
            self.dt, self.k_val, self.m_val, 9.81, batch_size
        )

        self.st.update({'x': nx, 'y': ny, 'vx': nvx, 'vy': nvy})
        self.full_x.extend(bx)
        self.full_y.extend(by)
        if mh_local > self.max_h: self.max_h = mh_local

        self.line.set_data(self.full_x, self.full_y)
        self.point.set_data([nx], [ny])
        self.ax.relim();
        self.ax.autoscale_view()
        self.canvas.draw_idle()

        if ny < 0:
            v_fin = np.sqrt(nvx ** 2 + nvy ** 2)
            self.update_table(self.dt, nx, self.max_h, v_fin)
            self.finish_sim()
        else:
            self.after(10, self.animate_loop)

    def update_table(self, dt, x, h, v):
        row = ctk.CTkFrame(self.scroll_table, fg_color="transparent")
        row.pack(fill="x", pady=1)
        vals = [f"{dt}", f"{x:.2f}", f"{h:.2f}", f"{v:.2f}"]
        for i, val in enumerate(vals):
            ctk.CTkLabel(row, text=val).place(relx=i * 0.25, rely=0.5, anchor="w", x=15)

    def finish_sim(self):
        self.is_running = False
        self.btn_run.configure(state="normal")
        self.point.set_visible(False)
        self.canvas.draw()
        self.run_index += 1

    def clear_all(self):
        if self.is_running: return
        self.ax.clear()
        self.setup_plot()
        for child in self.scroll_table.winfo_children(): child.destroy()
        self.run_index = 0


if __name__ == "__main__":
    App().mainloop()
