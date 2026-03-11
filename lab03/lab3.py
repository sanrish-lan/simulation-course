import customtkinter as ctk
import numpy as np
from PIL import Image, ImageTk

# Константы состояний клеток
EMPTY = 0
TREE = 1
BURNING = 2
ASH = 3
OBSTACLE = 4

# Цвета для отрисовки (RGB)
COLORS = np.array([
    [40, 30, 20],  # 0: Земля (Темно-коричневая)
    [34, 139, 34],  # 1: Дерево (Зеленый)
    [255, 255, 0],  # 2: Огонь (ЯРКО-ЖЕЛТЫЙ)
    [90, 90, 90],  # 3: Пепел (Темно-серый)
    [30, 144, 255]  # 4: Вода (Синий)
], dtype=np.uint8)


class UltimateForestFireApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Лесные пожары: Процедурная генерация и Роза ветров")
        self.geometry("1100x750")

        self.grid_size = 150
        self.is_running = False

        # Направление ветра (векторы)
        self.wind_dx = 0
        self.wind_dy = 0
        self.wind_buttons = {}

        self.grid = np.zeros((self.grid_size, self.grid_size), dtype=np.int8)
        self.burn_time = np.zeros((self.grid_size, self.grid_size), dtype=np.int8)

        # Переменные для инструментов рисования
        self.draw_mode = ctk.IntVar(value=BURNING)
        self.brush_size_var = ctk.IntVar(value=0)  # 0: 1x1, 1: 3x3, 2: 5x5

        self.setup_ui()
        self.generate_natural_terrain()
        self.update_canvas()

    def setup_ui(self):
        self.sidebar = ctk.CTkFrame(self, width=320, corner_radius=0)
        self.sidebar.pack(side="left", fill="y", padx=10, pady=10)

        ctk.CTkLabel(self.sidebar, text="Управление", font=ctk.CTkFont(size=18, weight="bold")).pack(pady=(10, 5))

        self.btn_start = ctk.CTkButton(self.sidebar, text="Старт", command=self.toggle_simulation)
        self.btn_start.pack(pady=5, fill="x", padx=10)

        self.btn_reset = ctk.CTkButton(self.sidebar, text="Сгенерировать новый мир",
                                       command=self.generate_natural_terrain)
        self.btn_reset.pack(pady=5, fill="x", padx=10)

        # Роза Ветров
        ctk.CTkLabel(self.sidebar, text="Направление ветра", font=ctk.CTkFont(weight="bold")).pack(pady=(15, 5))
        wind_frame = ctk.CTkFrame(self.sidebar, fg_color="transparent")
        wind_frame.pack(pady=5)

        wind_dirs = [
            ("↖", -1, -1), ("↑", 0, -1), ("↗", 1, -1),
            ("←", -1, 0), ("•", 0, 0), ("→", 1, 0),
            ("↙", -1, 1), ("↓", 0, 1), ("↘", 1, 1)
        ]

        for i, (symbol, dx, dy) in enumerate(wind_dirs):
            btn = ctk.CTkButton(wind_frame, text=symbol, width=40, height=40,
                                font=ctk.CTkFont(size=18, weight="bold"),
                                fg_color="#1F6AA5" if (dx == 0 and dy == 0) else "#3B3B3B",
                                hover_color="#144870",
                                command=lambda x=dx, y=dy: self.set_wind(x, y))
            btn.grid(row=i // 3, column=i % 3, padx=2, pady=2)
            self.wind_buttons[(dx, dy)] = btn

        ctk.CTkLabel(self.sidebar, text="Сила ветра").pack(pady=(5, 0))
        self.slider_wind = ctk.CTkSlider(self.sidebar, from_=0.0, to=0.5, number_of_steps=50)
        self.slider_wind.set(0.25)
        self.slider_wind.pack(fill="x", padx=10)

        ctk.CTkLabel(self.sidebar, text="Скорость роста (p)").pack(pady=(15, 0))
        self.slider_p = ctk.CTkSlider(self.sidebar, from_=0.0, to=0.05, number_of_steps=100)
        self.slider_p.set(0.005)
        self.slider_p.pack(fill="x", padx=10)

        ctk.CTkLabel(self.sidebar, text="Молнии (f)").pack(pady=(5, 0))
        self.slider_f = ctk.CTkSlider(self.sidebar, from_=0.0, to=0.0005, number_of_steps=100)
        self.slider_f.set(0.00005)
        self.slider_f.pack(fill="x", padx=10)

        # Выбор материала кисти
        ctk.CTkLabel(self.sidebar, text="Кисть:", font=ctk.CTkFont(weight="bold")).pack(pady=(15, 5))
        ctk.CTkRadioButton(self.sidebar, text="Огонь (Желтый)", variable=self.draw_mode, value=BURNING).pack(anchor="w",
                                                                                                             padx=20,
                                                                                                             pady=2)
        ctk.CTkRadioButton(self.sidebar, text="Вода/Озеро (Синий)", variable=self.draw_mode, value=OBSTACLE).pack(
            anchor="w", padx=20, pady=2)
        ctk.CTkRadioButton(self.sidebar, text="Лес (Зеленый)", variable=self.draw_mode, value=TREE).pack(anchor="w",
                                                                                                         padx=20,
                                                                                                         pady=2)
        ctk.CTkRadioButton(self.sidebar, text="Ластик (Земля)", variable=self.draw_mode, value=EMPTY).pack(anchor="w",
                                                                                                           padx=20,
                                                                                                           pady=2)

        # Выбор размера кисти
        ctk.CTkLabel(self.sidebar, text="Размер кисти:", font=ctk.CTkFont(weight="bold")).pack(pady=(15, 5))
        brush_frame = ctk.CTkFrame(self.sidebar, fg_color="transparent")
        brush_frame.pack(fill="x", padx=20)

        ctk.CTkRadioButton(brush_frame, text="1x1", variable=self.brush_size_var, value=0).pack(side="left",
                                                                                                padx=(0, 10))
        ctk.CTkRadioButton(brush_frame, text="3x3", variable=self.brush_size_var, value=1).pack(side="left",
                                                                                                padx=(0, 10))
        ctk.CTkRadioButton(brush_frame, text="5x5", variable=self.brush_size_var, value=2).pack(side="left")

        # Настройка холста (Канвас)
        self.canvas_frame = ctk.CTkFrame(self)
        self.canvas_frame.pack(side="right", expand=True, fill="both", padx=10, pady=10)

        self.canvas = ctk.CTkCanvas(self.canvas_frame, bg="black", highlightthickness=0)
        self.canvas.pack(expand=True, fill="both")

        self.canvas.bind("<B1-Motion>", self.paint)
        self.canvas.bind("<Button-1>", self.paint)

    def set_wind(self, dx, dy):
        self.wind_dx = dx
        self.wind_dy = dy
        for (kx, ky), btn in self.wind_buttons.items():
            if kx == dx and ky == dy:
                btn.configure(fg_color="#D35B58")
            else:
                btn.configure(fg_color="#3B3B3B")
                if kx == 0 and ky == 0 and dx != 0:
                    btn.configure(fg_color="#2B2B2B")
        if dx == 0 and dy == 0:
            self.wind_buttons[(0, 0)].configure(fg_color="#1F6AA5")

    def generate_natural_terrain(self):
        terrain = np.random.choice([0, 1], size=(self.grid_size, self.grid_size), p=[0.55, 0.45])
        for _ in range(4):
            neighbors = sum(np.roll(np.roll(terrain, i, axis=0), j, axis=1)
                            for i in (-1, 0, 1) for j in (-1, 0, 1)) - terrain
            terrain = np.where(neighbors >= 5, 1, 0)

        self.grid.fill(EMPTY)
        self.grid[terrain == 1] = OBSTACLE

        land_mask = (terrain == 0)
        tree_noise = np.random.rand(self.grid_size, self.grid_size)
        self.grid[land_mask & (tree_noise < 0.6)] = TREE

        self.burn_time.fill(0)
        self.update_canvas()

    def paint(self, event):
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        col = int(event.x / (w / self.grid_size))
        row = int(event.y / (h / self.grid_size))

        # Получаем радиус из радиокнопок: 0=1x1, 1=3x3, 2=5x5
        brush_radius = self.brush_size_var.get()
        val = self.draw_mode.get()

        if 0 <= row < self.grid_size and 0 <= col < self.grid_size:
            # Если радиус 0, рисуем одну клетку (1x1)
            if brush_radius == 0:
                self.grid[row, col] = val
                if val == BURNING:
                    self.burn_time[row, col] = 5
                else:
                    self.burn_time[row, col] = 0
            else:
                # Рисуем квадрат размером (2*radius + 1)
                r_start = max(0, row - brush_radius)
                r_end = min(self.grid_size, row + brush_radius + 1)
                c_start = max(0, col - brush_radius)
                c_end = min(self.grid_size, col + brush_radius + 1)

                self.grid[r_start:r_end, c_start:c_end] = val

                if val == BURNING:
                    self.burn_time[r_start:r_end, c_start:c_end] = 5
                else:
                    self.burn_time[r_start:r_end, c_start:c_end] = 0

            self.update_canvas()

    def toggle_simulation(self):
        self.is_running = not self.is_running
        if self.is_running:
            self.btn_start.configure(text="Пауза", fg_color="red")
            self.run_step()
        else:
            self.btn_start.configure(text="Старт", fg_color=["#3B8ED0", "#1F6AA5"])

    def run_step(self):
        if not self.is_running: return
        self.calculate_next_generation()
        self.update_canvas()
        self.after(60, self.run_step)

    def calculate_next_generation(self):
        p_growth = self.slider_p.get()
        f_lightning = self.slider_f.get()
        wind_str = self.slider_wind.get()

        new_grid = self.grid.copy()
        new_burn_time = self.burn_time.copy()

        is_burning = (self.grid == BURNING)
        rand_mat = np.random.rand(self.grid_size, self.grid_size)

        new_burn_time[is_burning] -= 1
        just_finished = is_burning & (new_burn_time <= 0)
        new_grid[just_finished] = ASH

        ash_decay = (self.grid == ASH) & (rand_mat < 0.02)
        new_grid[ash_decay] = EMPTY

        is_burn_int = is_burning.astype(float)

        # Соседи (Окрестность фон Неймана)
        N = np.roll(is_burn_int, 1, axis=0)
        N[0, :] = 0
        S = np.roll(is_burn_int, -1, axis=0)
        S[-1, :] = 0
        W = np.roll(is_burn_int, 1, axis=1)
        W[:, 0] = 0
        E = np.roll(is_burn_int, -1, axis=1)
        E[:, -1] = 0

        # Искры (Огонь на расстоянии 2 клеток)
        N2 = np.roll(is_burn_int, 2, axis=0)
        N2[0:2, :] = 0
        S2 = np.roll(is_burn_int, -2, axis=0)
        S2[-2:, :] = 0
        W2 = np.roll(is_burn_int, 2, axis=1)
        W2[:, 0:2] = 0
        E2 = np.roll(is_burn_int, -2, axis=1)
        E2[:, -2:] = 0

        # Базовое давление огня
        fire_pressure = (N + S + W + E) * 0.15
        spark_pressure = (N2 + S2 + W2 + E2) * 0.01

        # Вектор ветра
        if self.wind_dx > 0:
            fire_pressure += W * wind_str
            spark_pressure += W2 * wind_str * 0.5
        elif self.wind_dx < 0:
            fire_pressure += E * wind_str
            spark_pressure += E2 * wind_str * 0.5

        if self.wind_dy > 0:
            fire_pressure += N * wind_str
            spark_pressure += N2 * wind_str * 0.5
        elif self.wind_dy < 0:
            fire_pressure += S * wind_str
            spark_pressure += S2 * wind_str * 0.5

        total_ignition_prob = fire_pressure + spark_pressure + f_lightning

        ignited = (self.grid == TREE) & (rand_mat < total_ignition_prob)
        new_grid[ignited] = BURNING
        new_burn_time[ignited] = np.random.randint(3, 7, size=np.count_nonzero(ignited))

        grow_empty = (self.grid == EMPTY) & (rand_mat < p_growth)
        grow_ash = (self.grid == ASH) & (rand_mat < p_growth * 5.0)
        new_grid[grow_empty | grow_ash] = TREE

        new_grid[self.grid == OBSTACLE] = OBSTACLE

        self.grid = new_grid
        self.burn_time = new_burn_time

    def update_canvas(self):
        img_array = COLORS[self.grid]
        image = Image.fromarray(img_array, 'RGB')

        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        if w <= 1 or h <= 1: w, h = 600, 600

        image = image.resize((w, h), Image.NEAREST)
        self.photo = ImageTk.PhotoImage(image)
        self.canvas.create_image(0, 0, image=self.photo, anchor="nw")


if __name__ == "__main__":
    ctk.set_appearance_mode("dark")
    ctk.set_default_color_theme("blue")

    app = UltimateForestFireApp()
    app.mainloop()