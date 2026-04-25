import sys
import math
import random
import numpy as np
from scipy.stats import norm, chi2

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QLabel, QLineEdit, QPushButton, QTextEdit, QMessageBox
)
from PyQt6.QtCore import Qt

import matplotlib

matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure

MODERN_STYLE = """
QMainWindow, QWidget {
    background-color: #1e1e2e;
    color: #cdd6f4;
    font-family: 'Segoe UI', Arial, sans-serif;
    font-size: 14px;
}
QGroupBox {
    border: 1px solid #45475a;
    border-radius: 8px;
    margin-top: 14px;
    padding-top: 15px;
    font-weight: bold;
    color: #89b4fa;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 12px;
    padding: 0 5px;
}
QLineEdit {
    background-color: #313244;
    border: 1px solid #45475a;
    border-radius: 5px;
    padding: 6px 10px;
    color: #cdd6f4;
}
QLineEdit:focus {
    border: 1px solid #89b4fa;
    background-color: #3b3c50;
}
QLineEdit[readOnly="true"] {
    background-color: #181825;
    color: #a6adc8;
    border: 1px dashed #45475a;
}
QPushButton {
    background-color: #89b4fa;
    color: #11111b;
    border: none;
    border-radius: 5px;
    padding: 8px 12px;
    font-weight: bold;
}
QPushButton:hover { background-color: #b4befe; }
QPushButton:pressed { background-color: #74c7ec; }

QPushButton#BtnRandom {
    background-color: #45475a;
    color: #cdd6f4;
}
QPushButton#BtnRandom:hover { background-color: #585b70; }

QTextEdit {
    background-color: #181825;
    border: 1px solid #45475a;
    border-radius: 6px;
    padding: 8px;
    font-family: 'Consolas', 'Courier New', monospace;
    font-size: 13px;
}
QTabWidget::pane {
    border: 1px solid #45475a;
    border-radius: 6px;
    background-color: #1e1e2e;
}
QTabBar::tab {
    background-color: #313244;
    color: #a6adc8;
    padding: 10px 20px;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    margin-right: 2px;
}
QTabBar::tab:selected {
    background-color: #1e1e2e;
    color: #89b4fa;
    border-bottom: 2px solid #89b4fa;
    font-weight: bold;
}
QTabBar::tab:hover:!selected {
    background-color: #45475a;
}
"""

class ChartCanvas(Canvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(6, 4), tight_layout=True)
        super().__init__(self.fig)
        self.setParent(parent)
        self.axes = self.fig.add_subplot(111)

        self.fig.patch.set_facecolor('#1e1e2e')
        self.axes.set_facecolor('#181825')
        self.axes.tick_params(colors='#cdd6f4')
        for spine in self.axes.spines.values():
            spine.set_edgecolor('#45475a')
        self.axes.yaxis.label.set_color('#cdd6f4')
        self.axes.xaxis.label.set_color('#cdd6f4')
        self.axes.title.set_color('#cdd6f4')

    def plot_bars(self, x_labels, th_probs, emp_probs):
        self.axes.clear()
        x_indexes = np.arange(len(x_labels))
        width = 0.35

        self.axes.bar(x_indexes - width / 2, th_probs, width, label='Теоретические', color='#89b4fa', alpha=0.9)
        self.axes.bar(x_indexes + width / 2, emp_probs, width, label='Эмпирические', color='#f9e2af', alpha=0.9)

        self.axes.set_xticks(x_indexes)
        self.axes.set_xticklabels(x_labels)
        self.axes.set_ylabel("Вероятность")
        self.axes.legend(facecolor='#313244', edgecolor='#45475a', labelcolor='#cdd6f4')
        self.axes.grid(axis='y', color='#45475a', linestyle='--', alpha=0.6)
        self.draw()

    def plot_histogram(self, data, mean, var):
        self.axes.clear()

        if not data:
            self.draw()
            return

        std_dev = math.sqrt(var)
        n_size = len(data)

        if n_size <= 10:
            bins_count = n_size
        else:
            bins_count = min(100, int(math.sqrt(n_size)))

        min_val, max_val = min(data), max(data)
        if min_val == max_val:
            min_val -= (std_dev if std_dev > 0 else 1)
            max_val += (std_dev if std_dev > 0 else 1)

        self.axes.hist(data, bins=bins_count, density=True,
                       color='#89b4fa', edgecolor='#181825', alpha=0.7, label='Гистограмма')

        x_axis = np.linspace(min_val - std_dev, max_val + std_dev, 300)
        y_axis = norm.pdf(x_axis, mean, std_dev)
        self.axes.plot(x_axis, y_axis, linewidth=2.5, color='#a6e3a1', label='Теор. плотность')

        self.axes.set_ylabel("Плотность")
        self.axes.legend(facecolor='#313244', edgecolor='#45475a', labelcolor='#cdd6f4')
        self.axes.grid(axis='y', color='#45475a', linestyle='--', alpha=0.4)

        self.axes.set_xlim(min_val - std_dev * 0.5, max_val + std_dev * 0.5)

        self.draw()

class DiscreteTab(QWidget):
    def __init__(self):
        super().__init__()
        self.build_ui()

    def build_ui(self):
        layout = QHBoxLayout(self)

        panel_left = QGroupBox("Настройки распределения")
        form = QFormLayout(panel_left)
        panel_left.setFixedWidth(300)

        self.inputs_prob = []
        default_probs = [0.1, 0.2, 0.3, 0.3]
        for i in range(4):
            field = QLineEdit(str(default_probs[i]))
            field.textChanged.connect(self.recalc_auto_prob)
            self.inputs_prob.append(field)
            form.addRow(f"Prob {i + 1}:", field)

        self.field_auto = QLineEdit()
        self.field_auto.setReadOnly(True)
        form.addRow("Prob 5 (auto):", self.field_auto)

        self.btn_random = QPushButton("🎲 Случайные вероятности")
        self.btn_random.setObjectName("BtnRandom")
        self.btn_random.clicked.connect(self.randomize_probs)
        form.addRow(self.btn_random)

        form.addRow(QLabel(""), QLabel(""))

        self.field_n = QLineEdit("1000")
        form.addRow("Кол-во экспериментов:", self.field_n)

        self.btn_run = QPushButton("Запустить моделирование")
        self.btn_run.setMinimumHeight(45)
        self.btn_run.clicked.connect(self.run_simulation)
        form.addRow(self.btn_run)

        layout.addWidget(panel_left)

        panel_right = QVBoxLayout()
        self.canvas = ChartCanvas()
        panel_right.addWidget(self.canvas, stretch=2)

        self.lbl_result = QLabel("Запустите симуляцию для получения результатов...")
        self.lbl_result.setStyleSheet("font-size: 15px; line-height: 1.5;")
        self.lbl_result.setTextFormat(Qt.TextFormat.RichText)
        panel_right.addWidget(self.lbl_result)

        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        panel_right.addWidget(self.txt_log, stretch=1)

        layout.addLayout(panel_right)
        self.recalc_auto_prob()

    def randomize_probs(self):
        raw_vals = [random.uniform(0.1, 1.0) for _ in range(5)]
        total = sum(raw_vals)
        normalized = [round(v / total, 3) for v in raw_vals]
        for i in range(4):
            self.inputs_prob[i].setText(str(normalized[i]))

    def recalc_auto_prob(self):
        try:
            curr_sum = sum(float(f.text()) for f in self.inputs_prob if f.text())
            p5 = round(1.0 - curr_sum, 4)
            self.field_auto.setText(str(p5))
            if p5 < 0:
                self.field_auto.setStyleSheet("background-color: #4a1d23; border: 1px solid #f38ba8; color: #f38ba8;")
            else:
                self.field_auto.setStyleSheet("")
        except ValueError:
            pass

    def run_simulation(self):
        try:
            p_list = [float(f.text()) for f in self.inputs_prob]
            p_list.append(float(self.field_auto.text()))
            user_n = int(self.field_n.text())
        except ValueError:
            QMessageBox.warning(self, "Ошибка", "Проверьте корректность числовых данных.")
            return

        if any(p < 0 for p in p_list) or abs(sum(p_list) - 1.0) > 1e-4:
            QMessageBox.warning(self, "Ошибка", "Вероятности должны быть ≥ 0, а их сумма = 1.0")
            return

        states = [1, 2, 3, 4, 5]
        th_mean = sum(s * p for s, p in zip(states, p_list))
        th_var = sum(p * ((s - th_mean) ** 2) for s, p in zip(states, p_list))

        def simulate_logic(n_size):
            cum_probs = np.cumsum(p_list)
            sim_vals = []
            for _ in range(n_size):
                r = random.random()
                for idx, cp in enumerate(cum_probs):
                    if r <= cp:
                        sim_vals.append(states[idx])
                        break

            counts = {s: sim_vals.count(s) for s in states}
            emp_probs = [counts[s] / n_size for s in states]

            emp_mean = sum(s * p for s, p in zip(states, emp_probs))
            emp_var = sum(p * ((s - emp_mean) ** 2) for s, p in zip(states, emp_probs))

            err_mean = abs(emp_mean - th_mean) / th_mean * 100 if th_mean else 0
            err_var = abs(emp_var - th_var) / th_var * 100 if th_var else 0

            chi2_val = 0
            for s, p in zip(states, p_list):
                expected_f = n_size * p
                if expected_f > 0:
                    chi2_val += ((counts[s] - expected_f) ** 2) / expected_f

            dof = len(states) - 1
            chi2_crit = chi2.ppf(0.95, dof)
            passed = chi2_val <= chi2_crit
            return emp_probs, emp_mean, emp_var, err_mean, err_var, chi2_val, chi2_crit, passed

        log_data = "=== Проверка критерия Хи-квадрат при разных N ===\n"
        for test_n in [10, 100, 1000, 10000]:
            _, _, _, err_m, err_v, chi_v, chi_c, psd = simulate_logic(test_n)
            log_data += f"N={test_n:<6} | Ошибки (ср: {err_m:>4.1f}%, дисп: {err_v:>4.1f}%) | Chi2: {chi_v:>6.2f} <= {chi_c:.2f} ({psd})\n"
        self.txt_log.setText(log_data)

        emp_probs, emp_mean, emp_var, err_mean, err_var, chi2_val, chi2_crit, passed = simulate_logic(user_n)
        self.canvas.plot_bars(states, p_list, emp_probs)

        color = "#a6e3a1" if passed else "#f38ba8"
        txt_bool = "TRUE" if passed else "FALSE"
        sign = "&lt;=" if passed else "&gt;"

        self.lbl_result.setText(
            f"Average: <b>{emp_mean:.3f}</b> (отн. погр. = {err_mean:.1f}%)<br>"
            f"Variance: <b>{emp_var:.3f}</b> (отн. погр. = {err_var:.1f}%)<br>"
            f"Chi-squared: {chi2_val:.2f} {sign} {chi2_crit:.3f} is <span style='color:{color}; font-weight:bold;'>{txt_bool}</span>"
        )

class ContinuousTab(QWidget):
    def __init__(self):
        super().__init__()
        self.build_ui()

    def build_ui(self):
        layout = QHBoxLayout(self)

        panel_left = QGroupBox("Параметры нормальной СВ")
        form = QFormLayout(panel_left)
        panel_left.setFixedWidth(300)

        self.in_mean = QLineEdit("0")
        self.in_var = QLineEdit("1")
        self.in_size = QLineEdit("1000")

        form.addRow("Mean (Мат. ож.):", self.in_mean)
        form.addRow("Variance (Дисп.):", self.in_var)
        form.addRow("Sample size (N):", self.in_size)

        form.addRow(QLabel(""), QLabel(""))

        self.btn_run = QPushButton("Запустить моделирование")
        self.btn_run.setMinimumHeight(45)
        self.btn_run.clicked.connect(self.run_simulation)
        form.addRow(self.btn_run)

        layout.addWidget(panel_left)

        panel_right = QVBoxLayout()
        self.canvas = ChartCanvas()
        panel_right.addWidget(self.canvas, stretch=2)

        self.lbl_result = QLabel("Запустите симуляцию для получения результатов...")
        self.lbl_result.setStyleSheet("font-size: 15px; line-height: 1.5;")
        self.lbl_result.setTextFormat(Qt.TextFormat.RichText)
        panel_right.addWidget(self.lbl_result)

        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        panel_right.addWidget(self.txt_log, stretch=1)

        layout.addLayout(panel_right)

    def run_simulation(self):
        try:
            th_mean = float(self.in_mean.text())
            th_var = float(self.in_var.text())
            user_n = int(self.in_size.text())
        except ValueError:
            QMessageBox.warning(self, "Ошибка", "Введите корректные числа.")
            return

        if th_var <= 0:
            QMessageBox.warning(self, "Ошибка", "Дисперсия должна быть строго больше нуля.")
            return

        def simulate_logic(n_size):
            std_dev = math.sqrt(th_var)
            sim_data = []
            for _ in range((n_size + 1) // 2):
                u1, u2 = random.random(), random.random()
                if u1 == 0: u1 = 1e-9
                z0 = math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)
                z1 = math.sqrt(-2.0 * math.log(u1)) * math.sin(2.0 * math.pi * u2)
                sim_data.append(th_mean + z0 * std_dev)
                sim_data.append(th_mean + z1 * std_dev)

            sim_data = sim_data[:n_size]

            emp_mean = np.mean(sim_data)
            emp_var = np.var(sim_data, ddof=1) if n_size > 1 else 0

            if th_mean == 0:
                err_m_str = f"абс. погр. = {abs(emp_mean - th_mean):.3f}"
            else:
                err_m = abs(emp_mean - th_mean) / abs(th_mean) * 100
                err_m_str = f"отн. погр. = {err_m:.1f}%"

            err_v = abs(emp_var - th_var) / th_var * 100

            k_bins = max(3, int(math.log2(n_size) + 1))
            counts, edges = np.histogram(sim_data, bins=k_bins)

            chi2_val = 0
            for i in range(len(counts)):
                p_th = norm.cdf(edges[i + 1], th_mean, std_dev) - norm.cdf(edges[i], th_mean, std_dev)
                expected_f = n_size * p_th
                if expected_f > 0:
                    chi2_val += ((counts[i] - expected_f) ** 2) / expected_f

            dof = max(1, k_bins - 3)
            chi2_crit = chi2.ppf(0.95, dof)
            passed = chi2_val <= chi2_crit

            return sim_data, emp_mean, emp_var, err_m_str, err_v, chi2_val, chi2_crit, passed

        log_data = "=== Проверка критерия Хи-квадрат при разных N ===\n"
        for test_n in [10, 100, 1000, 10000]:
            _, _, _, e_m_str, e_v, c_v, c_c, p_s = simulate_logic(test_n)
            log_data += f"N={test_n:<6} | Ошибки (ср: {e_m_str.split('=')[1].strip():>5}, дисп: {e_v:>4.1f}%) | Chi2: {c_v:>6.2f} <= {c_c:.2f} ({p_s})\n"

        self.txt_log.setText(log_data)

        sim_data, emp_mean, emp_var, err_m_str, err_v, chi2_val, chi2_crit, passed = simulate_logic(user_n)

        self.canvas.plot_histogram(sim_data, th_mean, th_var)

        color = "#a6e3a1" if passed else "#f38ba8"
        txt_bool = "TRUE" if passed else "FALSE"
        sign = "&lt;=" if passed else "&gt;"

        self.lbl_result.setText(
            f"Average: <b>{emp_mean:.3f}</b> ({err_m_str})<br>"
            f"Variance: <b>{emp_var:.3f}</b> (отн. погр. = {err_v:.1f}%)<br>"
            f"Chi-squared: {chi2_val:.2f} {sign} {chi2_crit:.3f} is <span style='color:{color}; font-weight:bold;'>{txt_bool}</span>"
        )


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Лабораторная работа 6: Моделирование случайных величин")
        self.resize(900, 650)

        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)
        layout.setContentsMargins(10, 10, 10, 10)

        tabs = QTabWidget()
        tabs.addTab(DiscreteTab(), "Лаб 06-1 (Дискретная СВ)")
        tabs.addTab(ContinuousTab(), "Лаб 06-2 (Нормальная СВ)")

        layout.addWidget(tabs)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyleSheet(MODERN_STYLE)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())