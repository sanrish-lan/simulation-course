import sys
import math
from collections import Counter

import numpy as np

import matplotlib

matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
    QSpinBox,
    QDoubleSpinBox
)


# =====================================================================
# ОСТАВЛЕНО БЕЗ ИЗМЕНЕНИЙ (Вычислительная логика)
# =====================================================================

def simulate_poisson_stream(lmbda, T):
    """Моделирует один прогон. Возвращает список моментов времени поступления заявок."""
    t = 0.0
    arrivals = []
    while True:
        t += -np.log(np.random.rand()) / lmbda
        if t > T:
            break
        arrivals.append(t)
    return arrivals


def get_theoretical_probs(lmbda, T, max_k):
    """Возвращает словарь с теоретическими вероятностями Пуассона до k = max_k"""
    mu = lmbda * T
    probs = {}
    for k in range(max_k + 1):
        try:
            probs[k] = (math.pow(mu, k) * math.exp(-mu)) / math.factorial(k)
        except OverflowError:
            probs[k] = 0.0
    return probs


def calculate_empirical_stats(all_counts):
    """Считает эмпирическое среднее, дисперсию и частоты"""
    emp_mean = np.mean(all_counts)
    emp_var = np.var(all_counts, ddof=1) if len(all_counts) > 1 else 0.0
    freqs = dict(Counter(all_counts))
    return emp_mean, emp_var, freqs


# =====================================================================
# ТАБЛИЦА СТИЛЕЙ (QSS)
# =====================================================================

APP_STYLE = """
QMainWindow {
    background-color: #090c15;
}
QWidget {
    color: #f3f4f6;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    font-size: 13px;
}
/* Прозрачный фон для базовых элементов решает проблему черных прямоугольников */
QLabel {
    background-color: transparent;
}
QFrame#Card, QFrame#ConclusionCard {
    background-color: #101622;
    border: 1px solid #1d283d;
    border-radius: 12px;
}
QFrame#StatCard {
    background-color: #172030;
    border: 1px solid #26354e;
    border-radius: 8px;
}
QLabel#StatTitle {
    color: #94a3b8;
    font-size: 10px;
    font-weight: 700;
}
QLabel#StatValue {
    color: #3b82f6;
    font-size: 18px;
    font-weight: 800;
}
QPushButton#RunBtn {
    background-color: #3b82f6;
    color: #ffffff;
    border: none;
    border-radius: 8px;
    font-size: 13px;
    font-weight: 700;
    padding: 10px 16px;
}
QPushButton#RunBtn:hover {
    background-color: #2563eb;
}
QPushButton#RunBtn:pressed {
    background-color: #1d4ed8;
}
QPushButton#RunBtn:disabled {
    background-color: #1e293b;
    color: #64748b;
}

/* Стилизация встроенной панели навигации Matplotlib */
QToolBar#ChartToolbar {
    background-color: #101622;
    border: none;
    border-bottom: 1px solid #1d283d;
    spacing: 4px;
    padding: 2px;
}
QToolBar#ChartToolbar QToolButton {
    background-color: transparent;
    border: none;
    border-radius: 4px;
    padding: 4px;
}
QToolBar#ChartToolbar QToolButton:hover {
    background-color: #172030;
}
QToolBar#ChartToolbar QToolButton:pressed {
    background-color: #26354e;
}
QToolBar#ChartToolbar QLabel {
    color: #94a3b8;
    font-size: 11px;
    padding-left: 8px;
}
"""


# =====================================================================
# КОМПОНЕНТЫ ИНТЕРФЕЙСА
# =====================================================================

class StatCard(QFrame):
    def __init__(self, label):
        super().__init__()
        self.setObjectName("StatCard")
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(4)

        self.title_lbl = QLabel(label.upper())
        self.title_lbl.setObjectName("StatTitle")
        self.title_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.title_lbl)

        self.value_lbl = QLabel("—")
        self.value_lbl.setObjectName("StatValue")
        self.value_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.value_lbl)

    def set_value(self, val):
        self.value_lbl.setText(val)


class ParamField(QWidget):
    def __init__(self, label, min_val, max_val, default, is_float=False):
        super().__init__()
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 4, 0, 4)
        layout.setSpacing(10)

        lbl = QLabel(label)
        lbl.setStyleSheet("color: #94a3b8; font-weight: 500;")
        layout.addWidget(lbl)

        if is_float:
            self.spin = QDoubleSpinBox()
            self.spin.setDecimals(2)
            self.spin.setSingleStep(0.5)
        else:
            self.spin = QSpinBox()
            self.spin.setSingleStep(100)

        self.spin.setRange(min_val, max_val)
        self.spin.setValue(default)
        self.spin.setFixedWidth(90)
        self.spin.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        self.spin.setStyleSheet("""
            QDoubleSpinBox, QSpinBox {
                background-color: #172030;
                border: 1px solid #26354e;
                border-radius: 6px;
                color: #f3f4f6;
                padding: 4px 6px;
                font-weight: 600;
            }
            QDoubleSpinBox:focus, QSpinBox:focus {
                border: 1px solid #3b82f6;
            }
        """)
        layout.addWidget(self.spin)

    def get(self):
        return self.spin.value()


class PoissonChart(FigureCanvas):
    def __init__(self):
        self.fig = Figure(facecolor="#101622", edgecolor='none')
        self.ax = self.fig.add_subplot(111, facecolor="#101622")
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        # Предохранитель рекурсии и инициализация эталонных границ
        self._lock_active = False
        self._prev_xlim = (0.0, 1.0)
        self._prev_ylim = (0.0, 1.0)

        self._style_ax_empty()
        self.draw()

        # Связываем событие прокрутки мыши для зума
        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)

    def _style_base(self):
        self.ax.tick_params(colors="#94a3b8", labelsize=9)
        for spine in self.ax.spines.values():
            spine.set_color("#26354e")
            spine.set_linewidth(1.0)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.grid(axis="y", color="#172030", linewidth=1, linestyle=":")
        self.ax.xaxis.label.set_color("#94a3b8")
        self.ax.yaxis.label.set_color("#94a3b8")

    def _style_ax_empty(self):
        self.ax.clear()
        self.ax.set_title("Запустите симуляцию для построения распределения", color="#94a3b8", fontsize=11, pad=14)
        self.ax.set_xlabel("Число событий k", color="#64748b", fontsize=9)
        self.ax.set_ylabel("Вероятность P(X = k)", color="#64748b", fontsize=9)
        self._style_base()

        # Фиксируем границы по умолчанию
        self._prev_xlim = self.ax.get_xlim()
        self._prev_ylim = self.ax.get_ylim()

        # Подключаем слежение на пустом графике
        self.ax.callbacks.connect('xlim_changed', self.on_xlim_changed)
        self.ax.callbacks.connect('ylim_changed', self.on_ylim_changed)

    def draw_distribution(self, freqs, N, mu, lmbda, T):
        # 1. Очищаем холст (это также сбросит старые callbacks, предотвращая конфликты автоскейла)
        self.ax.clear()

        # 2. Строим распределение
        max_k = max(freqs.keys()) if freqs else 0
        ks = list(range(max_k + 2))
        emp = [freqs.get(k, 0) / N for k in ks]
        th_dict = get_theoretical_probs(mu, 1, max_k + 1)
        th = [th_dict.get(k, 0.0) for k in ks]

        self.ax.bar(
            ks,
            emp,
            width=0.6,
            color="#3b82f6",
            edgecolor="#60a5fa",
            linewidth=1.0,
            alpha=0.45,
            label="Эмпирическое",
            zorder=2,
        )
        self.ax.plot(
            ks,
            th,
            color="#f43f5e",
            linewidth=2.2,
            marker="o",
            markersize=5,
            markerfacecolor="#ffffff",
            markeredgecolor="#f43f5e",
            markeredgewidth=1.5,
            label="Теоретическое (Пуассон)",
            zorder=3,
        )

        # 3. Текстовое оформление и легенда
        self.ax.set_title(
            f"Распределение вероятностей событий  (λ={lmbda}, T={T}с, N={N})",
            color="#f3f4f6",
            fontsize=12,
            fontweight="bold",
            pad=16,
        )
        self.ax.set_xlabel("Число событий k", color="#94a3b8", fontsize=10, labelpad=8)
        self.ax.set_ylabel("Вероятность P(X = k)", color="#94a3b8", fontsize=10, labelpad=8)
        self._style_base()
        self.ax.legend(
            fontsize=9,
            facecolor="#172030",
            edgecolor="#26354e",
            labelcolor="#f3f4f6",
            framealpha=0.9,
            loc="upper right"
        )

        # 4. Корректируем автоматические отступы Matplotlib (срезаем уход в минус)
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        # Задаем жесткий нижний порог (-0.5 для X, чтобы нулевой столбец не обрезался, и 0.0 для Y)
        new_xmin = max(-0.5, xlim[0])
        new_ymin = max(0.0, ylim[0])

        self.ax.set_xlim(new_xmin, xlim[1])
        self.ax.set_ylim(new_ymin, ylim[1])

        # 5. Запоминаем эти идеальные стартовые границы как эталон
        self._prev_xlim = self.ax.get_xlim()
        self._prev_ylim = self.ax.get_ylim()

        # 6. Только теперь заново подключаем слежение за границами для интерактивного панорирования/зума
        self.ax.callbacks.connect('xlim_changed', self.on_xlim_changed)
        self.ax.callbacks.connect('ylim_changed', self.on_ylim_changed)

        # 7. Финальное обновление холста
        self.draw()

    def on_xlim_changed(self, event):
        """Блокирует смещение по оси X влево дальше нуля (с учетом красивого отступа -0.5)"""
        if self._lock_active:
            return
        xlim = self.ax.get_xlim()
        if xlim[0] < -0.5:
            self._lock_active = True
            self.ax.set_xlim(self._prev_xlim)
            self._lock_active = False
        else:
            self._prev_xlim = xlim

    def on_ylim_changed(self, event):
        """Блокирует смещение по оси Y вниз дальше абсолютного нуля"""
        if self._lock_active:
            return
        ylim = self.ax.get_ylim()
        if ylim[0] < -0.01:  # -0.01 — технический допуск для корректной прорисовки оси
            self._lock_active = True
            self.ax.set_ylim(self._prev_ylim)
            self._lock_active = False
        else:
            self._prev_ylim = ylim

    def on_scroll(self, event):
        """Плавное интерактивное масштабирование колесиком мыши"""
        if event.inaxes != self.ax:
            return

        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        base_scale = 1.25
        if event.button == 'up':
            scale_factor = 1 / base_scale
        elif event.button == 'down':
            scale_factor = base_scale
        else:
            scale_factor = 1.0

        new_width = (xlim[1] - xlim[0]) * scale_factor
        new_height = (ylim[1] - ylim[0]) * scale_factor

        xdata = event.xdata
        ydata = event.ydata

        rel_x = (xlim[1] - xdata) / (xlim[1] - xlim[0])
        rel_y = (ylim[1] - ydata) / (ylim[1] - ylim[0])

        self.ax.set_xlim([xdata - new_width * (1 - rel_x), xdata + new_width * rel_x])
        self.ax.set_ylim([ydata - new_height * (1 - rel_y), ydata + new_height * rel_y])
        self.draw()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Симулятор пуассоновского потока")
        self.resize(1150, 750)
        self.setMinimumSize(960, 660)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(16)

        # -----------------------------------------------------------------
        # ЛЕВАЯ ПАНЕЛЬ: УПРАВЛЕНИЕ И СТАТИСТИКА (Без лишних верхних заголовков)
        # -----------------------------------------------------------------
        sidebar = QVBoxLayout()
        sidebar.setContentsMargins(0, 0, 0, 0)
        sidebar.setSpacing(16)

        # Карточка параметров
        ctrl_card = QFrame()
        ctrl_card.setObjectName("Card")
        ctrl_layout = QVBoxLayout(ctrl_card)
        ctrl_layout.setContentsMargins(16, 16, 16, 16)
        ctrl_layout.setSpacing(12)

        ctrl_title = QLabel("ПАРАМЕТРЫ СИМУЛЯЦИИ")
        ctrl_title.setStyleSheet("font-size: 11px; font-weight: bold; color: #818cf8; letter-spacing: 0.5px;")
        ctrl_layout.addWidget(ctrl_title)

        self.p_lambda = ParamField("Интенсивность λ (зап/с)", 0.1, 500.0, 5.0, is_float=True)
        self.p_T = ParamField("Интервал T (секунд)", 0.1, 500.0, 3.0, is_float=True)
        self.p_N = ParamField("Число прогонов N", 10, 100000, 1000, is_float=False)

        ctrl_layout.addWidget(self.p_lambda)
        ctrl_layout.addWidget(self.p_T)
        ctrl_layout.addWidget(self.p_N)

        self.btn = QPushButton("Запустить симуляцию")
        self.btn.setObjectName("RunBtn")
        self.btn.clicked.connect(self.run)
        ctrl_layout.addWidget(self.btn)

        sidebar.addWidget(ctrl_card)

        # Карточка статистических результатов
        stats_card = QFrame()
        stats_card.setObjectName("Card")
        stats_layout = QVBoxLayout(stats_card)
        stats_layout.setContentsMargins(16, 16, 16, 16)
        stats_layout.setSpacing(12)

        stats_title = QLabel("ПОКАЗАТЕЛИ СХОДИМОСТИ")
        stats_title.setStyleSheet("font-size: 11px; font-weight: bold; color: #818cf8; letter-spacing: 0.5px;")
        stats_layout.addWidget(stats_title)

        stats_grid = QGridLayout()
        stats_grid.setSpacing(8)

        self.card_mu = StatCard("μ = λT (теор.)")
        self.card_mean = StatCard("Среднее (эмп.)")
        self.card_varth = StatCard("Дисперсия (теор.)")
        self.card_var = StatCard("Дисперсия (эмп.)")

        stats_grid.addWidget(self.card_mu, 0, 0)
        stats_grid.addWidget(self.card_mean, 0, 1)
        stats_grid.addWidget(self.card_varth, 1, 0)
        stats_grid.addWidget(self.card_var, 1, 1)

        stats_layout.addLayout(stats_grid)
        sidebar.addWidget(stats_card)

        sidebar.addStretch()
        main_layout.addLayout(sidebar, 0)

        # -----------------------------------------------------------------
        # ПРАВАЯ ПАНЕЛЬ: ГРАФИК И ТЕКСТОВЫЙ ВЫВОД
        # -----------------------------------------------------------------
        right_panel = QVBoxLayout()
        right_panel.setSpacing(16)

        # Контейнер для графика
        chart_frame = QFrame()
        chart_frame.setObjectName("Card")
        chart_layout = QVBoxLayout(chart_frame)
        chart_layout.setContentsMargins(12, 12, 12, 12)
        chart_layout.setSpacing(6)

        self.chart = PoissonChart()

        # Интеграция панели инструментов навигации Matplotlib (Zoom / Pan)
        self.toolbar = NavigationToolbar(self.chart, self)
        self.toolbar.setObjectName("ChartToolbar")

        chart_layout.addWidget(self.toolbar)
        chart_layout.addWidget(self.chart)
        right_panel.addWidget(chart_frame, 1)

        # Панель вывода аналитического заключения
        conclusion_frame = QFrame()
        conclusion_frame.setObjectName("ConclusionCard")
        conclusion_layout = QHBoxLayout(conclusion_frame)
        conclusion_layout.setContentsMargins(16, 14, 16, 14)
        conclusion_layout.setSpacing(12)

        self.status_dot = QLabel("●")
        self.status_dot.setStyleSheet("color: #4b5563; font-size: 16px;")
        conclusion_layout.addWidget(self.status_dot, 0, Qt.AlignmentFlag.AlignTop)

        self.conclusion = QLabel(
            "После запуска симуляции здесь появится аналитический вывод по результатам эксперимента.")
        self.conclusion.setStyleSheet("color: #94a3b8; font-size: 13px; line-height: 1.4;")
        self.conclusion.setWordWrap(True)
        conclusion_layout.addWidget(self.conclusion, 1)

        right_panel.addWidget(conclusion_frame)
        main_layout.addLayout(right_panel, 1)

    def run(self):
        lmbda = self.p_lambda.get()
        T = self.p_T.get()
        N = self.p_N.get()

        if lmbda <= 0 or T <= 0 or N <= 0:
            QMessageBox.warning(
                self,
                "Некорректные параметры",
                "Интенсивность λ, интервал T и число прогонов N должны быть больше нуля.",
            )
            return

        self.btn.setEnabled(False)
        self.btn.setText("Вычисление результатов...")
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        QApplication.processEvents()

        try:
            mu = lmbda * T
            all_counts = [len(simulate_poisson_stream(lmbda, T)) for _ in range(N)]
            emp_mean, emp_var, freqs = calculate_empirical_stats(all_counts)

            self.card_mu.set_value(f"{mu:.2f}")
            self.card_mean.set_value(f"{emp_mean:.2f}")
            self.card_varth.set_value(f"{mu:.2f}")
            self.card_var.set_value(f"{emp_var:.2f}")

            self.chart.draw_distribution(freqs, N, mu, lmbda, T)
            self._conclude(lmbda, T, N, mu, emp_mean, emp_var)
        finally:
            QApplication.restoreOverrideCursor()
            self.btn.setEnabled(True)
            self.btn.setText("Запустить симуляцию")

    def _conclude(self, lmbda, T, N, mu, emp_mean, emp_var):
        diff_m = abs(emp_mean - mu)
        diff_v = abs(emp_var - mu)
        ok = diff_m < 0.4 and diff_v < 0.8

        verdict = (
            "Результаты хорошо согласуются с теорией ✓"
            if ok
            else "Для лучшего соответствия теории рекомендуется увеличить число прогонов N."
        )

        text = (
            f"Успешно смоделировано {N} интервалов длительностью T = {T} с при интенсивности λ = {lmbda} соб/с.\n"
            f"Теоретическое среднее μ = λT = {mu:.2f}. "
            f"Получено эмпирическое среднее = {emp_mean:.2f} (отклонение: {diff_m:.2f}), "
            f"эмпирическая дисперсия = {emp_var:.2f} (отклонение: {diff_v:.2f}).\n"
            f"Для закона Пуассона математическое ожидание равно дисперсии (M[X] = D[X] = μ). {verdict}"
        )

        dot_color = "#10b981" if ok else "#f59e0b"
        self.status_dot.setStyleSheet(f"color: {dot_color}; font-size: 16px;")
        self.conclusion.setText(text)
        self.conclusion.setStyleSheet("color: #e5e7eb; font-size: 13px; line-height: 1.45;")


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setStyleSheet(APP_STYLE)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()