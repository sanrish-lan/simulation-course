import sys
import csv
import math
import random
import numpy as np
from datetime import datetime, timedelta
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QSlider, QGroupBox, QGridLayout,
    QTableWidget, QTableWidgetItem, QDoubleSpinBox, QMessageBox, QHeaderView
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal, QPointF, QRect
from PyQt6.QtGui import (
    QPainter, QPen, QBrush, QColor, QFont, QLinearGradient,
    QRadialGradient, QPainterPath, QPalette
)

COLORS = {
    "bg": "#0d1117",
    "surface": "#161b22",
    "surface2": "#21262d",
    "border": "#30363d",
    "accent": "#58a6ff",
    "accent2": "#f78166",
    "accent3": "#3fb950",
    "text": "#e6edf3",
    "text2": "#8b949e",
    "sunny": "#ffa726",
    "cloudy": "#78909c",
    "overcast": "#546e7a"
}

STATE_NAMES = ["Ясно", "Облачно", "Пасмурно"]
STATE_EMOJI = ["☀️", "⛅", "🌧"]
STATE_COLORS = [COLORS["sunny"], COLORS["cloudy"], COLORS["accent"]]
STATE_IDS = [1, 2, 3]


def stationary_distribution_ctmc(q_matrix):
    n = len(q_matrix)
    A = np.array(q_matrix).T.tolist()
    A[-1] = [1.0] * n
    b = [0.0] * (n - 1) + [1.0]
    try:
        pi, _, _, _ = np.linalg.lstsq(np.array(A), np.array(b), rcond=None)
        pi = pi.clip(0)
        return pi / np.sum(pi)
    except:
        return np.array([1 / 3, 1 / 3, 1 / 3])


def get_next_state(current, q_matrix):
    row = q_matrix[current]
    lambda_out = -row[current]
    if lambda_out <= 0.00001:
        return current
    r = random.random()
    cumulative = 0.0
    for j, rate in enumerate(row):
        if j == current: continue
        cumulative += (rate / lambda_out)
        if r <= cumulative:
            return j
    return current


class MarkovCanvas(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(320)
        self.q_matrix = [[0] * 3] * 3
        self.current_state = 0
        self.node_positions = []
        self._anim_progress = 0.0
        self._anim_from = 0
        self._anim_to = 0
        self._anim_timer = QTimer(self)
        self._anim_timer.timeout.connect(self._tick_anim)

    def set_matrix(self, m):
        self.q_matrix = m
        self.update()

    def set_state(self, s, prev=None):
        if prev is not None and prev != s:
            self._anim_from = prev
            self._anim_to = s
            self._anim_progress = 0.0
            self._anim_timer.start(16)
        self.current_state = s
        self.update()

    def _tick_anim(self):
        self._anim_progress = min(1.0, self._anim_progress + 0.06)
        self.update()
        if self._anim_progress >= 1.0:
            self._anim_timer.stop()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor(COLORS["surface"]))

        p.setPen(QColor(COLORS["text2"]))
        p.setFont(QFont("Consolas", 9))
        p.drawText(10, 18, "ГРАФ (ИНТЕНСИВНОСТИ ПЕРЕХОДОВ)")

        cx, cy, R = w // 2, h // 2 + 10, min(w, h) // 2 - 60
        n = 3
        angles = [math.pi / 2 + 2 * math.pi * i / n for i in range(n)]
        self.node_positions = [(cx + int(R * math.cos(a)), cy - int(R * math.sin(a))) for a in angles]

        self._draw_edges(p)
        if self._anim_timer.isActive() or self._anim_progress < 1.0:
            self._draw_anim_arrow(p)
        self._draw_nodes(p)

    def _draw_edges(self, painter):
        max_rate = max(0.001, max(max(r) for r in self.q_matrix))
        for i in range(3):
            for j in range(3):
                if i == j: continue
                rate = self.q_matrix[i][j]
                if rate < 0.001: continue

                x1, y1 = self.node_positions[i]
                x2, y2 = self.node_positions[j]

                ratio = rate / max_rate
                alpha = int(80 + 175 * ratio)
                col = QColor(STATE_COLORS[i])
                col.setAlpha(alpha)
                pen = QPen(col, 1 + ratio * 4)
                painter.setPen(pen)
                painter.setBrush(Qt.BrushStyle.NoBrush)

                dx, dy = x2 - x1, y2 - y1
                length = math.sqrt(dx * dx + dy * dy)
                mx = (x1 + x2) / 2 + (-dy / length) * 30
                my = (y1 + y2) / 2 + (dx / length) * 30

                path = QPainterPath()
                path.moveTo(x1, y1)
                path.quadTo(mx, my, x2, y2)
                painter.drawPath(path)

                t = 0.85
                ax = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * mx + t ** 2 * x2
                ay = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * my + t ** 2 * y2
                bx = (1 - (t + 0.05)) ** 2 * x1 + 2 * (1 - (t + 0.05)) * (t + 0.05) * mx + (t + 0.05) ** 2 * x2
                by = (1 - (t + 0.05)) ** 2 * y1 + 2 * (1 - (t + 0.05)) * (t + 0.05) * my + (t + 0.05) ** 2 * y2

                ang = math.atan2(by - ay, bx - ax)
                size = 8
                p1 = QPointF(bx - size * math.cos(ang - 0.4), by - size * math.sin(ang - 0.4))
                p2 = QPointF(bx - size * math.cos(ang + 0.4), by - size * math.sin(ang + 0.4))
                arr_path = QPainterPath()
                arr_path.moveTo(bx, by)
                arr_path.lineTo(p1)
                arr_path.lineTo(p2)
                arr_path.closeSubpath()
                painter.setBrush(QBrush(col))
                painter.setPen(Qt.PenStyle.NoPen)
                painter.drawPath(arr_path)

                mid_t = 0.5
                lx = (1 - mid_t) ** 2 * x1 + 2 * (1 - mid_t) * mid_t * mx + mid_t ** 2 * x2
                ly = (1 - mid_t) ** 2 * y1 + 2 * (1 - mid_t) * mid_t * my + mid_t ** 2 * y2
                painter.setPen(QColor(COLORS["text2"]))
                painter.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
                painter.drawText(int(lx) - 15, int(ly) - 4, f"{rate:.2f}")

    def _draw_anim_arrow(self, painter):
        if self._anim_from == self._anim_to: return
        i, j = self._anim_from, self._anim_to
        x1, y1 = self.node_positions[i]
        x2, y2 = self.node_positions[j]
        t = self._anim_progress
        px = x1 + (x2 - x1) * t
        py = y1 + (y2 - y1) * t
        glow = QRadialGradient(px, py, 14)
        glow.setColorAt(0, QColor(255, 255, 255, 200))
        glow.setColorAt(1, QColor(255, 255, 255, 0))
        painter.setBrush(QBrush(glow))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(int(px) - 14, int(py) - 14, 28, 28)

    def _draw_nodes(self, painter):
        node_r = 34
        for i, (nx, ny) in enumerate(self.node_positions):
            is_active = (i == self.current_state)
            col = QColor(STATE_COLORS[i])
            if is_active:
                glow = QRadialGradient(nx, ny, node_r * 2)
                gc = QColor(col)
                gc.setAlpha(80)
                glow.setColorAt(0, gc)
                glow.setColorAt(1, QColor(0, 0, 0, 0))
                painter.setBrush(QBrush(glow))
                painter.setPen(Qt.PenStyle.NoPen)
                painter.drawEllipse(nx - node_r * 2, ny - node_r * 2, node_r * 4, node_r * 4)

            bg = QColor(COLORS["surface2"])
            grad = QRadialGradient(nx - node_r // 3, ny - node_r // 3, node_r)
            grad.setColorAt(0, col.lighter(130) if is_active else col.darker(120))
            grad.setColorAt(1, bg)
            painter.setBrush(QBrush(grad))
            border_col = col if is_active else col.darker(150)
            border_col.setAlpha(200 if is_active else 120)
            painter.setPen(QPen(border_col, 2 if is_active else 1))
            painter.drawEllipse(nx - node_r, ny - node_r, node_r * 2, node_r * 2)

            painter.setPen(QColor(COLORS["text"]))
            painter.setFont(QFont("Segoe UI Emoji", 16))
            painter.drawText(QRect(nx - node_r, ny - node_r, node_r * 2, node_r - 2), Qt.AlignmentFlag.AlignCenter,
                             STATE_EMOJI[i])
            painter.setFont(QFont("Consolas", 8, QFont.Weight.Bold))
            painter.setPen(QColor(COLORS["text"]) if is_active else QColor(COLORS["text2"]))
            painter.drawText(QRect(nx - node_r, ny + 4, node_r * 2, node_r - 4), Qt.AlignmentFlag.AlignCenter,
                             STATE_NAMES[i])


class WeatherWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(260, 320)
        self.state = 0
        self.day = 0
        self.history = []
        self._pulse = 0.0
        self._pulse_dir = 1
        self._pulse_timer = QTimer(self)
        self._pulse_timer.timeout.connect(self._tick_pulse)
        self._pulse_timer.start(40)

    def set_state(self, state, day, history):
        self.state = state
        self.day = day
        self.history = history[-7:]
        self.update()

    def _tick_pulse(self):
        self._pulse += 0.04 * self._pulse_dir
        if self._pulse >= 1.0:
            self._pulse_dir = -1
        elif self._pulse <= 0.0:
            self._pulse_dir = 1
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()

        bg_colors = [(QColor("#1a2a4a"), QColor("#0d1117")), (QColor("#1a2535"), QColor("#0d1117")),
                     (QColor("#0f1a23"), QColor("#0d1117"))]
        grad = QLinearGradient(0, 0, 0, h)
        grad.setColorAt(0, bg_colors[self.state][0])
        grad.setColorAt(1, bg_colors[self.state][1])
        path = QPainterPath()
        path.addRoundedRect(0, 0, w, h, 20, 20)
        p.fillPath(path, QBrush(grad))

        border_col = QColor(STATE_COLORS[self.state])
        border_col.setAlpha(int(120 + 60 * self._pulse))
        p.setPen(QPen(border_col, 1.5))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawPath(path)

        p.setPen(QColor(COLORS["text2"]))
        p.setFont(QFont("Consolas", 11, QFont.Weight.Bold))
        p.drawText(0, 28, w, 20, Qt.AlignmentFlag.AlignCenter, f"День: {self.day}")

        emoji_size = int(56 + 8 * self._pulse * (1 if self.state == 0 else 0))
        p.setFont(QFont("Segoe UI Emoji", emoji_size))
        p.setPen(QColor(COLORS["text"]))
        p.drawText(QRect(0, 40, w, 100), Qt.AlignmentFlag.AlignCenter, STATE_EMOJI[self.state])

        p.setPen(QColor(STATE_COLORS[self.state]))
        p.setFont(QFont("Consolas", 18, QFont.Weight.Bold))
        p.drawText(QRect(0, 145, w, 30), Qt.AlignmentFlag.AlignCenter, STATE_NAMES[self.state])

        p.setPen(QColor(COLORS["text"]))
        p.setFont(QFont("Consolas", 28, QFont.Weight.Bold))
        p.drawText(QRect(0, 175, w, 44), Qt.AlignmentFlag.AlignCenter, f"{[22, 15, 10][self.state]}°C")

        if len(self.history) > 1:
            self._draw_mini_forecast(p, w, h)

    def _draw_mini_forecast(self, p, w, h):
        hist = self.history[-6:]
        n = len(hist)
        cell_w = (w - 40) // max(n, 1)
        p.setFont(QFont("Segoe UI Emoji", 13))
        for k, s in enumerate(hist):
            x = 20 + k * cell_w + cell_w // 2
            col = QColor(STATE_COLORS[s])
            col.setAlpha(255 if k == n - 1 else int(80 + 140 * (k + 1) / n))
            p.setPen(col)
            p.drawText(QRect(x - cell_w // 2, 228, cell_w, 26), Qt.AlignmentFlag.AlignCenter, STATE_EMOJI[s])


class DistributionChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(180)
        self.empirical = [1 / 3, 1 / 3, 1 / 3]
        self.theoretical = [1 / 3, 1 / 3, 1 / 3]
        self._anim = [0.0, 0.0, 0.0]
        self._target = [1 / 3, 1 / 3, 1 / 3]
        self._anim_timer = QTimer(self)
        self._anim_timer.timeout.connect(self._tick)
        self._anim_timer.start(20)

    def update_data(self, emp, theo):
        self.empirical = emp
        self.theoretical = theo
        self._target = emp[:]

    def _tick(self):
        changed = False
        for i in range(3):
            diff = self._target[i] - self._anim[i]
            if abs(diff) > 0.001:
                self._anim[i] += diff * 0.12
                changed = True
        if changed: self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor(COLORS["surface"]))

        p.setPen(QColor(COLORS["text2"]))
        p.setFont(QFont("Consolas", 9))
        p.drawText(10, 16, "РАСПРЕДЕЛЕНИЕ (доля времени)")

        margin_l, margin_r, margin_t, margin_b = 14, 14, 26, 36
        chart_w, chart_h = w - margin_l - margin_r, h - margin_t - margin_b
        group_w, bar_w = chart_w // 3, (chart_w // 3) // 3
        max_val = max(max(self.empirical), max(self.theoretical), 0.01)

        for i in range(3):
            gx = margin_l + i * group_w
            theo_h = int(self._target[i] / max_val * chart_h * 0.9)
            col = QColor(STATE_COLORS[i])
            col.setAlpha(80)
            p.setBrush(QBrush(col))
            p.setPen(Qt.PenStyle.NoPen)
            p.drawRoundedRect(gx + bar_w // 2, margin_t + chart_h - theo_h, bar_w, theo_h, 3, 3)

            emp_h = int(self._anim[i] / max_val * chart_h * 0.9)
            ex = gx + bar_w + bar_w // 2 + 2
            p.setBrush(QBrush(QColor(STATE_COLORS[i])))
            p.drawRoundedRect(ex, margin_t + chart_h - emp_h, bar_w, emp_h, 3, 3)

            p.setPen(QColor(COLORS["text"]))
            p.setFont(QFont("Consolas", 8))
            p.drawText(gx, h - margin_b + 14, group_w, 16, Qt.AlignmentFlag.AlignCenter,
                       STATE_EMOJI[i] + " " + STATE_NAMES[i])

            p.setPen(QColor(COLORS["text2"]))
            p.setFont(QFont("Consolas", 7))
            if emp_h > 12:
                p.drawText(ex - 4, margin_t + chart_h - emp_h - 2, bar_w + 8, 12, Qt.AlignmentFlag.AlignCenter,
                           f"{self._anim[i]:.2f}")


class MatrixEditor(QGroupBox):
    matrix_changed = pyqtSignal(list)

    def __init__(self, parent=None):
        super().__init__("Интенсивности переходов (Q-матрица)", parent)
        self.setStyleSheet(
            f"QGroupBox {{ color: {COLORS['text2']}; font-family: Consolas; border: 1px solid {COLORS['border']}; border-radius: 8px; margin-top: 10px; padding-top: 10px; }}")
        layout = QGridLayout(self)

        self.spins = []
        default_q = [
            [-0.8, 0.6, 0.2],
            [0.4, -0.7, 0.3],
            [0.1, 0.5, -0.6],
        ]

        for j in range(3):
            lbl = QLabel(f"{STATE_EMOJI[j]}")
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(lbl, 0, j + 1)

        for i in range(3):
            lbl = QLabel(f"{STATE_NAMES[i]}")
            lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:10px;")
            layout.addWidget(lbl, i + 1, 0)
            row_spins = []
            for j in range(3):
                sp = QDoubleSpinBox()
                sp.setRange(-100.0, 100.0)
                sp.setSingleStep(0.1)
                sp.setDecimals(2)
                sp.setValue(default_q[i][j])
                sp.setStyleSheet(
                    f"QDoubleSpinBox {{ background: {COLORS['surface2']}; color: {COLORS['text']}; font-family: Consolas; font-size: 11px; }} QDoubleSpinBox:disabled {{ color: {COLORS['accent2']}; }}")
                if i == j:
                    sp.setDisabled(True)
                else:
                    sp.valueChanged.connect(self._on_changed)
                layout.addWidget(sp, i + 1, j + 1)
                row_spins.append(sp)
            self.spins.append(row_spins)

        self._status = QLabel("q_ii = -sum(q_ij)")
        self._status.setStyleSheet(f"color:{COLORS['accent3']}; font-family:Consolas; font-size:9px;")
        layout.addWidget(self._status, 4, 0, 1, 4)

    def _on_changed(self):
        for i in range(3):
            row_sum = sum(self.spins[i][j].value() for j in range(3) if i != j)
            self.spins[i][i].blockSignals(True)
            self.spins[i][i].setValue(-row_sum)
            self.spins[i][i].blockSignals(False)
        self.matrix_changed.emit(self.get_matrix())

    def get_matrix(self):
        return [[sp.value() for sp in row] for row in self.spins]


class StatsTable(QTableWidget):
    def __init__(self, parent=None):
        super().__init__(4, 4, parent)
        self.setHorizontalHeaderLabels(["", "Ясно ☀️", "Облачно ⛅", "Пасмурно 🌧"])
        self.verticalHeader().setVisible(False)
        for i, text in enumerate(["Всего дн.", "Эмпир.", "Теор.", "Отклон."]):
            item = QTableWidgetItem(text)
            item.setForeground(QColor(COLORS["text2"]))
            self.setItem(i, 0, item)
        self.horizontalHeader().setStretchLastSection(True)
        self.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self.setStyleSheet(
            f"QTableWidget {{ background: {COLORS['surface']}; color: {COLORS['text']}; gridline-color: {COLORS['border']}; font-family: Consolas; font-size: 10px; }} QHeaderView::section {{ background: {COLORS['surface2']}; color: {COLORS['text2']}; }}")

    def update_stats(self, exact_time_in_states, total_time, theoretical):
        for j in range(3):
            emp = exact_time_in_states[j] / total_time if total_time > 0 else 0.0
            theo = theoretical[j]
            diff = emp - theo
            col = QColor(COLORS["accent3"]) if abs(diff) < 0.05 else QColor(COLORS["accent2"])

            items = [f"{exact_time_in_states[j]:.1f}", f"{emp:.3f}", f"{theo:.3f}", f"{diff:+.3f}"]
            for r, txt in enumerate(items):
                item = QTableWidgetItem(txt)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                if r == 3: item.setForeground(col)
                self.setItem(r, j + 1, item)


class HistoryChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(90)
        self.history = []

    def add_state(self, s):
        self.history.append(s)
        if len(self.history) > 100:
            self.history = self.history[-100:]
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor(COLORS["surface"]))
        p.setPen(QColor(COLORS["text2"]))
        p.setFont(QFont("Consolas", 9))
        p.drawText(10, 16, "ИСТОРИЯ ПО ДНЯМ")

        if not self.history: return
        ml, mr, mt, mb = 10, 10, 22, 20
        cw, ch = w - ml - mr, h - mt - mb
        n = len(self.history)
        step = cw / max(n - 1, 1)

        for i in range(3):
            gy = mt + ch - int((i / 2) * ch)
            p.setPen(QPen(QColor(COLORS["border"]), 1, Qt.PenStyle.DotLine))
            p.drawLine(ml, gy, w - mr, gy)

        if n > 1:
            for k in range(n - 1):
                x1, y1 = ml + int(k * step), mt + ch - int(self.history[k] / 2 * ch)
                x2, y2 = ml + int((k + 1) * step), mt + ch - int(self.history[k + 1] / 2 * ch)
                col = QColor(STATE_COLORS[self.history[k]])
                col.setAlpha(200)
                p.setPen(QPen(col, 2))
                p.drawLine(x1, y1, x2, y2)

        lx, ly = ml + int((n - 1) * step), mt + ch - int(self.history[-1] / 2 * ch)
        p.setBrush(QBrush(QColor(STATE_COLORS[self.history[-1]])))
        p.setPen(QPen(QColor(COLORS["surface"]), 2))
        p.drawEllipse(lx - 4, ly - 4, 8, 8)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CTMC Модель погоды (Единица времени - 1 день)")
        self.setMinimumSize(1100, 760)

        self.q_matrix = [
            [-0.8, 0.6, 0.2],
            [0.4, -0.7, 0.3],
            [0.1, 0.5, -0.6],
        ]

        self.current_state = 0
        self.current_day = 0
        self.exact_time = 0.0
        self.running = False

        self.exact_time_in_states = [0.0, 0.0, 0.0]
        self.history = []
        self.sim_data = []

        self._step_timer = QTimer(self)
        self._step_timer.timeout.connect(self._step)
        self.speed_ms = 800

        self._update_theoretical()
        self._build_ui()
        self._apply_styles()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(10, 10, 10, 10)

        left = QVBoxLayout()
        title = QLabel("МОДЕЛИРОВАНИЕ ПОГОДЫ ПО ДНЯМ (CTMC)")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet(
            f"color: {COLORS['accent']}; font-family: Consolas; font-size: 13px; font-weight: bold; border-bottom: 1px solid {COLORS['border']};")
        left.addWidget(title)

        self.markov_canvas = MarkovCanvas()
        left.addWidget(self.markov_canvas)

        self.matrix_editor = MatrixEditor()
        self.matrix_editor.matrix_changed.connect(self._on_matrix_changed)
        left.addWidget(self.matrix_editor)

        root.addLayout(left, 38)

        center = QVBoxLayout()
        self.weather_widget = WeatherWidget()
        self.weather_widget.setFixedWidth(280)
        wwrap = QHBoxLayout()
        wwrap.addStretch()
        wwrap.addWidget(self.weather_widget)
        wwrap.addStretch()
        center.addLayout(wwrap)

        ctrl = QHBoxLayout()
        self.start_btn = QPushButton("▶ СТАРТ")
        self.start_btn.clicked.connect(self._toggle_run)
        self.start_btn.setStyleSheet(self._btn_style(COLORS["accent3"], COLORS["bg"], bold=True))
        self.step_btn = QPushButton("⏭ СЛЕД. ДЕНЬ")
        self.step_btn.clicked.connect(self._manual_step)
        self.step_btn.setStyleSheet(self._btn_style(COLORS["surface2"], COLORS["accent"]))
        self.reset_btn = QPushButton("↺ СБРОС")
        self.reset_btn.clicked.connect(self._reset)
        self.reset_btn.setStyleSheet(self._btn_style(COLORS["surface2"], COLORS["accent2"]))
        self.save_btn = QPushButton("💾 CSV")
        self.save_btn.clicked.connect(self._save_csv)
        self.save_btn.setStyleSheet(self._btn_style(COLORS["surface2"], COLORS["accent3"]))
        for b in [self.start_btn, self.step_btn, self.reset_btn, self.save_btn]: ctrl.addWidget(b)
        center.addLayout(ctrl)

        spd_row = QHBoxLayout()
        spd_lbl = QLabel("Скорость:")
        spd_lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:10px;")
        self.speed_slider = QSlider(Qt.Orientation.Horizontal)
        self.speed_slider.setRange(50, 2000)
        self.speed_slider.setValue(800)
        self.speed_slider.setInvertedAppearance(True)
        self.speed_slider.valueChanged.connect(self._on_speed)
        self.speed_slider.setStyleSheet(
            f"QSlider::groove:horizontal {{ background:{COLORS['border']}; height:4px; border-radius:2px; }} QSlider::handle:horizontal {{ background:{COLORS['accent']}; width:14px; height:14px; margin:-5px 0; border-radius:7px; }}")
        self.speed_val_lbl = QLabel("0.8 с / день")
        self.speed_val_lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:10px; width:75px;")
        spd_row.addWidget(spd_lbl)
        spd_row.addWidget(self.speed_slider)
        spd_row.addWidget(self.speed_val_lbl)
        center.addLayout(spd_row)

        root.addLayout(center, 30)

        right = QVBoxLayout()
        self.dist_chart = DistributionChart()
        right.addWidget(self.dist_chart)
        self.hist_chart = HistoryChart()
        right.addWidget(self.hist_chart)

        stats_lbl = QLabel("СТАТИСТИЧЕСКАЯ ОБРАБОТКА (ИДЕАЛЬНАЯ ДОЛЯ ВРЕМЕНИ)")
        stats_lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:9px;")
        right.addWidget(stats_lbl)

        self.stats_table = StatsTable()
        self.stats_table.verticalHeader().setDefaultSectionSize(18)
        self.stats_table.horizontalHeader().setFixedHeight(20)
        self.stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.stats_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.stats_table.setFixedHeight(96)
        right.addWidget(self.stats_table)

        root.addLayout(right, 32)

        self.markov_canvas.set_matrix(self.q_matrix)
        self.markov_canvas.set_state(self.current_state)
        self.weather_widget.set_state(self.current_state, self.current_day, [self.current_state])
        self._update_stats()

    def _btn_style(self, bg, border, bold=False):
        w = "bold" if bold else "normal"
        return f"QPushButton {{ background: {bg}; color: {border}; border: 1px solid {border}; border-radius: 6px; padding: 7px 12px; font-family: Consolas; font-size: 10px; font-weight: {w}; }} QPushButton:hover {{ background: {border}22; }}"

    def _apply_styles(self):
        self.setStyleSheet(f"QMainWindow, QWidget {{ background: {COLORS['bg']}; color: {COLORS['text']}; }}")

    def _update_theoretical(self):
        self.theoretical = stationary_distribution_ctmc(self.q_matrix).tolist()

    def _on_matrix_changed(self, m):
        self.q_matrix = m
        self.markov_canvas.set_matrix(m)
        self._update_theoretical()
        self._update_stats()

    def _toggle_run(self):
        self.running = not self.running
        if self.running:
            self.start_btn.setText("⏸ ПАУЗА")
            self.start_btn.setStyleSheet(self._btn_style(COLORS["accent2"], COLORS["bg"], bold=True))
            self._step_timer.start(self.speed_ms)
        else:
            self.start_btn.setText("▶ СТАРТ")
            self.start_btn.setStyleSheet(self._btn_style(COLORS["accent3"], COLORS["bg"], bold=True))
            self._step_timer.stop()

    def _manual_step(self):
        if not self.running: self._step()

    def _step(self):
        target_time = self.current_day + 1.0
        prev_state = self.current_state

        while self.exact_time < target_time:
            rate_out = -self.q_matrix[self.current_state][self.current_state]

            if rate_out <= 0.00001:
                dt = target_time - self.exact_time
                self.exact_time_in_states[self.current_state] += dt
                self.exact_time = target_time
                break

            dt = random.expovariate(rate_out)

            if self.exact_time + dt > target_time:
                self.exact_time_in_states[self.current_state] += (target_time - self.exact_time)
                self.exact_time = target_time
            else:
                self.exact_time_in_states[self.current_state] += dt
                self.exact_time += dt
                self.current_state = get_next_state(self.current_state, self.q_matrix)

        self.current_day += 1
        self.history.append(self.current_state)

        self.sim_data.append({
            "День": self.current_day,
            "Состояние_ID": STATE_IDS[self.current_state],
            "Состояние": STATE_NAMES[self.current_state]
        })

        self.markov_canvas.set_state(self.current_state, prev_state)
        self.weather_widget.set_state(self.current_state, self.current_day, self.history)
        self.hist_chart.add_state(self.current_state)
        self._update_stats()

    def _update_stats(self):
        total_time = sum(self.exact_time_in_states)
        emp = [t / total_time if total_time > 0 else 1 / 3 for t in self.exact_time_in_states]
        self.dist_chart.update_data(emp, self.theoretical)
        self.stats_table.update_stats(self.exact_time_in_states, total_time, self.theoretical)

    def _reset(self):
        if self.running: self._toggle_run()
        self.current_state = 0
        self.current_day = 0
        self.exact_time = 0.0
        self.exact_time_in_states = [0.0, 0.0, 0.0]
        self.history = []
        self.sim_data = []
        self.markov_canvas.set_state(0)
        self.weather_widget.set_state(0, 0, [0])
        self.hist_chart.history.clear()
        self.hist_chart.update()
        self._update_stats()

    def _on_speed(self, v):
        self.speed_ms = v
        self.speed_val_lbl.setText(f"{v / 1000:.1f} с / день")
        if self.running: self._step_timer.setInterval(v)

    def _save_csv(self):
        if not self.sim_data:
            QMessageBox.warning(self, "Нет данных", "Сначала запустите симуляцию.")
            return
        fname = f"weather_days_ctmc_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        try:
            with open(fname, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["День", "Состояние_ID", "Состояние"])
                writer.writeheader()
                writer.writerows(self.sim_data)
            QMessageBox.information(self, "Сохранено", f"Экспорт завершен:\n{fname}")
        except Exception as e:
            QMessageBox.critical(self, "Ошибка", str(e))


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor(COLORS["bg"]))
    palette.setColor(QPalette.ColorRole.WindowText, QColor(COLORS["text"]))
    palette.setColor(QPalette.ColorRole.Base, QColor(COLORS["surface"]))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(COLORS["surface2"]))
    palette.setColor(QPalette.ColorRole.Button, QColor(COLORS["surface2"]))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor(COLORS["text"]))
    palette.setColor(QPalette.ColorRole.Highlight, QColor(COLORS["accent"]))
    app.setPalette(palette)

    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()