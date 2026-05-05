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
    "bg":        "#0d1117",
    "surface":   "#161b22",
    "surface2":  "#21262d",
    "border":    "#30363d",
    "accent":    "#58a6ff",
    "accent2":   "#f78166",
    "accent3":   "#3fb950",
    "text":      "#e6edf3",
    "text2":     "#8b949e",
    "sunny":     "#ffa726",
    "cloudy":    "#78909c",
    "overcast":  "#546e7a",
    "sunny_bg":  "#1a1400",
    "cloudy_bg": "#0d1520",
    "ovrcast_bg":"#0a0f14",
}

STATE_NAMES  = ["Ясно", "Облачно", "Пасмурно"]
STATE_EMOJI  = ["☀️",   "⛅",      "🌧"]
STATE_COLORS = [COLORS["sunny"], COLORS["cloudy"], COLORS["accent"]]
STATE_IDS    = [1, 2, 3]

def stationary_distribution(trans_matrix):
    n = len(trans_matrix)
    A = (np.array(trans_matrix).T - np.eye(n)).tolist()
    A.append([1.0] * n)
    b = [0.0] * n + [1.0]
    A = np.array(A)
    b = np.array(b)
    pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return pi.clip(0)


def next_state(current, trans_matrix):
    row = trans_matrix[current]
    r = random.random()
    cumulative = 0.0
    for j, p in enumerate(row):
        cumulative += p
        if r <= cumulative:
            return j
    return len(row) - 1


class MarkovCanvas(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(320)
        self.trans_matrix = [[1/3, 1/3, 1/3]] * 3
        self.current_state = 0
        self.node_positions = []
        self._anim_progress = 0.0
        self._anim_from = 0
        self._anim_to   = 0

        self._anim_timer = QTimer(self)
        self._anim_timer.timeout.connect(self._tick_anim)

    def set_matrix(self, m):
        self.trans_matrix = m
        self.update()

    def set_state(self, s, prev=None):
        if prev is not None and prev != s:
            self._anim_from = prev
            self._anim_to   = s
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
        f = QFont("Consolas", 9)
        p.setFont(f)
        p.drawText(10, 18, "ГРАФ МАРКОВСКОЙ ЦЕПИ")

        cx, cy, R = w // 2, h // 2 + 10, min(w, h) // 2 - 60
        n = 3
        angles = [math.pi / 2 + 2 * math.pi * i / n for i in range(n)]
        self.node_positions = [
            (cx + int(R * math.cos(a)), cy - int(R * math.sin(a)))
            for a in angles
        ]

        self._draw_edges(p)
        if self._anim_timer.isActive() or self._anim_progress < 1.0:
            self._draw_anim_arrow(p)
        self._draw_nodes(p)

    def _draw_edges(self, painter):
        for i in range(3):
            for j in range(3):
                prob = self.trans_matrix[i][j]
                if prob < 0.001:
                    continue
                x1, y1 = self.node_positions[i]
                x2, y2 = self.node_positions[j]
                if i == j:
                    self._draw_self_loop(painter, x1, y1, prob, STATE_COLORS[i])
                    continue
                alpha = int(60 + 180 * prob)
                col = QColor(STATE_COLORS[i])
                col.setAlpha(alpha)
                pen = QPen(col, 1 + prob * 3)
                pen.setStyle(Qt.PenStyle.SolidLine)
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
                ax = (1-t)**2 * x1 + 2*(1-t)*t * mx + t**2 * x2
                ay = (1-t)**2 * y1 + 2*(1-t)*t * my + t**2 * y2
                bx = (1-(t+0.05))**2 * x1 + 2*(1-(t+0.05))*(t+0.05) * mx + (t+0.05)**2 * x2
                by = (1-(t+0.05))**2 * y1 + 2*(1-(t+0.05))*(t+0.05) * my + (t+0.05)**2 * y2
                self._draw_arrowhead(painter, ax, ay, bx, by, col)

                mid_t = 0.5
                lx = (1-mid_t)**2 * x1 + 2*(1-mid_t)*mid_t * mx + mid_t**2 * x2
                ly = (1-mid_t)**2 * y1 + 2*(1-mid_t)*mid_t * my + mid_t**2 * y2
                painter.setPen(QColor(COLORS["text2"]))
                painter.setFont(QFont("Consolas", 8))
                painter.drawText(int(lx) - 15, int(ly) - 4, f"{prob:.2f}")

    def _draw_self_loop(self, painter, x, y, prob, color):
        col = QColor(color)
        col.setAlpha(int(80 + 150 * prob))
        pen = QPen(col, 1 + prob * 2)
        painter.setPen(pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        r = 22
        painter.drawEllipse(x - r, y - 50, r * 2, r * 2)
        painter.setFont(QFont("Consolas", 7))
        painter.setPen(QColor(COLORS["text2"]))
        painter.drawText(x - 12, y - 55, f"{prob:.2f}")

    def _draw_arrowhead(self, painter, x1, y1, x2, y2, color):
        ang = math.atan2(y2 - y1, x2 - x1)
        size = 8
        p1 = QPointF(x2 - size * math.cos(ang - 0.4), y2 - size * math.sin(ang - 0.4))
        p2 = QPointF(x2 - size * math.cos(ang + 0.4), y2 - size * math.sin(ang + 0.4))
        path = QPainterPath()
        path.moveTo(x2, y2)
        path.lineTo(p1)
        path.lineTo(p2)
        path.closeSubpath()
        painter.setBrush(QBrush(color))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawPath(path)

    def _draw_anim_arrow(self, painter):
        if self._anim_from == self._anim_to:
            return
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
                painter.drawEllipse(nx - node_r * 2, ny - node_r * 2,
                                    node_r * 4, node_r * 4)

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
            painter.drawText(QRect(nx - node_r, ny - node_r, node_r * 2, node_r - 2),
                             Qt.AlignmentFlag.AlignCenter, STATE_EMOJI[i])
            painter.setFont(QFont("Consolas", 8, QFont.Weight.Bold))
            painter.setPen(QColor(COLORS["text"]) if is_active else QColor(COLORS["text2"]))
            painter.drawText(QRect(nx - node_r, ny + 4, node_r * 2, node_r - 4),
                             Qt.AlignmentFlag.AlignCenter, STATE_NAMES[i])


class WeatherWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(260, 320)
        self.state = 0
        self.day = 0
        self.date = datetime.now()
        self.history = []
        self._pulse = 0.0
        self._pulse_dir = 1
        self._pulse_timer = QTimer(self)
        self._pulse_timer.timeout.connect(self._tick_pulse)
        self._pulse_timer.start(40)

    def set_state(self, state, day, history):
        self.state = state
        self.day = day
        self.date = datetime.now() + timedelta(days=day)
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

        bg_colors = [
            (QColor("#1a2a4a"), QColor("#0d1117")),
            (QColor("#1a2535"), QColor("#0d1117")),
            (QColor("#0f1a23"), QColor("#0d1117")),
        ]
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
        p.setFont(QFont("Consolas", 9))
        p.drawText(0, 28, w, 20, Qt.AlignmentFlag.AlignCenter,
                   f"День {self.day + 1}  •  {self.date.strftime('%d %b %Y')}")

        emoji_size = int(56 + 8 * self._pulse * (1 if self.state == 0 else 0))
        p.setFont(QFont("Segoe UI Emoji", emoji_size))
        p.setPen(QColor(COLORS["text"]))
        p.drawText(QRect(0, 40, w, 100), Qt.AlignmentFlag.AlignCenter, STATE_EMOJI[self.state])

        col = QColor(STATE_COLORS[self.state])
        p.setPen(col)
        p.setFont(QFont("Consolas", 18, QFont.Weight.Bold))
        p.drawText(QRect(0, 145, w, 30), Qt.AlignmentFlag.AlignCenter, STATE_NAMES[self.state])

        temps = [22, 15, 10]
        base_t = temps[self.state]
        p.setPen(QColor(COLORS["text"]))
        p.setFont(QFont("Consolas", 28, QFont.Weight.Bold))
        p.drawText(QRect(0, 175, w, 44), Qt.AlignmentFlag.AlignCenter, f"{base_t}°C")

        if len(self.history) > 1:
            self._draw_mini_forecast(p, w, h)

        p.setPen(QPen(QColor(COLORS["border"]), 1))
        p.drawLine(20, h - 55, w - 20, h - 55)

        p.setPen(QColor(COLORS["text2"]))
        p.setFont(QFont("Consolas", 8))
        p.drawText(QRect(0, h - 48, w // 2, 20), Qt.AlignmentFlag.AlignCenter, "💧 Влажность")
        p.drawText(QRect(w // 2, h - 48, w // 2, 20), Qt.AlignmentFlag.AlignCenter, "💨 Ветер")
        humid = [40, 65, 85][self.state]
        wind  = [5, 12, 20][self.state]
        p.setPen(QColor(COLORS["text"]))
        p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
        p.drawText(QRect(0, h - 30, w // 2, 20), Qt.AlignmentFlag.AlignCenter, f"{humid}%")
        p.drawText(QRect(w // 2, h - 30, w // 2, 20), Qt.AlignmentFlag.AlignCenter, f"{wind} м/с")

    def _draw_mini_forecast(self, p, w, h):
        hist = self.history[-6:]
        n = len(hist)
        if n == 0:
            return
        cell_w = (w - 40) // max(n, 1)
        y_base = 228
        p.setFont(QFont("Segoe UI Emoji", 13))
        for k, s in enumerate(hist):
            x = 20 + k * cell_w + cell_w // 2
            col = QColor(STATE_COLORS[s])
            if k == n - 1:
                col.setAlpha(255)
            else:
                col.setAlpha(int(80 + 140 * (k + 1) / n))
            p.setPen(col)
            p.drawText(QRect(x - cell_w // 2, y_base, cell_w, 26),
                       Qt.AlignmentFlag.AlignCenter, STATE_EMOJI[s])
            if k == n - 1:
                p.setBrush(QBrush(col))
                p.setPen(Qt.PenStyle.NoPen)
                p.drawEllipse(x - 3, y_base + 28, 6, 6)


class DistributionChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(180)
        self.empirical = [1/3, 1/3, 1/3]
        self.theoretical = [1/3, 1/3, 1/3]
        self._anim = [0.0, 0.0, 0.0]
        self._target = [1/3, 1/3, 1/3]
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
        if changed:
            self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor(COLORS["surface"]))

        p.setPen(QColor(COLORS["text2"]))
        p.setFont(QFont("Consolas", 9))
        p.drawText(10, 16, "РАСПРЕДЕЛЕНИЕ СОСТОЯНИЙ")

        margin_l, margin_r, margin_t, margin_b = 14, 14, 26, 36
        chart_w = w - margin_l - margin_r
        chart_h = h - margin_t - margin_b

        group_w = chart_w // 3
        bar_w   = group_w // 3

        max_val = max(max(self.empirical), max(self.theoretical), 0.01)

        for i in range(3):
            gx = margin_l + i * group_w

            theo_h = int(self._target[i] / max_val * chart_h * 0.9)
            tx = gx + bar_w // 2
            col = QColor(STATE_COLORS[i])
            col.setAlpha(80)
            p.setBrush(QBrush(col))
            p.setPen(Qt.PenStyle.NoPen)
            p.drawRoundedRect(tx, margin_t + chart_h - theo_h, bar_w, theo_h, 3, 3)

            emp_h = int(self._anim[i] / max_val * chart_h * 0.9)
            ex = gx + bar_w + bar_w // 2 + 2
            ecol = QColor(STATE_COLORS[i])
            p.setBrush(QBrush(ecol))
            p.setPen(Qt.PenStyle.NoPen)
            p.drawRoundedRect(ex, margin_t + chart_h - emp_h, bar_w, emp_h, 3, 3)

            p.setPen(QColor(COLORS["text"]))
            p.setFont(QFont("Consolas", 8))
            lbl_y = h - margin_b + 14
            p.drawText(gx, lbl_y, group_w, 16,
                       Qt.AlignmentFlag.AlignCenter, STATE_EMOJI[i] + " " + STATE_NAMES[i])

            p.setPen(QColor(COLORS["text2"]))
            p.setFont(QFont("Consolas", 7))
            if emp_h > 12:
                p.drawText(ex - 4, margin_t + chart_h - emp_h - 2, bar_w + 8, 12,
                           Qt.AlignmentFlag.AlignCenter, f"{self._anim[i]:.2f}")

        p.setFont(QFont("Consolas", 8))
        lx = w - 120
        col2 = QColor(STATE_COLORS[0])
        col2.setAlpha(80)
        p.fillRect(lx, 4, 10, 10, col2)
        p.setPen(QColor(COLORS["text2"]))
        p.drawText(lx + 13, 14, "Теор.")
        p.fillRect(lx + 55, 4, 10, 10, QColor(STATE_COLORS[0]))
        p.drawText(lx + 68, 14, "Эмп.")


class MatrixEditor(QGroupBox):
    matrix_changed = pyqtSignal(list)

    def __init__(self, parent=None):
        super().__init__("Матрица переходов", parent)
        self.setStyleSheet(f"""
            QGroupBox {{
                color: {COLORS['text2']};
                font-family: Consolas;
                font-size: 10px;
                border: 1px solid {COLORS['border']};
                border-radius: 8px;
                margin-top: 10px;
                padding-top: 10px;
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 4px;
            }}
        """)
        layout = QGridLayout(self)
        layout.setSpacing(4)

        self.spins = []
        default = [
            [0.6, 0.3, 0.1],
            [0.2, 0.5, 0.3],
            [0.1, 0.3, 0.6],
        ]

        header_style = f"color:{COLORS['text2']}; font-family:Consolas; font-size:9px;"
        for j in range(3):
            lbl = QLabel(f"{STATE_EMOJI[j]} {STATE_NAMES[j]}")
            lbl.setStyleSheet(header_style)
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(lbl, 0, j + 1)
        for i in range(3):
            lbl = QLabel(f"{STATE_EMOJI[i]} {STATE_NAMES[i]}")
            lbl.setStyleSheet(header_style)
            lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            layout.addWidget(lbl, i + 1, 0)

            row_spins = []
            for j in range(3):
                sp = QDoubleSpinBox()
                sp.setRange(0.0, 1.0)
                sp.setSingleStep(0.05)
                sp.setDecimals(2)
                sp.setValue(default[i][j])
                sp.setStyleSheet(f"""
                    QDoubleSpinBox {{
                        background: {COLORS['surface2']};
                        color: {COLORS['text']};
                        border: 1px solid {COLORS['border']};
                        border-radius: 4px;
                        padding: 2px;
                        font-family: Consolas;
                        font-size: 11px;
                    }}
                    QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {{
                        width: 14px;
                        background: {COLORS['surface']};
                    }}
                """)
                sp.valueChanged.connect(self._on_changed)
                layout.addWidget(sp, i + 1, j + 1)
                row_spins.append(sp)
            self.spins.append(row_spins)

        self._status = QLabel("✓ Сумма строк = 1.0")
        self._status.setStyleSheet(f"color:{COLORS['accent3']}; font-family:Consolas; font-size:9px;")
        layout.addWidget(self._status, 4, 0, 1, 4)

    def _on_changed(self):
        m = self.get_matrix()
        ok = all(abs(sum(row) - 1.0) < 0.01 for row in m)
        if ok:
            self._status.setText("✓ Сумма строк = 1.0")
            self._status.setStyleSheet(f"color:{COLORS['accent3']}; font-family:Consolas; font-size:9px;")
            self.matrix_changed.emit(m)
        else:
            errs = []
            for i, row in enumerate(m):
                s = sum(row)
                if abs(s - 1.0) >= 0.01:
                    errs.append(f"строка {i+1}: {s:.2f}")
            self._status.setText("⚠ " + ", ".join(errs))
            self._status.setStyleSheet(f"color:{COLORS['accent2']}; font-family:Consolas; font-size:9px;")

    def get_matrix(self):
        return [[sp.value() for sp in row] for row in self.spins]

    def normalize(self):
        for row in self.spins:
            s = sum(sp.value() for sp in row)
            if s > 0:
                for sp in row:
                    sp.blockSignals(True)
                    sp.setValue(sp.value() / s)
                    sp.blockSignals(False)
        self._on_changed()


class StatsTable(QTableWidget):
    def __init__(self, parent=None):
        super().__init__(4, 4, parent)
        headers_h = ["", "Ясно ☀️", "Облачно ⛅", "Пасмурно 🌧"]
        self.setHorizontalHeaderLabels(headers_h)

        self.verticalHeader().setVisible(False)

        row_labels = ["Кол-во дней", "Эмпир. доля", "Теор. доля", "Отклонение"]
        for i, text in enumerate(row_labels):
            item = QTableWidgetItem(text)
            item.setForeground(QColor(COLORS["text2"]))
            self.setItem(i, 0, item)

        self.horizontalHeader().setStretchLastSection(True)
        self.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self.setStyleSheet(f"""
            QTableWidget {{
                background: {COLORS['surface']};
                color: {COLORS['text']};
                gridline-color: {COLORS['border']};
                border: 1px solid {COLORS['border']};
                font-family: Consolas;
                font-size: 10px;
            }}
            QHeaderView::section {{
                background: {COLORS['surface2']};
                color: {COLORS['text2']};
                border: 1px solid {COLORS['border']};
                padding: 4px;
                font-family: Consolas;
                font-size: 9px;
            }}
            QTableWidget::item {{ padding: 4px; }}
        """)

    def update_stats(self, counts, total, theoretical):
        for j in range(3):
            emp = counts[j] / total if total > 0 else 0.0
            theo = theoretical[j]
            diff = emp - theo
            col = QColor(COLORS["accent3"]) if abs(diff) < 0.05 else QColor(COLORS["accent2"])

            items = [str(counts[j]), f"{emp:.3f}", f"{theo:.3f}", f"{diff:+.3f}"]
            for r, txt in enumerate(items):
                item = QTableWidgetItem(txt)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                if r == 3:
                    item.setForeground(col)
                self.setItem(r, j + 1, item)

class HistoryChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(90)
        self.history = []

    def add_state(self, s):
        self.history.append(s)
        if len(self.history) > 200:
            self.history = self.history[-200:]
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor(COLORS["surface"]))
        p.setPen(QColor(COLORS["text2"]))
        p.setFont(QFont("Consolas", 9))
        p.drawText(10, 16, "ИСТОРИЯ СОСТОЯНИЙ")

        if not self.history:
            return

        ml, mr, mt, mb = 10, 10, 22, 20
        cw = w - ml - mr
        ch = h - mt - mb
        n = len(self.history)
        step = cw / max(n - 1, 1)

        for i in range(3):
            gy = mt + ch - int((i / 2) * ch)
            p.setPen(QPen(QColor(COLORS["border"]), 1, Qt.PenStyle.DotLine))
            p.drawLine(ml, gy, w - mr, gy)
            p.setPen(QColor(COLORS["text2"]))
            p.setFont(QFont("Consolas", 7))
            p.drawText(w - mr - 2, gy + 4, STATE_NAMES[i][0])

        if n > 1:
            for k in range(n - 1):
                x1 = ml + int(k * step)
                y1 = mt + ch - int(self.history[k] / 2 * ch)
                x2 = ml + int((k + 1) * step)
                y2 = mt + ch - int(self.history[k + 1] / 2 * ch)
                col = QColor(STATE_COLORS[self.history[k]])
                col.setAlpha(200)
                p.setPen(QPen(col, 2))
                p.drawLine(x1, y1, x2, y2)

        lx = ml + int((n - 1) * step)
        ly = mt + ch - int(self.history[-1] / 2 * ch)
        col = QColor(STATE_COLORS[self.history[-1]])
        p.setBrush(QBrush(col))
        p.setPen(QPen(QColor(COLORS["surface"]), 2))
        p.drawEllipse(lx - 4, ly - 4, 8, 8)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Марковская модель погоды")
        self.setMinimumSize(1100, 760)

        self.trans_matrix = [
            [0.6, 0.3, 0.1],
            [0.2, 0.5, 0.3],
            [0.1, 0.3, 0.6],
        ]
        self.current_state = 0
        self.day = 0
        self.running = False
        self.counts = [0, 0, 0]
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
        root.setSpacing(10)
        root.setContentsMargins(10, 10, 10, 10)

        left = QVBoxLayout()
        left.setSpacing(8)

        title = QLabel("МАРКОВСКАЯ МОДЕЛЬ ПОГОДЫ")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet(f"""
            color: {COLORS['accent']};
            font-family: Consolas;
            font-size: 14px;
            font-weight: bold;
            letter-spacing: 2px;
            padding: 8px;
            border-bottom: 1px solid {COLORS['border']};
        """)
        left.addWidget(title)

        self.markov_canvas = MarkovCanvas()
        left.addWidget(self.markov_canvas)

        self.matrix_editor = MatrixEditor()
        self.matrix_editor.matrix_changed.connect(self._on_matrix_changed)
        left.addWidget(self.matrix_editor)

        norm_btn = QPushButton("⚖ Нормализовать строки")
        norm_btn.clicked.connect(self.matrix_editor.normalize)
        norm_btn.setStyleSheet(self._btn_style(COLORS["surface2"], COLORS["accent"]))
        left.addWidget(norm_btn)

        root.addLayout(left, 38)

        center = QVBoxLayout()
        center.setSpacing(8)

        self.weather_widget = WeatherWidget()
        self.weather_widget.setFixedWidth(280)
        weather_wrap = QHBoxLayout()
        weather_wrap.addStretch()
        weather_wrap.addWidget(self.weather_widget)
        weather_wrap.addStretch()
        center.addLayout(weather_wrap)

        ctrl = QHBoxLayout()
        self.start_btn = QPushButton("▶  СТАРТ")
        self.start_btn.clicked.connect(self._toggle_run)
        self.start_btn.setStyleSheet(self._btn_style(COLORS["accent3"], COLORS["bg"], bold=True))
        self.step_btn = QPushButton("⏭  ШАГ")
        self.step_btn.clicked.connect(self._manual_step)
        self.step_btn.setStyleSheet(self._btn_style(COLORS["surface2"], COLORS["accent"]))
        self.reset_btn = QPushButton("↺  СБРОС")
        self.reset_btn.clicked.connect(self._reset)
        self.reset_btn.setStyleSheet(self._btn_style(COLORS["surface2"], COLORS["accent2"]))
        self.save_btn = QPushButton("💾  CSV")
        self.save_btn.clicked.connect(self._save_csv)
        self.save_btn.setStyleSheet(self._btn_style(COLORS["surface2"], COLORS["accent3"]))
        for b in [self.start_btn, self.step_btn, self.reset_btn, self.save_btn]:
            ctrl.addWidget(b)
        center.addLayout(ctrl)

        spd_row = QHBoxLayout()
        spd_lbl = QLabel("Скорость:")
        spd_lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:10px;")
        self.speed_slider = QSlider(Qt.Orientation.Horizontal)
        self.speed_slider.setRange(50, 2000)
        self.speed_slider.setValue(800)
        self.speed_slider.setInvertedAppearance(True)
        self.speed_slider.valueChanged.connect(self._on_speed)
        self.speed_slider.setStyleSheet(f"""
            QSlider::groove:horizontal {{ background:{COLORS['border']}; height:4px; border-radius:2px; }}
            QSlider::handle:horizontal {{
                background:{COLORS['accent']}; width:14px; height:14px;
                margin:-5px 0; border-radius:7px;
            }}
            QSlider::sub-page:horizontal {{ background:{COLORS['accent']}; border-radius:2px; }}
        """)
        self.speed_val_lbl = QLabel("0.8 с/шаг")
        self.speed_val_lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:10px; width:60px;")
        spd_row.addWidget(spd_lbl)
        spd_row.addWidget(self.speed_slider)
        spd_row.addWidget(self.speed_val_lbl)
        center.addLayout(spd_row)

        self.day_lbl = QLabel("День: 0  |  Шагов: 0")
        self.day_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.day_lbl.setStyleSheet(f"color:{COLORS['accent']}; font-family:Consolas; font-size:11px;")
        center.addWidget(self.day_lbl)

        root.addLayout(center, 30)

        right = QVBoxLayout()
        right.setSpacing(8)

        self.dist_chart = DistributionChart()
        right.addWidget(self.dist_chart)

        self.hist_chart = HistoryChart()
        right.addWidget(self.hist_chart)

        right.setSpacing(2)
        right.setContentsMargins(0, 0, 0, 0)

        stats_lbl = QLabel("СТАТИСТИКА")
        stats_lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:9px;")
        right.addWidget(stats_lbl)

        self.stats_table = StatsTable()

        self.stats_table.verticalHeader().setDefaultSectionSize(18)
        self.stats_table.horizontalHeader().setFixedHeight(20)

        self.stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.stats_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)

        self.stats_table.setFixedHeight(96)

        self.stats_table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.stats_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        right.addWidget(self.stats_table)


        trans_lbl = QLabel("МАТРИЦА ЧАСТОТ ПЕРЕХОДОВ")
        trans_lbl.setStyleSheet(f"color:{COLORS['text2']}; font-family:Consolas; font-size:9px;")
        right.addWidget(trans_lbl)

        self.trans_count_table = QTableWidget(3, 4, self)
        self.trans_count_table.setHorizontalHeaderLabels(["", "→Ясно", "→Обл.", "→Пасм."])

        self.trans_count_table.verticalHeader().setVisible(False)
        self.trans_count_table.verticalHeader().setDefaultSectionSize(18)
        self.trans_count_table.horizontalHeader().setFixedHeight(20)

        self.trans_count_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.trans_count_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)

        self.trans_count_table.setFixedHeight(78)

        self.trans_count_table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.trans_count_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        self.trans_count_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.trans_count_table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self.trans_count_table.setStyleSheet(self.stats_table.styleSheet() + "font-size: 9px; padding: 0px;")

        right.addWidget(self.trans_count_table)

        self.trans_counts = [[0] * 3 for _ in range(3)]

        root.addLayout(right, 32)

        self.markov_canvas.set_matrix(self.trans_matrix)
        self.markov_canvas.set_state(self.current_state)
        self.weather_widget.set_state(self.current_state, self.day, [self.current_state])
        self._update_stats()

    def _btn_style(self, bg, border, bold=False):
        w = "bold" if bold else "normal"
        return f"""
            QPushButton {{
                background: {bg};
                color: {border};
                border: 1px solid {border};
                border-radius: 6px;
                padding: 7px 12px;
                font-family: Consolas;
                font-size: 10px;
                font-weight: {w};
            }}
            QPushButton:hover {{ background: {border}22; }}
            QPushButton:pressed {{ background: {border}44; }}
        """

    def _apply_styles(self):
        self.setStyleSheet(f"""
            QMainWindow, QWidget {{ background: {COLORS['bg']}; color: {COLORS['text']}; }}
            QSplitter::handle {{ background: {COLORS['border']}; }}
        """)

    def _update_theoretical(self):
        self.theoretical = stationary_distribution(self.trans_matrix).tolist()

    def _on_matrix_changed(self, m):
        self.trans_matrix = m
        self.markov_canvas.set_matrix(m)
        self._update_theoretical()
        self._update_stats()

    def _toggle_run(self):
        self.running = not self.running
        if self.running:
            self.start_btn.setText("⏸  ПАУЗА")
            self.start_btn.setStyleSheet(self._btn_style(COLORS["accent2"], COLORS["bg"], bold=True))
            self._step_timer.start(self.speed_ms)
        else:
            self.start_btn.setText("▶  СТАРТ")
            self.start_btn.setStyleSheet(self._btn_style(COLORS["accent3"], COLORS["bg"], bold=True))
            self._step_timer.stop()

    def _manual_step(self):
        if not self.running:
            self._step()

    def _step(self):
        prev_state = self.current_state
        new_state = next_state(self.current_state, self.trans_matrix)

        self.trans_counts[prev_state][new_state] += 1

        self.current_state = new_state
        self.counts[new_state] += 1
        self.history.append(new_state)
        self.day += 1

        self.sim_data.append({
            "day": self.day,
            "state": STATE_IDS[new_state],
            "state_name": STATE_NAMES[new_state],
            "from_state": STATE_IDS[prev_state],
            "timestamp": datetime.now().isoformat()
        })

        self.markov_canvas.set_state(new_state, prev_state)
        self.weather_widget.set_state(new_state, self.day - 1, self.history)
        self.hist_chart.add_state(new_state)
        self.day_lbl.setText(f"День: {self.day}  |  Шагов: {sum(self.counts)}")
        self._update_stats()

    def _update_stats(self):
        total = sum(self.counts)
        emp = [c / total if total > 0 else 1/3 for c in self.counts]
        self.dist_chart.update_data(emp, self.theoretical)
        self.stats_table.update_stats(self.counts, max(total, 1), self.theoretical)

        for i in range(3):
            row_total = max(sum(self.trans_counts[i]), 1)
            name_item = QTableWidgetItem(f"{STATE_EMOJI[i]} {STATE_NAMES[i]}")
            name_item.setForeground(QColor(STATE_COLORS[i]))
            self.trans_count_table.setItem(i, 0, name_item)
            for j in range(3):
                freq = self.trans_counts[i][j] / row_total
                item = QTableWidgetItem(f"{self.trans_counts[i][j]} ({freq:.2f})")
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.trans_count_table.setItem(i, j + 1, item)

    def _reset(self):
        was_running = self.running
        if self.running:
            self._toggle_run()
        self.current_state = 0
        self.day = 0
        self.counts = [0, 0, 0]
        self.history = []
        self.sim_data = []
        self.trans_counts = [[0]*3 for _ in range(3)]
        self.markov_canvas.set_state(0)
        self.weather_widget.set_state(0, 0, [0])
        self.hist_chart.history.clear()
        self.hist_chart.update()
        self.day_lbl.setText("День: 0  |  Шагов: 0")
        self._update_stats()

    def _on_speed(self, v):
        self.speed_ms = v
        self.speed_val_lbl.setText(f"{v/1000:.1f} с/шаг")
        if self.running:
            self._step_timer.setInterval(v)

    def _save_csv(self):
        if not self.sim_data:
            QMessageBox.warning(self, "Нет данных", "Запустите симуляцию перед сохранением.")
            return
        fname = f"weather_simulation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        total = sum(self.counts)
        emp = [c / total if total > 0 else 1/3 for c in self.counts]

        try:
            with open(fname, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["day","state","state_name","from_state","timestamp"])
                writer.writeheader()
                writer.writerows(self.sim_data)

            QMessageBox.information(self, "Сохранено", f"Файл сохранён:\n{fname}")
        except Exception as e:
            QMessageBox.critical(self, "Ошибка", str(e))

    def on_speed(self, v):
        self._on_speed(v)



def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window,        QColor(COLORS["bg"]))
    palette.setColor(QPalette.ColorRole.WindowText,    QColor(COLORS["text"]))
    palette.setColor(QPalette.ColorRole.Base,          QColor(COLORS["surface"]))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(COLORS["surface2"]))
    palette.setColor(QPalette.ColorRole.Button,        QColor(COLORS["surface2"]))
    palette.setColor(QPalette.ColorRole.ButtonText,    QColor(COLORS["text"]))
    palette.setColor(QPalette.ColorRole.Highlight,     QColor(COLORS["accent"]))
    app.setPalette(palette)

    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()