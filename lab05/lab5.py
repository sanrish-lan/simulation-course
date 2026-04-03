import sys
import random
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QPushButton, QSpinBox,
                             QGridLayout, QGroupBox, QFrame, QTabWidget)
from PyQt6.QtCore import QTimer, Qt, QUrl, QRect, QEvent
from PyQt6.QtGui import QFont, QColor, QLinearGradient, QPainter, QBrush, QPalette
from PyQt6.QtMultimedia import QMediaPlayer, QAudioOutput

SYMBOLS = ["🍒", "🍋", "🍊", "🍇", "💎", "🎰"]
PAYOUTS = {"🍒": 5, "🍋": 10, "🍊": 15, "🍇": 20, "💎": 50, "🎰": 100}

STYLE = """
QMainWindow { background-color: #0d0d0d }
QTabWidget::pane { border: 2px solid #5e4d1a; background: #121212; border-radius: 10px }
QTabBar::tab { background: #1a1a1a; color: #8a732e; padding: 10px 25px; border-top-left-radius: 8px; border-top-right-radius: 8px; margin-right: 2px }
QTabBar::tab:selected { background: #d4af37; color: black; font-weight: bold }

QPushButton#SpinButton {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #d4af37, stop:1 #8a732e);
    color: black; border-radius: 12px; font-size: 24px; font-weight: bold; border: 1px solid #ffecb3; min-height: 55px
}
QPushButton#SpinButton:hover { background: #f1c40f }
QPushButton#SpinButton:pressed { background: #aa8a2e }
QPushButton#SpinButton:disabled { background: #333; color: #666; border: 1px solid #222 }

QPushButton#DepButton { background: #27ae60; color: white; border-radius: 5px; font-weight: bold; padding: 5px }

QLabel#BalanceLabel { color: #f1c40f; font-size: 22px; font-weight: bold }
QLabel#BetLabel { color: #fff; font-size: 14px; font-weight: bold }

QGroupBox { color: #d4af37; font-weight: bold; border: 1px solid #3d3d3d; margin-top: 15px }
QSpinBox { background: #1a1a1a; color: #d4af37; border: 1px solid #5e4d1a; border-radius: 4px; padding: 4px; font-size: 14px }
"""


class Reel(QFrame):
    def __init__(self):
        super().__init__()
        self.setFixedSize(130, 210)
        self.setStyleSheet("background-color: #f5f5f5; border-radius: 10px; border: 3px solid #d4af37")
        self.font = QFont("Segoe UI Emoji", 80)
        self.current_sym = "🎰"
        self.next_sym = random.choice(SYMBOLS)
        self.offset = 0
        self.state = "stopped"
        self.final_sym = ""
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_animation)

    def set_final_symbol(self, sym):
        self.final_sym = sym
        self.state = "stopping"

    def update_animation(self):
        if self.state == "stopped": return
        self.offset += 35
        if self.offset >= 210:
            self.offset = 0
            self.current_sym = self.next_sym
            if self.state == "finishing":
                self.state = "stopped"
            elif self.state == "stopping":
                self.next_sym = self.final_sym
                self.state = "finishing"
            else:
                self.next_sym = random.choice(SYMBOLS)
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        painter.setFont(self.font)
        painter.drawText(QRect(0, self.offset, 130, 210), Qt.AlignmentFlag.AlignCenter, self.current_sym)
        painter.drawText(QRect(0, self.offset - 210, 130, 210), Qt.AlignmentFlag.AlignCenter, self.next_sym)
        painter.setPen(Qt.PenStyle.NoPen)

        top_grad = QLinearGradient(0, 0, 0, 80)
        top_grad.setColorAt(0.0, QColor(0, 0, 0, 230))
        top_grad.setColorAt(1.0, QColor(0, 0, 0, 0))
        painter.setBrush(QBrush(top_grad))
        painter.drawRect(0, 0, 130, 80)

        bot_grad = QLinearGradient(0, 130, 0, 210)
        bot_grad.setColorAt(0.0, QColor(0, 0, 0, 0))
        bot_grad.setColorAt(1.0, QColor(0, 0, 0, 230))
        painter.setBrush(QBrush(bot_grad))
        painter.drawRect(0, 130, 130, 80)

        glass = QLinearGradient(0, 0, 130, 0)
        glass.setColorAt(0, QColor(255, 255, 255, 0))
        glass.setColorAt(0.5, QColor(255, 255, 255, 30))
        glass.setColorAt(1, QColor(255, 255, 255, 0))
        painter.setBrush(QBrush(glass))
        painter.drawRect(0, 0, 130, 210)


class CasinoTab(QWidget):
    def __init__(self):
        super().__init__()
        self.balance = 5000
        layout = QVBoxLayout(self)

        balance_panel = QHBoxLayout()
        self.lbl_balance = QLabel(f"БАЛАНС: {self.balance} $")
        self.lbl_balance.setObjectName("BalanceLabel")
        btn_dep = QPushButton("ДЕПОЗИТ +1000")
        btn_dep.setObjectName("DepButton")
        btn_dep.clicked.connect(self.deposit)
        balance_panel.addWidget(self.lbl_balance)
        balance_panel.addStretch()
        balance_panel.addWidget(btn_dep)
        layout.addLayout(balance_panel)

        m_border = QFrame()
        m_border.setStyleSheet(
            "background-color: #1a0505; border: 5px solid #d4af37; border-radius: 20px; padding: 10px")
        m_layout = QVBoxLayout(m_border)
        title = QLabel("КОРОЛЕВСКИЙ КУШ")
        title.setFont(QFont("Impact", 32))
        title.setStyleSheet("color: #d4af37; border: none")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        m_layout.addWidget(title)

        reels_layout = QHBoxLayout()
        self.reels = [Reel(), Reel(), Reel()]
        for r in self.reels: reels_layout.addWidget(r)
        m_layout.addLayout(reels_layout)

        self.lbl_status = QLabel("ПРОБЕЛ / ENTER — ИГРАТЬ")
        self.lbl_status.setFont(QFont("Arial Black", 12))
        self.lbl_status.setStyleSheet("color: #fff; border: none")
        self.lbl_status.setAlignment(Qt.AlignmentFlag.AlignCenter)
        m_layout.addWidget(self.lbl_status)
        layout.addWidget(m_border)

        bet_panel = QHBoxLayout()
        self.spin_bet = QSpinBox()
        self.spin_bet.setRange(10, 100000000)
        self.spin_bet.setSingleStep(10)
        self.spin_bet.setValue(100)
        self.spin_bet.setButtonSymbols(QSpinBox.ButtonSymbols.NoButtons)
        self.spin_bet.setMinimumWidth(120)

        lbl_bet = QLabel("ВАША СТАВКА:")
        lbl_bet.setObjectName("BetLabel")
        bet_panel.addStretch()
        bet_panel.addWidget(lbl_bet)
        bet_panel.addWidget(self.spin_bet)
        bet_panel.addStretch()
        layout.addLayout(bet_panel)

        self.btn_spin = QPushButton("К Р У Т И Т Ь")
        self.btn_spin.setObjectName("SpinButton")
        self.btn_spin.clicked.connect(self.start_spin)
        layout.addWidget(self.btn_spin)

        w_box = QGroupBox("НАСТРОЙКА ВЕРОЯТНОСТЕЙ")
        w_grid = QGridLayout(w_box)
        self.inputs = []
        self.p_labels = []
        for i, s in enumerate(SYMBOLS):
            lb = QLabel(f"{s} (x{PAYOUTS[s]})")
            lb.setFont(QFont("Arial", 11))
            sb = QSpinBox()
            sb.setRange(1, 100)
            sb.setValue(1)
            # УБИРАЕМ СТРЕЛОЧКИ
            sb.setButtonSymbols(QSpinBox.ButtonSymbols.NoButtons)
            sb.valueChanged.connect(self.update_percents)
            pl = QLabel("0.0%")
            pl.setStyleSheet("color: #888; font-weight: bold")
            w_grid.addWidget(lb, i // 3, (i % 3) * 3)
            w_grid.addWidget(sb, i // 3, (i % 3) * 3 + 1)
            w_grid.addWidget(pl, i // 3, (i % 3) * 3 + 2)
            self.inputs.append(sb)
            self.p_labels.append(pl)
        layout.addWidget(w_box)
        self.update_percents()

        # Настройка звуков
        self.spin_player = QMediaPlayer()
        self.spin_audio = QAudioOutput()
        self.spin_player.setAudioOutput(self.spin_audio)
        self.spin_audio.setVolume(0.5)
        self.spin_player.setSource(QUrl.fromLocalFile(r"gambling.MP3"))
        self.spin_player.setLoops(QMediaPlayer.Loops.Infinite)

        self.jack_player = QMediaPlayer()
        self.jack_audio = QAudioOutput()
        self.jack_player.setAudioOutput(self.jack_audio)
        self.jack_audio.setVolume(0.8)
        self.jack_player.setSource(QUrl.fromLocalFile(r"moneyjackpot.mp3"))

        self.timer = QTimer()
        self.timer.timeout.connect(self.animate_tick)
        self.ticks = 0

    def update_percents(self):
        total = sum(sb.value() for sb in self.inputs)
        for i, sb in enumerate(self.inputs):
            p = (sb.value() / total) * 100
            self.p_labels[i].setText(f"{p:.1f}%")

    def deposit(self):
        self.balance += 1000
        self.lbl_balance.setText(f"БАЛАНС: {self.balance} $")

    def start_spin(self):
        if not self.btn_spin.isEnabled(): return
        bet = self.spin_bet.value()
        if self.balance < bet:
            self.lbl_status.setText("НЕДОСТАТОЧНО СРЕДСТВ!")
            return
        self.balance -= bet
        self.lbl_balance.setText(f"БАЛАНС: {self.balance} $")
        self.btn_spin.setEnabled(False)
        self.lbl_status.setText("ВРАЩЕНИЕ...")

        total = sum(sb.value() for sb in self.inputs)

        def calc():
            a, c = random.random(), 0
            for k, sb in enumerate(self.inputs):
                c += sb.value() / total
                if a < c: return SYMBOLS[k]
            return SYMBOLS[-1]

        self.res = [calc() for _ in range(3)]
        for r in self.reels:
            r.state = "spinning"
            r.timer.start(16)
        self.spin_player.play()
        self.jack_player.stop()
        self.ticks = 0
        self.timer.start(16)

    def animate_tick(self):
        self.ticks += 1
        # Остановка 1: 2.5с (156 тиков)
        if self.ticks == 156:
            self.reels[0].set_final_symbol(self.res[0])
        # Остановка 2: 4.0с (250 тиков)
        elif self.ticks == 250:
            self.reels[1].set_final_symbol(self.res[1])
        # Остановка 3: 5.5с (344 тика)
        elif self.ticks == 344:
            self.reels[2].set_final_symbol(self.res[2])

        if all(r.state == "stopped" for r in self.reels) and self.ticks > 344:
            self.timer.stop()
            self.spin_player.stop()
            self.final_check()

    def final_check(self):
        if self.res[0] == self.res[1] == self.res[2]:
            win = self.spin_bet.value() * PAYOUTS[self.res[0]]
            self.balance += win
            self.lbl_status.setText(f"ВЫИГРЫШ: {win} $ !!!")
            self.jack_player.play()
        else:
            self.lbl_status.setText("НЕ ПОВЕЗЛО...")
        self.lbl_balance.setText(f"БАЛАНС: {self.balance} $")
        self.btn_spin.setEnabled(True)


class YesNoTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        t = QLabel("ОРАКУЛ")
        t.setFont(QFont("Impact", 45))
        t.setStyleSheet("color: #d4af37")
        t.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(t)

        self.lbl_res = QLabel("?")
        self.lbl_res.setFont(QFont("Arial Black", 85))
        self.lbl_res.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.lbl_res)

        self.btn_ask = QPushButton("СПРОСИТЬ")
        self.btn_ask.setObjectName("SpinButton")
        self.btn_ask.clicked.connect(self.start_anim)
        layout.addWidget(self.btn_ask)

        self.timer = QTimer()
        self.timer.timeout.connect(self.step)
        self.ticks = 0

    def set_color(self, c):
        p = self.lbl_res.palette()
        p.setColor(QPalette.ColorRole.WindowText, QColor(c))
        self.lbl_res.setPalette(p)

    def start_anim(self):
        if not self.btn_ask.isEnabled(): return
        self.btn_ask.setEnabled(False)
        self.ticks = 0
        self.set_color("#444")
        self.timer.start(100)

    def step(self):
        self.ticks += 1
        if self.ticks < 15:
            self.lbl_res.setText(random.choice(["ДА", "НЕТ"]))
        else:
            self.timer.stop()
            res = "ДА" if random.random() > 0.5 else "НЕТ"
            self.lbl_res.setText(res)
            self.set_color("#2ecc71" if res == "ДА" else "#e74c3c")
            self.btn_ask.setEnabled(True)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Royal Casino V3")
        self.setFixedSize(720, 800)
        self.setStyleSheet(STYLE)
        self.tabs = QTabWidget()
        self.t1, self.t2 = YesNoTab(), CasinoTab()
        self.tabs.addTab(self.t1, "🔮 ОРАКУЛ")
        self.tabs.addTab(self.t2, "🎰 СЛОТЫ")
        self.setCentralWidget(self.tabs)
        # Фильтр событий для управления с клавиатуры
        QApplication.instance().installEventFilter(self)

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Type.KeyPress:
            if event.key() in (Qt.Key.Key_Space, Qt.Key.Key_Return, Qt.Key.Key_Enter):
                # Снимаем фокус с полей ввода, чтобы данные применились
                focused = QApplication.focusWidget()
                if isinstance(focused, QSpinBox):
                    focused.clearFocus()

                # Запуск действия в зависимости от вкладки
                if self.tabs.currentIndex() == 1:
                    self.t2.start_spin()
                else:
                    self.t1.start_anim()
                return True
        return super().eventFilter(obj, event)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())