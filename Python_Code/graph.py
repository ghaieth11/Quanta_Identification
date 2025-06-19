# === imports===
import sys
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QHBoxLayout, 
    QWidget, QLabel, QStackedWidget, QDoubleSpinBox, QLineEdit, QComboBox
)
from PyQt5.QtWidgets import QVBoxLayout, QWidget
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl

from PyQt5.QtCore import QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import qutip as qt
from PyQt5.QtCore import Qt
import random

#===  pdf  =====
import os
import webbrowser


# ===modules ===
from system import System
from qubit import Qubit
from control import Control
from estimator import Estimator
from utils import BlochStates



# === Bloch sphere canvas class ===
class BlochSphereCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure()
        self.ax = fig.add_subplot(111, projection='3d')
        super().__init__(fig)
        self.setParent(parent)
        self._init_sphere()
        self.trajectory = []

    def _init_sphere(self):
        self.ax.clear()
        b = qt.Bloch(fig=self.figure, axes=self.ax)

        b.vector_color = ['blue']
        b.add_annotation([1, 0, 0], r'|+⟩')
        b.add_annotation([0, 1, 0], r'|i⟩')
        b.add_annotation([0, 0, 1], r'Z')
        b.render()
    
        
        self.draw()

    def update_vector(self, vec):
        self._init_sphere()
        self.trajectory.append(vec)
        if len(self.trajectory) > 1:
            traj = np.array(self.trajectory)
            self.ax.plot3D(traj[:, 1], -traj[:, 0], traj[:, 2], color='red', linewidth=1)
        self.ax.quiver(0, 0, 0, vec[1],-vec[0],vec[2], color='blue', linewidth=2)
        self.draw()

    def reset_trajectory(self):
        self.trajectory = []

# === Main GUI class ===
class QubitInterface(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Qubit Bloch Sphere Simulator")
        self.setMinimumSize(900, 700)

        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)

        self.init_welcome_page()
        self.init_setup_page()
        self.init_simulation_page()
        self.init_estimation_page()
        self.init_doc_page()

        self.qubit = None
        self.system = None
        self.control = None
        self.estimator = None
        self.t = 0.0

    # === Pages will be implemented in detail next ===
    
   #--------------Page 1 -------------------------- 
    def init_welcome_page(self):
        self.welcome_page = QWidget()

        # Main layout for the welcome page (vertical)
        main_layout = QVBoxLayout()
        main_layout.setSpacing(40)  # Space between title and buttons
        main_layout.setContentsMargins(50, 50, 50, 50)

        # Center title
        title = QLabel("Welcome to the Qubit Simulator")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("font-size: 28px; font-weight: bold;")
        main_layout.addWidget(title)

        # Add a spacer above and below to center vertically
        main_layout.addStretch()

        # Horizontal layout for the three buttons
        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(30)

        
        learn_btn = QPushButton("Learn")
        start_btn = QPushButton("Start")
        doc_btn = QPushButton("Documentation")

        for btn in [learn_btn, start_btn, doc_btn]:
            btn.setFixedSize(140, 50)
            btn.setStyleSheet("font-size: 16px;")
            btn_layout.addWidget(btn, alignment=Qt.AlignCenter)

        # Connect buttons to respective pages
        learn_btn.clicked.connect(lambda : self.open_pdf("Identification.pdf"))
        start_btn.clicked.connect(lambda: self.stack.setCurrentWidget(self.setup_page))
        doc_btn.clicked.connect(lambda: self.open_pdf("functional_doc.pdf"))

        # Add the horizontal button layout to the main layout
        main_layout.addLayout(btn_layout)
        main_layout.addStretch()

        self.welcome_page.setLayout(main_layout)
        self.stack.addWidget(self.welcome_page)
        
# --------------Page 2 ----------------------------------------------        
    
    def init_setup_page(self):
        self.setup_page = QWidget()
        
        layout = QVBoxLayout()
        
        self.canvas = BlochSphereCanvas(self)
        
        layout.addWidget(self.canvas)
        param_layout = QHBoxLayout()

        self.gamma1_input = QDoubleSpinBox()
        self.gamma1_input.setPrefix("γ₁: ")
        self.gamma1_input.setRange(0.0, 10.0)

        self.gamma2_input = QDoubleSpinBox()
        self.gamma2_input.setPrefix("γ₂: ")
        self.gamma2_input.setRange(0.0, 10.0)

        self.kappa_input = QDoubleSpinBox()
        self.kappa_input.setPrefix("κ: ")
        self.kappa_input.setRange(0.0, 10.0)

        self.omega_input = QDoubleSpinBox()
        self.omega_input.setPrefix("ω: ")
        self.omega_input.setRange(0.0, 10.0)

        self.control_input = QDoubleSpinBox()
        self.control_input.setPrefix("u: ")
        self.control_input.setRange(0.0, 10000.0)

        self.time_input = QDoubleSpinBox()
        self.time_input.setPrefix("t: ")
        self.time_input.setRange(0.0, 1000.0)

        for spinbox in [
        self.gamma1_input, self.gamma2_input, self.kappa_input,
        self.omega_input, self.control_input, self.time_input
            ]:
            spinbox.setSingleStep(0.1)
            param_layout.addWidget(spinbox)
        
        # Choix de l'état initial
        self.initial_state_combo = QComboBox()
        self.initial_state_combo.addItems(BlochStates.keys())
        self.initial_state_combo.setToolTip("Choisissez l'état initial du qubit")
        layout.addWidget(QLabel("État initial :"))
        layout.addWidget(self.initial_state_combo)

        layout.addLayout(param_layout)
        
        random_btn = QPushButton("🎲 Randomize Parameters")
        random_btn.clicked.connect(self.randomize_parameters)
        layout.addWidget(random_btn)
        
        button_layout = QHBoxLayout()

        back_btn = QPushButton("← Back to Menu")
        back_btn.clicked.connect(lambda: self.stack.setCurrentWidget(self.welcome_page))
        button_layout.addWidget(back_btn)

        run_btn = QPushButton("▶ Run Simulation")
        run_btn.clicked.connect(self.run_simulation)  # à définir
        button_layout.addWidget(run_btn)

        layout.addLayout(button_layout)
        
        self.setup_page.setLayout(layout)
        self.stack.addWidget(self.setup_page)
        
# ----------------- Page 3: Simulation Page -----------------
    def init_simulation_page(self):
        self.simulation_page = QWidget()
        layout = QVBoxLayout()

        self.canvas_result = BlochSphereCanvas(self)
        layout.addWidget(self.canvas_result)

        # Boutons de navigation
        button_layout = QHBoxLayout()

        reset_btn = QPushButton("← Reset")
        reset_btn.clicked.connect(lambda: self.stack.setCurrentWidget(self.setup_page))
        button_layout.addWidget(reset_btn)

        estimate_btn = QPushButton("▶ Estimate")
        estimate_btn.clicked.connect(lambda: self.stack.setCurrentWidget(self.estimate_page))
        button_layout.addWidget(estimate_btn)

        layout.addLayout(button_layout)

        self.simulation_page.setLayout(layout)
        self.stack.addWidget(self.simulation_page)
        
        
    def run_simulation(self) : 
        self.stack.setCurrentWidget(self.simulation_page)
        
        # 1. Récupération des paramètres depuis les spinboxes de la page 2
        gamma1 = self.gamma1_input.value()
        gamma2 = self.gamma2_input.value()
        kappa = self.kappa_input.value()
        omega = self.omega_input.value()
        u = self.control_input.value()
        t_final = self.time_input.value()
        
        # Lire le choix de l'état initial
        selected_label = self.initial_state_combo.currentText()
        initial_vector = BlochStates[selected_label]
        
        # 2. Création des objets physiques
        self.qubit = Qubit(omega,kappa,gamma1,gamma2)
        self.system = System(*initial_vector)
        self.control = Control(u)
        self.estimator = Estimator(self.system)
        
        # 3. Initialisation du timer de simulation
        self.sim_duration = t_final
        self.t = 0.0
        self.canvas_result.reset_trajectory()  # canvas de la page 3
        self.timer = QTimer()
        self.timer.setInterval(10)
        self.timer.timeout.connect(self.step_simulation)
        self.timer.start()
    
    def step_simulation(self):
        dt = 0.01
        self.t += dt

        self.system.evolve(self.qubit, self.control, dt)
        x, y, z = self.system.get_coordinates() 
        self.canvas_result.update_vector([x, y, z])

        if self.t >= self.sim_duration:
            self.timer.stop()
                    
    
# ------------------ page 4 -------------------------------------------    
    def init_estimation_page(self):
        
        gamma1 = self.gamma1_input.value()
        gamma2 = self.gamma2_input.value()
        kappa = self.kappa_input.value()
        omega = self.omega_input.value()
        u = self.control_input.value()
        t_final = self.time_input.value()
        
        # 2. Création des objets physiques
        self.qubit = Qubit(omega,kappa,gamma1,gamma2)
        self.control = Control(u)
        
        self.estimate_page = QWidget()
        layout = QVBoxLayout()

        self.canvas_estimation_page = BlochSphereCanvas(self)
        layout.addWidget(self.canvas_estimation_page)

        # Boutons de navigation
        button_layout = QHBoxLayout()
        
        reset_btn = QPushButton("← Go Back")
        reset_btn.clicked.connect(lambda: self.stack.setCurrentWidget(self.setup_page))
        button_layout.addWidget(reset_btn)

        reset_btn = QPushButton("Estimate ω")
        reset_btn.clicked.connect(self.estimate_omega)
        button_layout.addWidget(reset_btn)
        
        reset_btn = QPushButton("Estimate κ")
        reset_btn.clicked.connect(self.estimate_kappa)
        button_layout.addWidget(reset_btn)
        
        reset_btn = QPushButton("Estimate γ₁")
        reset_btn.clicked.connect(self.estimate_gamma1)
        button_layout.addWidget(reset_btn)
        
        reset_btn = QPushButton("Estimate γ₂")
        reset_btn.clicked.connect(self.estimate_gamma2)
        button_layout.addWidget(reset_btn)
        

        layout.addLayout(button_layout)
        
        self.label = QLabel("Résultat de l’estimation apparaîtra ici")
        layout.addWidget(self.label)

        self.estimate_page.setLayout(layout)
        self.stack.addWidget(self.estimate_page)
        
        
    def estimate_gamma1(self):
        S = System(0, 0, -1)
        E = Estimator(S)
        gamma1_est = E.estimate_gamma1(self.qubit, n=100000)
        true_val = self.qubit.get_gamma1()
        self.label.setText(f"γ₁: true={true_val:.4f}, estimated={gamma1_est:.4f}")
        self.animate_trajectory(S, self.qubit, self.control, t_final=self.time_input.value())
        

    def estimate_gamma2(self):
        S = System(0, 1, 0)
        E = Estimator(S)
        gamma2_est = E.estimate_gamma2(self.qubit, n=100000)
        true_val = self.qubit.get_gamma2()
        self.label.setText(f"γ₂: true={true_val:.4f}, estimated={gamma2_est:.4f}")
        self.animate_trajectory(S, self.qubit, self.control, t_final=self.time_input.value())
        

    def estimate_kappa(self):
        S = System(0, 0, -1)
        S.initialize()  # Nécessaire pour lancer le protocole
        E = Estimator(S)
        kappa_est = E.estimate_kappa(self.qubit, n=100000)
        true_val = self.qubit.get_kappa()
        self.label.setText(f"κ: true={true_val:.4f}, estimated={kappa_est:.4f}")
        self.animate_trajectory(S, self.qubit, self.control, t_final=self.time_input.value())
        

    def estimate_omega(self):
        S = System(0, 1, 0)
        S.initialize()
        E = Estimator(S)
        omega_est = E.estimate_omega(self.qubit, n=100000)
        true_val = self.qubit.get_omega()
        self.label.setText(f"ω: true={true_val:.4f}, estimated={omega_est:.4f}")
        self.animate_trajectory(S, self.qubit, self.control, t_final=self.time_input.value())
        
        
    
    def animate_trajectory(self, system, qubit, control, t_final):
        if self.stack.currentWidget() == self.estimate_page:
            canvas = self.canvas_estimation_page
        else:
            canvas = self.canvas_result
        
        canvas.reset_trajectory()
        self.t = 0.0
        self.sim_duration = t_final
        self.temp_system = system
        self.temp_qubit = qubit
        self.temp_control = control
        self.temp_canvas = canvas

        self.timer_estimation = QTimer()
        self.timer_estimation.setInterval(10)
        self.timer_estimation.timeout.connect(self.step_estimation_animation)
        self.timer_estimation.start()
    
    
    def step_estimation_animation(self):
        dt = 0.01
        self.t += dt

        self.temp_system.evolve(self.temp_qubit, self.temp_control, dt)
        x, y, z = self.temp_system.get_coordinates()
        self.temp_canvas.update_vector([x, y, z])

        if self.t >= self.sim_duration:
            self.timer_estimation.stop()
        
        
    
    def randomize_parameters(self):
        self.qubit = Qubit.random()

        self.gamma1_input.setValue(self.qubit.get_gamma1())
        self.gamma2_input.setValue(self.qubit.get_gamma2())
        self.kappa_input.setValue(self.qubit.get_kappa())
        self.omega_input.setValue(self.qubit.get_omega())
        self.control_input.setValue(random.uniform(0, 10))
        self.time_input.setValue(random.uniform(0.5, 10.0))
    
    
    def init_doc_page(self):
        self.doc_page = QWidget()
        layout = QVBoxLayout()

        label = QLabel("Documentation Viewer")
        label.setStyleSheet("font-size: 20px; font-weight: bold; margin-bottom: 10px;")

        layout.addWidget(label)
        self.doc_page.setLayout(layout)
        self.stack.addWidget(self.doc_page)

    def open_pdf(self, pdf_path):
        if not os.path.exists(pdf_path):
            print(f"PDF file not found: {pdf_path}")
            return

        if pdf_path == "Identification.pdf":
            print("Opening Identification document...")
        elif pdf_path == "functional_doc.pdf":
            print("Opening functional documentation...")
        elif pdf_path == "technical_doc.pdf":
            print("Opening technical documentation...")
        else:
            print(f"Opening unknown PDF file: {pdf_path}")

        os.system(f"open '{pdf_path}'")  # Use appropriate command for your OS
    
    
    
if __name__ == "__main__": 
    app = QApplication(sys.argv)
    window = QubitInterface()
    window.show()
    sys.exit(app.exec_())