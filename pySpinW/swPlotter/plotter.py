import pyvista
import vtk
# This will be removed when I'm sensible
from PyQt5 import Qt, QtWidgets
from PyQt5.Qt import QMainWindow
from .baseplotter import plotter
from .plotterTools import plotterTools

class showPlotter:
    def __init__(self, spinWobj, interface=None, **kwargs):

        if interface is None:
            self.interface = pyvista.Plotter()
        else:
            self.interface = interface

        self.plotter = plotter(spinWobj, self.interface, **kwargs)
        self.interface.camera.ParallelProjectionOn()

        self.logo = plotterTools.textBase(self.interface)
        self.logo.textColor = 'white'
        self.logo.font_size = 12
        self.logo.addText('PySpinW - Unreleased', 0.95, 0.95, 'rt')

    def setextent(self, extent):
        self.plotter.extent = extent
        if hasattr(self, 'highlighter'):
            self.highlighter.addActors(self.plotter.componentsShow)

    def add_lattice(self):
        self.plotter.add_lattice()

    def rm_lattice(self):
        self.plotter.rm_lattice()

    def add_atoms(self):
        self.plotter.add_atoms()

    def rm_atoms(self):
        self.plotter.rm_atoms()

    def add_bonds(self):
        self.plotter.add_bonds()

    def rm_bonds(self):
        self.plotter.rm_bonds()

    def along_a(self):
        self.plotter.along_a()

    def along_b(self):
        self.plotter.along_b()

    def along_c(self):
        self.plotter.along_c()

    def rotate(self, angle1, angle2):
        self.plotter.rotate(angle1, angle2)

    def removeLogo(self):
        self.logo.removeText()


class winPlotter(showPlotter):
    def __init__(self, spinWobj, interface=None, **kwargs):
        showPlotter.__init__(self, spinWobj, interface, **kwargs)
        self.highlighter = plotterTools.MouseClicker(self, ['highlightClick', 'shiftLeftClick'])
        self.highlighter.addActors(self.plotter.componentsShow)
        self.highlighter.SetDefaultRenderer(self.interface.renderer)
        self.interface.iren.SetInteractorStyle(self.highlighter)
        self.logo.textColor = 'white'

    def plot(self):
        self.highlighter.addActors(self.plotter.componentsShow)
        #Do we have a plotting window?
        if hasattr(self.interface, 'ren_win'):
            self.interface.show(interactive=True)
        else:
            self.interface.plot(interactive=True)

    def setextent(self, extent):
        showPlotter.setextent(self, extent)


class nbPlotter(showPlotter):
    def __init__(self, spinWobj, interface=pyvista.Plotter(off_screen=True), **kwargs):
        import panel as pn
        showPlotter.__init__(self, spinWobj, interface, **kwargs)
        pn.extension('vtk')
        self.VTK = pn.pane.VTK
        self.removeLogo()
        self.__plotsize = [300, 300]

    def setSize(self, width=300, height=300):
        self.__plotsize = [width, height]

    def plot(self):
        return self.VTK(self.interface.ren_win, width=self.__plotsize[0], height=self.__plotsize[1])


class qtWindow(showPlotter, QMainWindow):
    def __init__(self, spinWobj, parent=None, show=True):

        # Create a window
        QMainWindow.__init__(self, parent)

        # create the frame
        self.frame = Qt.QFrame()
        # Create the plotting interface
        showPlotter.__init__(self, spinWobj, pyvista.QtInteractor(self.frame))

        vlayout = Qt.QVBoxLayout()
        # Button Layout
        btn_layout = Qt.QHBoxLayout()
        btn_layout.addItem(Qt.QSpacerItem(0, 0, Qt.QSizePolicy.Expanding, Qt.QSizePolicy.Minimum))
        btn1 = Qt.QPushButton("Button 1")
        btn1.clicked.connect(self.along_c)
        btn2 = Qt.QPushButton("Extend")
        btn2.clicked.connect(self.extentDlg)
        btn_layout.addWidget(btn1)
        btn_layout.addWidget(btn2)

        # add the vtki interactor object
        vlayout.addWidget(self.interface)
        vlayout.addLayout(btn_layout)

        self.frame.setLayout(vlayout)
        self.setCentralWidget(self.frame)

        # simple menu to demo functions
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        exitButton = Qt.QAction('Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)

        self.interface.camera.ParallelProjectionOn()

        # Add elements menu
        meshMenu = mainMenu.addMenu('Add Elements')
        self.add_lattice_action = Qt.QAction('Lattice', self)
        self.add_lattice_action.triggered.connect(self.add_lattice)
        meshMenu.addAction(self.add_lattice_action)
        self.add_atoms_action = Qt.QAction('Atoms', self)
        self.add_atoms_action.triggered.connect(self.add_atoms)
        meshMenu.addAction(self.add_atoms_action)
        self.add_bonds_action = Qt.QAction('Bonds', self)
        self.add_bonds_action.triggered.connect(self.add_bonds)
        meshMenu.addAction(self.add_bonds_action)

        # Remove elements menu
        RMmeshMenu = mainMenu.addMenu('Remove Elements')
        self.rm_lattice_action = Qt.QAction('Lattice', self)
        self.rm_lattice_action.triggered.connect(self.rm_lattice)
        RMmeshMenu.addAction(self.rm_lattice_action)
        self.rm_atoms_action = Qt.QAction('Atoms', self)
        self.rm_atoms_action.triggered.connect(self.rm_atoms)
        RMmeshMenu.addAction(self.rm_atoms_action)
        self.rm_bonds_action = Qt.QAction('Bonds', self)
        self.rm_bonds_action.triggered.connect(self.rm_bonds)
        RMmeshMenu.addAction(self.rm_bonds_action)

        # add the mouse interactor style
        self.highlighter = plotterTools.MouseClicker(self, ['highlightClick', 'shiftLeftClick'])
        self.highlighter.addActors(self.plotter.componentsShow)
        self.highlighter.SetDefaultRenderer(self.interface.renderer)
        self.highlighter.tooltipRenderer.font_size = 12
        self.interface.iren.SetInteractorStyle(self.highlighter)

        self.interface.background_color = 'white'
        self.logo.textColor = 'black'

        if show:
            self.show()

    def extentDlg(self):
        dlg = extentDlg(self.plotter.extent)
        if dlg.exec_():
            self.setextent([int(float(dlg.extentA.text())), int(float(dlg.extentB.text())), int(float(dlg.extentC.text()))])


class extentDlg(QtWidgets.QDialog):

    def __init__(self, currExt):
        super(extentDlg, self).__init__()
        self.extentA = QtWidgets.QLineEdit()
        self.extentA.setText('{:0.02f}'.format(currExt[0]))
        self.extentB = QtWidgets.QLineEdit()
        self.extentB.setText('{:0.02f}'.format(currExt[1]))
        self.extentC = QtWidgets.QLineEdit()
        self.extentC.setText('{:0.02f}'.format(currExt[2]))
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        layout = QtWidgets.QFormLayout()
        layout.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        layout.addRow('Extent A', self.extentA)
        layout.addRow('Extent B', self.extentB)
        layout.addRow('Extent C', self.extentC)
        layout.addWidget(self.button_box)

        self.setLayout(layout)
        self.setWindowTitle("Set Extent")
        self.setMinimumWidth(200)