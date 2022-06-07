from PyQt5 import QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
from pyqtgraph.dockarea import Dock
import scipy.signal.windows as windows
import scipy.fft as fft


from pyqtgraph.Qt import QtCore
import pyqtgraph as pg

from scipy.interpolate import interp1d
from pymodaq.daq_utils.daq_utils import linspace_step
from pymodaq.pid.utils import PIDModelGeneric, OutputToActuator, InputFromDetector
from pymodaq.daq_utils.plotting.viewer1D.viewer1D_main import Viewer1D
from pymodaq.daq_utils.daq_utils import Axis
from pymodaq.daq_utils.daq_utils import   DataFromPlugins
from pymodaq.daq_utils.math_utils import ft,ift



class PIDModelSpetralInterferometry(PIDModelGeneric):
    limits = dict(max=dict(state=False, value=1),
                  min=dict(state=False, value=-1),)
    konstants = dict(kp=1.0, ki=0.0, kd=0.0)

    Nsetpoint = 1
    setpoint_ini = [2. for ind in range(Nsetpoint)]
    # actuators = ['Move 01']
    actuators_name = ['Move 01']
    # detector = 'Spectro'
    detectors_name = ['CameraTIS']
    params = [
        {'title': 'Spectrum ROI', 'name': 'spectrum', 'type': 'group', 'expanded': True, 'visible': True,
         'children': [
            {'title': 'Omega min', 'name': 'omega_min', 'type': 'float', 'value':0.0 },
            {'title': 'Omega max', 'name': 'omega_max', 'type': 'float', 'value': 1.0}]},
        {'title': 'Inverse Fourier', 'name':'inverse_fourier','type': 'group', 'expanded': True, 'visible': True,
        'children': [
            {'title': 'N sampling (power of 2)', 'name': 'N_samp_power', 'type': 'float', 'value':13 },
            {'title': 'Centering', 'name': 'centering', 'type': 'bool', 'value':True },
            {'title': 'ROI', 'name': 'ifft_roi', 'type': 'group', 'expanded':True,'visible':True,
            'children':[
                {'title': 't min','name':'t_min','type' : 'float','value' : 0.0},
                {'title': 't max','name':'t_max','type' : 'float','value' : 1.0}]}]},
        {'title': 'Wavelength (nm)','name' : 'lambda_0','type':'float','value' : 633},
        {'title': 'Nbr of shots for RMS','name' : 'N_rms','type':'int','value' : 50},
        {'title': 'Show plots', 'name': 'show_plots', 'type': 'group', 'expanded': True, 'visible': True,
         'children': [
             {'title': 'Show camera data', 'name': 'show_camera', 'type': 'bool', 'value':True },
             {'title': 'Show ROI', 'name': 'show_roi', 'type': 'bool', 'value':True },
             {'title': 'Show FFT', 'name': 'show_fft', 'type': 'bool', 'value':True },
             {'title': 'Show phase', 'name': 'show_phase', 'type': 'bool', 'value':True }]}
         
    ]


    def __init__(self, pid_controller):
        super().__init__(pid_controller)


        self.d1 = Dock("Camera data")
        widget_calib = QtWidgets.QWidget()
        print('Init Widget: There should be a dock')
        self.pid_controller.dock_area.addDock(self.d1)
        self.img1a = pg.ImageItem()
        self.img1a.setImage(np.random.normal(size = (100,100)))
        w1 = pg.PlotWidget(title = "Camera Plot")
        w1.addItem(self.img1a)
        self.d1.addWidget(w1)

        self.roi = pg.RectROI([20,20],[20,20],pen=(5,30))
        self.roi.addRotateHandle([1,0],[0.5,0.5])
        w1.addItem(self.roi)

        
        # self.d2 = Dock("Camera data")
        # widget_calib = QtWidgets.QWidget()
        # print('Init Widget: There should be a dock')
        # self.pid_controller.dock_area.addDock(self.d2)
        # self.img1b = pg.ImageItem()
        # self.img1b.setImage(np.random.normal(size = (100,100)))
        # w2 = pg.PlotWidget(title = "Camera ROI")
        # w2.addItem(self.img1b)
        # self.d2.addWidget(w2)


        



        # w1 = pg.ImageView()
        # self.plotImage = w1.plot()
        # self.plotImage.setImage(np.random.normal(size=(100,100)))
        # self.d1.addWidget(w1)




        ##Plot des franges
        self.dock_calib = Dock('Fringes')
        widget_calib = QtWidgets.QWidget()
        print('Init Widget: There should be a dock')
        self.pid_controller.dock_area.addDock(self.dock_calib)
        self.u = pg.PlotWidget(title="Moyenne signal roi")
        self.plotfringes = self.u.plot()
        self.dock_calib.addWidget(self.u)
        self.plotfringes.setData(y = np.random.normal(size=100))
       

        ##Plot de la TF des franges
        self.dock_calib = Dock('Fourier transform')
        widget_calib = QtWidgets.QWidget()
        self.pid_controller.dock_area.addDock(self.dock_calib)
        self.r = pg.PlotWidget(title="Fourier transform ")
        self.plottf = self.r.plot()
        self.dock_calib.addWidget(self.r)
        self.plottf.setData(y = np.random.normal(size=100))
       


        #ROI sur la tf
        self.lr = pg.LinearRegionItem([0,1])
        self.lr.setZValue(-10)
        self.r.addItem(self.lr)



        ##Plot de la phase
        self.dock_calib = Dock('Phase')
        widget_calib = QtWidgets.QWidget()
        print('Init Widget: There should be a dock')
        self.pid_controller.dock_area.addDock(self.dock_calib)
        self.ph = pg.PlotWidget(title="Phase")
        self.plotph = self.ph.plot()
        self.dock_calib.addWidget(self.ph)
        test_signal = np.random.normal(size = 100)
        self.plotph.setData(y =np.random.normal(size = 100))

        self.input = np.array([0,0])
        # self.lbda = self.settings.child('lambda_0').value()
        self.lbda = 633 #nm
        self.c = 0.299792 #nm per as
        


        self.N_rms = self.settings.child('N_rms').value()
        self.delays = np.zeros(self.N_rms)
        self.i_loop = 0
        







    def updateGUI(self,param):
        self.viewer_calib.show_data(param)


    def update_settings(self, param):
        """
        Get a parameter instance whose value has been modified by a user on the UI
        Parameters
        ----------
        param: (Parameter) instance of Parameter object
        """
        if param.name() == 'show_camera':
            if  not param.value():
                self.d1.hide()
            else :
                self.d1.show()
        if param.name() == 'show_roi':
            if  not param.value():
                self.u.hide()
            else :
                self.u.show()
        if param.name() == 'show_fft':
            if  not param.value():
                self.r.hide()
            else :
                self.r.show()
        if param.name() == 'show_phase':
            if  not param.value():
                self.ph.hide()
            else :
                self.ph.show()

            

    def ini_model(self):
        super().ini_model()
    
    def convert_input(self, measurements):
        """
        Convert the image of the camera into x and y positions of the center of the beam.
        Parameters
        ----------
        measurements: (Ordereddict) Data from the camera

        Returns
        -------
        tuple: the coordinate of the center of the beam
        """

        key = list(measurements['CameraTIS']['data2D'].keys())[0]  # so it can also be used from another plugin having another key
        image = np.array(measurements['CameraTIS']['data2D'][key]['data'])
        
        # self.img1b.setImage(self.roi.getArrayRegion(image, self.img1a), levels=(0, image.max()))
        self.fringes = np.mean(self.roi.getArrayRegion(image, self.img1a),axis = 0)

        S = ft(self.fringes)

        
        x_min = int(self.lr.getRegion()[0])
        x_max = int(self.lr.getRegion()[1])
        # phase_roi = np.unwrap(np.angle(S[x_min:x_max]))   

        phase_roi = np.angle(S)[x_max]
     
        
        if self.settings.child('show_plots','show_camera').value():
            self.img1a.setImage(image)
        if self.settings.child('show_plots','show_roi').value():
            self.plotfringes.setData(y = self.fringes)
        if self.settings.child('show_plots','show_fft').value():
            self.plottf.setData(y = np.abs(S))
        # if self.settings.child('show_plots','show_phase').value():
        #     self.plotph.setData(y = phase_roi)

        # polyfit_T = np.polynomial.polynomial.polyfit(np.arange(len(phase_roi)), phase_roi, 1)
        # self.input[1] = polyfit_T[1]
        self.input[1] = phase_roi
        self.input = np.unwrap(self.input)
        self.input[0] = self.input[1]
        self.phi = self.input[1]

        phi_distance = self.phi*self.lbda/2*np.pi # nm
        self.phi_time = phi_distance/self.c # as

        if self.i_loop%(self.N_rms) == 0:
            self.i_loop = 0
            print(np.std(self.delays))
        else:
            self.i_loop += 1
            
        self.delays[self.i_loop] = self.phi_time

        # self.delays[self.i_loop%(self.N_rms - 1)] = self.phi_time
        # if self.i_loop > self.N_rms - 1:
        #     print(np.std(self.delays))
        # self.i_loop += 1
        return InputFromDetector([self.phi])


        # return InputFromDetector([0])


    def convert_output(self, output, dt, stab=True):
        """
        Convert the output of the PID in units to be fed into the actuator
        Parameters
        ----------
        output: (float) output value from the PID from which the model extract
         a value of the same units as the actuator

        Returns
        -------
        list: the converted output as a list (if there are a few actuators)

        """
        self.curr_phase = output[0]
        # self.distance_order = 10**3*self.curr_phase*self.lbda/(2*np.pi*self.c) #in mum
        self.distance_order = 10**-3*self.curr_phase*self.lbda/(2*2*np.pi)
        return OutputToActuator(mode='rel', values= [self.distance_order])





def main():
    from pymodaq.dashboard import DashBoard
    from pymodaq.daq_utils.daq_utils import get_set_preset_path
    from pymodaq.daq_utils import gui_utils as gutils
    from pathlib import Path
    from PyQt5 import QtWidgets
    from pymodaq.pid.pid_controller import DAQ_PID

    import sys
    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QMainWindow()
    area = gutils.DockArea()
    win.setCentralWidget(area)
    win.resize(1000, 500)
    win.setWindowTitle('PyMoDAQ Dashboard')




    dashboard = DashBoard(area)
    file = Path(get_set_preset_path()).joinpath("BeamSteering.xml")
    if file.exists():
        dashboard.set_preset_mode(file)
        # prog.load_scan_module()
        pid_area = gutils.DockArea()
        pid_window = QtWidgets.QMainWindow()
        pid_window.setCentralWidget(pid_area)

        prog = DAQ_PID(pid_area, dashboard.modules_manager)
        pid_window.show()
        pid_window.setWindowTitle('PidController333')
        QtWidgets.QApplication.processEvents()


    else:
        msgBox = QtWidgets.QMessageBox()
        msgBox.setText(f"The default file specified in the configuration file does not exists!\n"
                       f"{file}\n"
                       f"Impossible to load the DAQ_PID Module")
        msgBox.setStandardButtons(msgBox.Ok)
        ret = msgBox.exec()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()



