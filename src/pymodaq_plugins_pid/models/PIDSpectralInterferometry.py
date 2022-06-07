from PyQt5 import QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
from pyqtgraph.dockarea import Dock



from pyqtgraph.Qt import QtCore
import pyqtgraph as pg

from scipy.interpolate import interp1d
from pymodaq.daq_utils.daq_utils import linspace_step
from pymodaq.pid.utils import PIDModelGeneric, OutputToActuator, InputFromDetector
from pymodaq.daq_utils.plotting.viewer1D.viewer1D_main import Viewer1D
from pymodaq.daq_utils.daq_utils import Axis
from pymodaq.daq_utils.daq_utils import DataFromPlugins
from pymodaq.daq_utils.math_utils import ft,ift



class PIDModelSpetralInterferometry(PIDModelGeneric):
    limits = dict(max=dict(state=False, value=1),
                  min=dict(state=False, value=-1),)
    konstants = dict(kp=1.0, ki=0.0, kd=0.0)

    Nsetpoint = 1
    setpoint_ini = [2. for ind in range(Nsetpoint)]
    actuators = ['SmarActMCS']
    actuators_name = ['SmarAct']
    # detector = 'Spectro'
    detectors_name = ['Spectro']
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
        {'title': 'Wavelength (nm)','name' : 'lambda_0','type':'float','value':800}
    ]


    def __init__(self, pid_controller):
        super().__init__(pid_controller)

        ##Plot du spectre
        self.dock_calib = Dock('Spectrum flame')   #fenetre
        widget_calib = QtWidgets.QWidget()     #widget
        print('Init Widget: There should be a dock')
        # NEW DOCK
        self.pid_controller.dock_area.addDock(self.dock_calib)     #ajout dans pid controller
        # ## NEW WIDGET
        w = pg.PlotWidget(title="Dock 6 plot")       #affichage
        # ## NEW PLOT
        self.plotData = w.plot()
        # ## ADD WIDGET TO DOCK
        self.dock_calib.addWidget(w)
        # ## SET DATA TO PLOT
        self.plotData.setData(y = np.random.normal(size=100))
        w.setLabel('bottom','Frequency',units = 'PHz')


        ##Plot de la TFI
        self.dock_calib = Dock('Inverse Fourier transform')
        widget_calib = QtWidgets.QWidget()
        print('Init Widget: There should be a dock')
        self.pid_controller.dock_area.addDock(self.dock_calib)
        u = pg.PlotWidget(title="Inverse Fourier transform ")
        self.plotitf = u.plot()
        self.dock_calib.addWidget(u)
        self.plotitf.setData(y = np.random.normal(size=100))
        u.setLabel('bottom','Time',units = 'fs ')

        ##Plot de la TF
        self.dock_calib = Dock('Fourier transform')
        widget_calib = QtWidgets.QWidget()
        self.pid_controller.dock_area.addDock(self.dock_calib)
        r = pg.PlotWidget(title="Fourier transform ")
        self.plottf = r.plot()
        self.dock_calib.addWidget(r)
        self.plottf.setData(y = np.random.normal(size=100))
        r.setLabel('bottom','Time',units = 'PHz ')


        #ROI sur la tf
        self.lr = pg.LinearRegionItem([10,20])
        self.lr.setZValue(-10)
        r.addItem(self.lr)



        ##Plot de la phase
        self.dock_calib = Dock('Phase')
        widget_calib = QtWidgets.QWidget()
        print('Init Widget: There should be a dock')
        self.pid_controller.dock_area.addDock(self.dock_calib)
        ph = pg.PlotWidget(title="Phase")
        self.plotph = ph.plot()
        self.dock_calib.addWidget(ph)
        test_signal = np.random.normal(size = 100)
        self.plotph.setData(y =np.random.normal(size = 100))
        
        self.isNotCalibrated = True
        self.wavelength = []
        self.frequency = []






    def updateGUI(self,param):
        self.viewer_calib.show_data(param)


    def update_settings(self, param):
        """
        Get a parameter instance whose value has been modified by a user on the UI
        Parameters
        ----------
        param: (Parameter) instance of Parameter object
        """
        if param.name() == '':
            pass

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

        self.lbda = self.settings.child('lambda_0').value()
        key = list(measurements['Spectro']['data1D'].keys())[0]  # so it can also be used from another plugin having another key
        wavelength_spectrum = np.array(measurements['Spectro']['data1D'][key]['data'])
        print(self.lbda)
        
        if self.isNotCalibrated:       ### MAKE FIT TO CONVERT SPECTRUM FROM WAVELENGTH TO FREQUENCY

            self.wavelength = np.array(measurements['Spectro']['data1D'][key]['x_axis']['data'])
            self.omega_init= self.wavelengthToFrequency(self.wavelength*10**6)
            self.factor_sampling = 1
            self.d_omega = np.max(np.abs(np.diff(self.omega_init))) #choice of the largest frequence resolution
            self.n_sampling = int(abs(self.omega_init[-1] - self.omega_init[0])/self.d_omega)*self.factor_sampling
            print('fit done')

            self.omega = np.linspace(self.omega_init[0], self.omega_init[-1], num=self.n_sampling, endpoint=True)
            print('linear omega done')
            self.isNotCalibrated = False
        

        
        self.fit = interp1d(self.omega_init, wavelength_spectrum, kind='linear')
        self.omega_spectrum = self.fit(self.omega)
        self.plotData.setData( x =self.omega,y = self.omega_spectrum )
 

        ##Prints used to find how to call spectro data
        #options = list(measurements['Spectro']['data1D'][key].keys())
        # [print(measurements['Spectro']['data1D'][key][option]) for option in options]
        # print('key name:')
        # [print(option) for option in options]
        # print('data:')
        # [print(measurements['Spectro']['data1D'][key][option]) for option in options]

        


        # N_samp_power = self.settings.child('inverse_fourier','N_samp_power').value()
        # print(N_samp_power)
        S = ift(wavelength_spectrum)
        #t,S = self.do_inverse_Fourier(self.omega,self.omega_spectrum,2**N_samp_power)
        E = ft(S)


        t_min = int(self.lr.getRegion()[0])
        t_max = int(self.lr.getRegion()[1])

        phase_roi = np.unwrap(np.angle(E[t_min:t_max]))
        polyfit_T = np.polynomial.polynomial.polyfit(np.arange(len(phase_roi)), phase_roi, 1)


        self.plotitf.setData(y = np.abs(S))
        self.plottf.setData(y = np.abs(E))
        self.plotph.setData(y = phase_roi)

        return InputFromDetector([polyfit_T[1]])



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
        #print('output converted')

        self.curr_output = output
        # print(self.modules_manager.actuators_name)
        # print(self.modules_manager.detectors_name)
        # self.paused = False
        # self.modules_manager.move_actuators(output,'abs',polling=False)
        

        # central_frequency = omega_c
        ## CONVERSION FROM RADIAN TO TIME
        omega_0 = 2*np.pi*3/(self.lbda*10**-9)   
        phase_time = output[0]/self.lbda

        ## CONVERION FROM TIME TO DISTANCE
        phase_distance = phase_time*3*10**11/(2*omega_0) # in mm
        # output = 1.0
        # print(phase_distance)
        return OutputToActuator(mode='rel', values=[phase_distance])
        # return OutputToActuator(mode='abs', values=output)



    def ApplyFilter(x,y,start=False,end=False):
        if not(isinstance(start, bool)):
            mask = x >= start
        else:
            mask = np.ones_like(x).astype(bool)
        if not(isinstance(end, bool)):
            mask = np.logical_and(mask,x <= end)
        x = x[mask]
        y = y[mask]
        return x, y
    def wavelengthToFrequency(self,wavelength):
        c = 0.299792458  # in mum/fs
        omega = 2 * np.pi * c / wavelength # 2*pi*c/lambda in PHz
        return omega
    def Jacobian(self,spectra,omega):
        spectra = spectra / omega ** 2 #Jacobian
        order = np.argsort(omega_init) # Reorder frequency axis
        spectra = spectra[:, order] # Apply ordering
        omega_init = omega_init[order]
        return omega,spectra 
    def doInterpolate(self,spectra,omega_init,factor_sampling = 2):
        n_sampling = np.round(omega_init.shape[0], decimals=-3) * factor_sampling
        fit = interp1d(omega_init, spectra, kind='linear')
        omega = np.linspace(omega_init[0], omega_init[-1], num=n_sampling, endpoint=True)
        spectra = fit(omega)
        return omega,spectra
    # Going fo frequency
    def switchToFrequency(self,spectra,wavelength):
        omega_init = self.wavelengthToFrequency(wavelength) * 1e3
        omega_init,spectra = self.Jacobian(spectra,omega_init)
        spectra = spectra - spectra[0]# Remove first value to reduce background
        omega,spectra = self.doInterpolate(self,spectra,omega_init,factor_sampling = 2)
        return omega,spectra




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



