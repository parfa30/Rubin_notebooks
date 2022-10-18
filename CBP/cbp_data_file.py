from astropy.io import fits
import numpy as np
import pandas as pd
import os

class CBPDataFile(object):
    def __init__(self, conf_file, ConfData, SpectroData, CBPData, SCData):
        """
        conf_file: has general configuration information and is used to configure the system/CSCs 
        and general information about the test
        ConfData: data collected from CSCs on configuration of system when test is performed
        SpectroData: data collected from Fiber Spectrograph
        CBPData: data collected from photodiode on integrating sphere feeding CBP with Keithley 6517B
        SCData: data collected from solar cell using Keithley 6517B
        """
        self.conf_file = conf_file
        self.save_dir = '/home/parfa30/DATA/CBP'
        self.hdulist = fits.HDUList()
        self.data_tables = {'SPECTRO':SpectroData, 'CBP':CBPData, 'SOLARCELL':SCData}
        
        
        
    def create_header(self):
        hdr = fits.Header()
HIERARCH spectroexptime = 1.0                                                   
HIERARCH keithleydischarge_time = 2                                             
HIERARCH lasernburst = 5                                                        
HIERARCH keithleyrange = 2E-07                                                  
HIERARCH keithleynsamples = 142                                                 
HIERARCH keithleyrate = 1                                                       
HIERARCH keithleyfunc = 'CHAR    '                                              
HIERARCH keithleytimer = 0.02                                                   
HIERARCH laserwavelength = 859.0                                                
HIERARCH lasernpulses = 31                                                      
HIERARCH laserdelay_before = 0.2                                                
HIERARCH laserdelay_after = 0.231                                               
LASERQSW= 'MAX     '                                                            
HIERARCH laserdelay_start = 0.5                                                 
HIERARCH laserwheelfilter = 'red532  '                                          
HIERARCH cbpmountazalt = '(10, 6) '                                             
HIERARCH keysightmode = 'CHAR    '                                              
HIERARCH keysightrang = 2E-06                                                   
HIERARCH keysightnplc = 0.1                                                     
HIERARCH keysightnsamples = 1427                                                
HIERARCH keysightinterval = 0.002                                               
HIERARCH keysightdelay = 0.0                                                    
HIERARCH focuserposition = 8                                                    
EXPNUM  =               111113                                                  
EXPTYPE = 'NOIMG   '                                                            
HIERARCH cbppinhole = '5mm     '                                                
PD_NAME = 'SM05PD3A'                                                            
SC_NAME = 'SC_18_1.7'                                                           
DISTANCE= 'short   ' 

        d_scheduler['keithley.func'] = 'CHAR'
d_scheduler['keithley.timer'] = 0.02
d_scheduler['keithley.rate'] = 1
d_scheduler['keithley.range'] = 2e-7

d_scheduler['laser.nburst'] = 5
d_scheduler['laser.delay_start'] = 0.5

{'keithley': False,
               'keithley.range': 'IDLE',
               'keithley.nsamples': 'IDLE',
               'keithley.rate': 1,
               'keithley.func': 'CHAR',
               'keithley.timer': 0.1,
               'keithley.discharge_time': 2,
               'spectro': False,
               'spectro.exptime': 1.0,
               'laser': False,
               'laser.wavelength': 520.0,
               'laser.npulses': 40,
               'laser.nburst': 1,
               'laser.delay_before': 2,
               'laser.delay_after': 2,
               'laser.qsw': 'MAX',
               'laser.delay_start': 1,
               'laserwheel': False,
               'laserwheel.filter': 'EMPTY',
               'cbpmount': True,
               'cbpmount.azalt': 'IDLE',
               'digitalanalyser': False,
               'digitalanalyser.duration': 1,
               'keysight': False,
               'keysight.mode': 'CHAR',
               'keysight.rang': '2e-6',
               'keysight.nplc': 0.1,
               'keysight.nsamples': 100,
               'keysight.delay': 0,
               'keysight.interval': 0.002,
               'filterwheel': True,
               'filterwheel.filter': 'PINHOLE',
               'camera': True,
               'camera.exptime': 1,
               'camera.temperature': -70,
               'camera.gain': 2,
               'camera.hspeed': 2,
               'camera.vspeed': 1,
               'camera.shutter_mode': 'AUTO',
               'camera.trigdelay': 0,
               'camera.closing_time': 40,
               'camera.opening_time': 20,
               'ps': True,
               'ps.volet': 'ON',
               'focuser': True,
               'focuser.position': 8,
               'mount': True,
               'mount.az_alt':'IDLE'}

niter = stardice.update_header.ask_user_niter()
pinhole = stardice.update_header.ask_user_pinhole()
SC_name = stardice.update_header.ask_user_solarcell()
PD_name = stardice.update_header.ask_user_photodiode()
distance = stardice.update_header.ask_user_length()
irectory = '/home/dice/stardice/data/cbp_bench2/golden_sample'
basename = os.path.splitext(os.path.basename(__file__))[0]
outdir = stardice.laser_lib.create_directory(directory, basename, user_subdir)
        for key, val in info.items():
        try:
            hdr[key.upper()] = val
        except:
            if key == 'end':
                hdr['STOP'] = val.isot
            elif key == 'start':
                hdr[key.upper()] = val.isot
            else:
                hdr[key.upper()] = str(val)
        if key == 'elec':
            for elec in info[key]:
                hdr['{}_OFFSET'.format(elec)] = x_offset[elec]
                
        
        empty_primary = fits.PrimaryHDU(header=hdr)
        self.hdulist.append(empty_primary)
        
    def create_BinTable(self,data, name):
        t = Table.from_pandas(data)

        self.hdulist.append(fits.BinTableHDU(t.as_array(), name=name))
        
    def data_file(self):
        self.create_header()
        for name, data_table in self.data_tables.items():
            if data_table != None:
                self.create_BinTable(data_table, name)
                
        time = Time.now().isot
        test_name = 'CBP_calibration_{}'.format(time)
        self.hdulist.writeto(f'{os.path.join(self.save_dir,test_name)}.fits', overwrite=True)
                
        #Save file somewhere

        
        
        
    
    
