from astropy.io import fits
import astropy.stats
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import csv
import yaml

from lsst.ts import salobj
from lsst_efd_client import EfdClient, resample
from astropy.time import Time, TimeDelta





async def Setup():
    d = open('_config.yaml','r')
    config = yaml.safe_load(d)
    csc_ = {}
    domain = salobj.Domain()
    if 'SC_ELEC' in config["GENERAL"]["CSC_INCLUDE"]:
        config_name = 'SC_ELEC'
        sc_electrometer = salobj.Remote(name=config[config_name]['NAME'], domain=domain, index=config[config_name]['NAME']) 
        csc[config_name] = sc_electrometer
        
    elif 'CBP_ELEC' in config["GENERAL"]["CSC_INCLUDE"]:
        config_name = 'CBP_ELEC'
        cbp_electrometer = salobj.Remote(name=config[config_name]['NAME'], domain=domain, index=config[config_name]['NAME']) 
        csc[config_name] = cbp_electrometer

        
    elif 'CBP' in config["GENERAL"]["CSC_INCLUDE"]:
        config_name = 'CBP'
        cbp = salobj.Remote(name=config[config_name]['NAME'], domain=domain, index=config[config_name]['NAME']) 
        csc[config_name] = cbp
        
    elif 'SPECTRO' in config["GENERAL"]["CSC_INCLUDE"]:
        config_name = 'SPECTRO'
        fiber_spectro = salobj.Remote(name=config[config_name]['NAME'], domain=domain, index=config[config_name]['NAME']) 
        csc[config_name] = fiber_spectro
        
    elif 'LASER' in config["GENERAL"]["CSC_INCLUDE"]:
        config_name = 'LASER'
        laser = salobj.Remote(name=config[config_name]['NAME'], domain=domain, index=config[config_name]['NAME']) 
        csc[config_name] = laser
        
    for name, csc in csc_.items():
        await csc.start_task
        await salobj.set_summary_state(csc, salobj.State.ENABLED, override=config[name]['NAME'], timeout=20)
        
    for elec_name in ['SC_ELEC','CBP_ELEC']:
        elec = csc_[elec_name]
        await elec.cmd_setIntegrationTime.set_start(intTime=config[elec_name]['INTTIME'])
        await elec.cmd_setMode.set_start(mode=config[elec_name]['MODE'])
        await elec.cmd_setRange.set_start(setRange=config[elec_name]['RANGE'])
        await electrometer.cmd_performZeroCalib.set_start(timeout=10)
await electrometer.cmd_setDigitalFilter.set_start(activateFilter=False, activateAvgFilter=False, activateMedFilter=False, timeout=10)
        
    # laser configuration
    # fiber spectrograph config?

async def start_elec_measure(electrometers):
    if isinstance(electrometers, list):
        to_run = []
        for elec in electrometers:
            #await elec_id[elec].cmd_performZeroCalib.set_start(timeout=10)
            to_run.append(csc_[elec].cmd_startScan.set_start(timeout=10))
        await asyncio.gather(*to_run)

async def stop_elec_measure(electrometers):
    if isinstance(electrometers, list):
        to_run = []
        for elec in electrometers:
            to_run.append(csc_[elec].cmd_stopScan.set_start(timeout=10))
        await asyncio.gather(*to_run)
         
async def take_elec_measurement(exp_time, electrometers):
    await start_measure(electrometers)
    await asyncio.sleep(exp_time)
    await stop_measure(electrometers)
    
async def spect_meas(exp_time):
    FiberSpectrograph.evt_largeFileObjectAvailable.flush()
    tmp1 = await FiberSpectrograph.cmd_expose.set_start(duration=exp_time, numExposures=1)
    
    return filename
    
async def TakeData():
    start_time
    #Turn on light source
    await asyncio.gather([take_elec_measurement(config['SC_ELEC']['EXPTIME'], ['SC_ELEC','CBP_ELEC']), spect_meas(config['SPECTRO']['EXPTIME'])])
    
    #Turn off the light source
    end_time
    
def transfer_file(filen):
    bashCommand = "scp saluser@toonie.tu.lsst.org:/home/saluser/develop/electrometerFitsFiles/{} {}".format(filen, dir_)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    
async def get_elec_data(start_time, end_time):
    # Electrometer Data
    start_log_topic = 'lsst.sal.Electrometer.command_startScan'
    start_df =  await client.select_time_series(start_log_topic, ['salIndex'], start=start_time, end=end_time)
    start_df['message'] = 'startScan Completed'
    start_df['functionName'] = 'do_startScan'
    msg_log_topic = 'lsst.sal.Electrometer.logevent_logMessage'
    msg_df = await client.select_time_series(msg_log_topic,['salIndex','functionName','message'], start=start_time, end=end_time)
    elec_df = pd.concat([start_df, msg_df])
    elec_df.sort_index(inplace=True)
    elec_df = pd.concat([start_df, msg_df])
    elec_df.sort_index(inplace=True)
    
    data = {}
    for e_id, elec in elec_id.items():
        scans = elec_df[elec_df.salIndex == e_id]
        scans.reset_index(inplace=True)
        if len(scans) > 0:
            for i, row in scans.iterrows():
                if row['functionName'] == 'write_fits_file':
                    try:
                        file_row = scans.iloc[i]
                        filen = file_row['message'].split(' ')[-1]
                        summary_row = scans.iloc[i+1]
                        if 'Scan Summary' in summary_row['message']:
                            x = summary_row['message']
                            x = x.split(':')[1].split(',')
                            mean_ = float(x[0].strip(' ').strip('['))
                            median_ = float(x[1].strip(' '))
                            std_ = float(x[2].strip(' ').strip(']'))
                            data[e_id] = [filen, mean_, median_, std_]
                            break
                    except:
                            print('error')
                            
    for e_id, d in data.items():
        filen = d[0]
        transfer_file(filen)
        start_time, end_time = os.path.splitext(filen)[0].split('_')
        ff = os.path.join(dir_, filen)
        hdu = fits.open(ff)
        print(hdu.info)
        edata = fits.open(ff)[0].data
        values = edata[1]
        time_delta = edata[0]
        meas_time = [Time(t + float(start_time) - 37, scale='tai',format='unix_tai') for t in time_delta]
        times = [x.datetime64 for x in meas_time]
        data[e_id].append(values)
        data[e_id].append(time_delta)
        data[e_id].append(pd.to_datetime(times))
        
async def get_spectro_data(start_time, end_time):
    spectro_dir = os.path.join(config['GENERAL']['DATA_DIR'],'fiberSpec_files')
    spectro_log_topic = 'lsst.sal.FiberSpectrograph.logevent_largeFileObjectAvailable'
    spectro_df =  await client.select_time_series(spectro_log_topic, ['salIndex','url'], start=start_time, end=end_time)
    this_spectro_df = spectro_df[spectro_df['salIndex'] == config['SPECTRO']['ID']]
    spectra = []
    for url in this_spectro_df['url']:
        filename=lfa.url.split('FiberSpectrograph')[-1]
        fs_file = os.path.join(spectro_dir, f)
        os.system('curl {} --output {}'.format(url, fs_file))
        hdu = fits.open(fs_file)
        wave = hdu[1].data['wavelength'][0].flatten()
        spectra.append(hdu[0].data)
    waves = np.array(wave, dtype=[('Wavelength', '<f8')])
    return spectra, waves
        
        
        
async def GetData(start_time, end_time):
    await get_elec_data(start_time, end_time)
    await get_spectro_data(start_time, end_time)

def SaveData(config_file, ConfigData, ):
    #This is in cbp_data_file
    
    