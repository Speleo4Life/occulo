from pywinauto.findwindows    import find_window
from pywinauto.win32functions import SetForegroundWindow
import win32gui
import socket
import numpy as np
import msvcrt
import csv
import os
import pylink
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
from ConsoleTimer import consoleTimer
from psychopy import visual, core, event, sound, prefs
from pylsl import StreamInfo, StreamOutlet
from datetime import datetime
from eyelink_commands import *
from ss_trial_h import *
from load_config_csv import load_cond1_config_file, load_cond2_config_file

#prefs.general['audioLib'] = ['sounddevice']

### List of Cyton (EOG) channels to record
### You can use any subset of channels [1,2,3,4,5,6,7,8]
### N2P, N4P, UP (EOG), DOWN (REFERENCE)
eog_channels = [2,4]

### Color codes for bg / drawing stims
GRY = (0,0,0)
BLK = (-1,-1,-1)
RED = (1,-1,-1)
GRN = (-1,1,-1)

### Get stim locations based on angle
SCREEN_RES = (1920, 1080) #(1920, 1080)
SCREEN_WIDTH_CM = 120.65
SCREEN_DIAG_CM = 138.68
EYE_TO_SCREEN_CM = 50.00 # 96.5 (old setup)
SM = 5.25
SMMD = 7.75
MD = 10.50
LG = 21.00
TOFFSETX_SM = (EYE_TO_SCREEN_CM * np.tan(np.deg2rad(SM)))
TOFFSETX_MD = (EYE_TO_SCREEN_CM * np.tan(np.deg2rad(MD)))
TOFFSETX_LG = (EYE_TO_SCREEN_CM * np.tan(np.deg2rad(LG)))
TOFFSETX_CAL_HORZ_1 = -1*(EYE_TO_SCREEN_CM * np.tan(np.deg2rad(LG)))
TOFFSETX_CAL_HORZ_2 = -1*(EYE_TO_SCREEN_CM * np.tan(np.deg2rad(MD)))
TOFFSETX_CAL_HORZ_3 = -1*(EYE_TO_SCREEN_CM * np.tan(np.deg2rad(SM)))
TOFFSETX_CAL_HORZ_4 = (EYE_TO_SCREEN_CM * np.tan(np.deg2rad(SM)))
TOFFSETX_CAL_HORZ_5 = (EYE_TO_SCREEN_CM * np.tan(np.deg2rad(MD)))
TOFFSETX_CAL_HORZ_6 = (EYE_TO_SCREEN_CM * np.tan(np.deg2rad(LG)))
FG_COLOUR = GRY
BG_COLOUR = BLK


PPCM = np.sqrt(SCREEN_RES[0]**2. + SCREEN_RES[1]**2.) / SCREEN_DIAG_CM

EL_SAMPLE_RATE = 250
isi = 1.0 / EL_SAMPLE_RATE

# Return to fixation minimum gaze duration (ms)
FIX_RETURN_GAZE_DUR = 250
FIX_RETURN_GAZE_FRAMES = int(FIX_RETURN_GAZE_DUR / (1000 * isi))

UDP_IP = "127.0.0.1"
UDP_PORT_TO_EOG  = 50999
UDP_PORT_FROM_EOG  = 51000
UDP_PORT_FROM_LR = 51001

lab_recorder_exe = 'C:\\labstreaminglayer\\Apps\\LabRecorder\\build\\Release\\LabRecorder.exe'

# This will ensure that the Command Prompt window will remain selected after PyGame loads
# Makes the program less annoying to use.
def window_enum_handler(hwnd, resultList):
    if win32gui.IsWindowVisible(hwnd) and win32gui.GetWindowText(hwnd) != '':
        resultList.append((hwnd, win32gui.GetWindowText(hwnd)))

def get_app_list(handles=[]):
    mlst=[]
    win32gui.EnumWindows(window_enum_handler, handles)
    for handle in handles:
        mlst.append(handle)
    return mlst

appwindows = get_app_list()



try: 
    handle = win32gui.FindWindow(0, "C:\\WINDOWS\\system32\\CMD.exe - python2  ss_run_h.py")
    win32gui.SetForegroundWindow(handle)
except:
    try: 
        handle = win32gui.FindWindow(0, 'C:\\WINDOWS\\system32\\cmd.exe - python2  C:\\Users\\Visionlab\\Desktop\\tochi_working_version\\ss_run_h.py')
        win32gui.SetForegroundWindow(handle)
    except:
        pass
    
# Get participant info
print('')

# got_pname = False
got_pname = True
pname = str(999)
while not got_pname:
    pname = raw_input("Enter participant number: ")
    if len(pname) > 4:
        print("Input Error! Participant number may be 4 digits max.")
    else:
        got_pname = True

# ipd = raw_input("Enter participant's interpupilary distance (mm): ")
ipd = 0
fix_shift = 0 # float(ipd) / 20. # Divide by two and convert to centimeters
# cond = 0
cond = 1
# Get condition info, display reminders for proper screen/Eyelink setup
print("""
Please select the condition that you wish to run:
1. Eyes open, Markers on & off (alternating), Luminance
2. Eyes open & closed (alternating), Markers off, No Luminance
""")
while cond == 0:
    try:
        user_cond = int(raw_input("Enter condition number: "))
        if (user_cond > 0 and user_cond < 3)  : #and (user_cond != 2) :
            cond = user_cond
        else:
            raise Exception
    except :
        print("Input Error! Enter a number corresponding to option 1 or 2")

# Display condition-specific instructions for experimenter
# if cond == 1 :
#     markers_on = True
#     user_prac = raw_input("Run calibration trials? (y/n) ")
#     if user_prac.startswith('y') or user_prac.startswith('Y') :
#         run_calib_trials = True
#     else :
#         print("Skip calibration trials for Condition 1? Press ENTER to continue, ESC to quit.") 
#         while True :
#             if msvcrt.kbhit() :
#                 if msvcrt.getch() == b'\x1b' :
#                     quit() 
#                 else :
#                     run_calib_trials = False
#                     break
# else :
#     markers_on = False
#     run_calib_trials = False
run_calib_trials = True
edfn = "P" + pname + "_C" + str(cond) + ".edf"
edf_out_path = "C:\\Users\\Visionlab\\Desktop\\ClosedEyes\\Data\\EDF\\"


# Randomly choose true or false (to be updated with CSV specs)
markersOnFirstBlock = np.random.rand() > 0.5

### Create trial order depending on condition
if cond == 1 :

    print("\nPlease select the config file for Condition 1...")

    trialTargetOrder, trialBlockNumber, blockNumbers, blockMarkersOn, blockTargVert = load_cond1_config_file()

    if not all(np.size(x) > 0 for x in [trialTargetOrder, trialBlockNumber, blockNumbers, blockMarkersOn, blockTargVert]) :
        print('')
        print("Config Error! Could not load or parse CSV file successfully. Exiting...")
        exit()

    n_trials = np.size(trialTargetOrder, axis=0)
    num_blocks = np.size(np.unique(trialBlockNumber), axis=0)
    block_length = n_trials/num_blocks

    numCalibTargets = 6
    numTrialsPerTarget_HorzCalib = 1
    caliOrderHorz = [3,2,1,4,5,6]
    # caliOrderHorz = [3,2,1,4,5,6, # T,R,S,Q,W,X,Y,Z
    #                  3,2,1,4,5,6,
    #                  3,2,1,4,5,6]
else : 
    print("\nPlease select the config file for Condition 2...")

    trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen, blockTargVert = load_cond2_config_file()

    if not all(np.size(x) > 0 for x in [trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen, blockTargVert]) :
        print('')
        print("Config Error! Could not load or parse CSV file successfully. Exiting...")
        exit()

    n_trials = np.size(trialTargetOrder, axis=0)
    num_blocks = np.size(np.unique(trialBlockNumber), axis=0)
    block_length = n_trials/num_blocks

    numCalibTargets = 3
    numTrialsPerTarget_HorzCalib = 1
    caliOrderHorz = [3,2,1,4,5,6] # 1,2,3,4,5,6
                    # 3,2,1,4,5,6,
                    # 3,2,1,4,5,6]
# elif cond == 3:
#     numTargets = 4
#     numTrialsPerLetter = 12
#     trialLetterOrder = np.repeat(np.arange(numTargets), numTrialsPerLetter)
#     np.random.shuffle(trialLetterOrder)

#     numMarkerTrials = 5
#     d = np.floor_divide(numMarkerTrials,numTargets)
#     r = np.remainder(numMarkerTrials,numTargets)
#     practiceMarkerTrials = np.tile(np.arange(numTargets), d)
#     r_ints = np.random.randint(0, high=numTargets, size=r)
#     practiceMarkerTrials = np.concatenate((practiceMarkerTrials, r_ints))
#     np.random.shuffle(practiceMarkerTrials)

#     trialLetterOrder = np.concatenate((practiceMarkerTrials, trialLetterOrder))
#     n_trials = np.size(trialLetterOrder,0)

#     num_blocks = 1
#     block_length = n_trials
# else: # cond == 4
#     numMarkerTrials = 0
#     numTargets = 4
#     numTrialsPerLetter = 12
#     trialLetterOrder = np.repeat(np.arange(numTargets), numTrialsPerLetter)
#     np.random.shuffle(trialLetterOrder)
#     n_trials = np.size(trialLetterOrder,0)

#     num_blocks = 1
#     block_length = n_trials


# Output demographics information to CSV
exp_date = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
targets_str = ' '.join(str(l) for l in trialTargetOrder)
if cond == 1:
    block_markers_str = ' '.join(str(b) for b in blockMarkersOn)
elif cond == 2:
    block_markers_str = ' '.join(str(b) for b in blockEyesOpen)
else:
    block_markers_str = ''

dems_out = "C:\\Users\\Visionlab\\Desktop\\ClosedEyes\\Data\\ClosedEyes_Demographics.csv"
dems_exists = os.path.isfile(dems_out)
dems_header = ["Participant","Condition","IPD","Date","Block Size","Block Markers","N Trials"]
dems_row = [pname, str(cond), str(ipd), exp_date, str(block_length), block_markers_str, str(n_trials)]

if not dems_exists :
    with open(dems_out, 'wb') as csvfile :
        demswriter = csv.writer(csvfile, delimiter=',')
        demswriter.writerow(dems_header)

with open(dems_out, 'ab') as csvfile :
    demswriter = csv.writer(csvfile, delimiter=',')
    demswriter.writerow(dems_row)


# eog_to_UDP_str = '['+str(UDP_IP)+','+str(UDP_PORT_TO_EOG)+']'
# eog_from_UDP_str = '['+str(UDP_IP)+','+str(UDP_PORT_FROM_EOG)+']'
eog_UDP_str = '['+str(UDP_IP)+','+str(UDP_PORT_FROM_EOG)+']'
lr_UDP_str = '['+str(UDP_IP)+','+str(UDP_PORT_FROM_LR)+']'
print('')
print("Opening UDP ports for EOG: send=" + eog_UDP_str + ", receive=" + lr_UDP_str)
# # Open UDP socket to send record on/off messages to OpenBCI_LSL
# socket_eog_to = socket.socket(socket.AF_INET, # Internet
#                               socket.SOCK_DGRAM) # UDP
# socket_eog_to.sendto("e", (UDP_IP, UDP_PORT_TO_EOG)) # turn off LSL streaming until first trial

# Open UDP socket to send record on/off messages to OpenBCI_LSL
socket_eog = socket.socket(socket.AF_INET, # Internet
                           socket.SOCK_STREAM) # TCP
socket_eog.bind((UDP_IP, UDP_PORT_FROM_EOG))
socket_eog.setblocking(0) # set to non-blocking port (for the while loop)
socket_eog.listen(1)


# # Open UDP socket to receive stream start message from OpenBCI_LSL
# socket_eog_from = socket.socket(socket.AF_INET, # Internet
#                                 socket.SOCK_DGRAM) # UDP
# socket_eog_from.bind((UDP_IP, UDP_PORT_FROM_EOG))
# socket_eog_from.setblocking(0) # set to non-blocking port (for the while loop)

# Open UDP socket to receive exp start/stop messages from LabRecorder
print('')
print('')
print("Opening UDP ports for Lab Recorder: receive=" + lr_UDP_str)
print('')
socket_lr = socket.socket(socket.AF_INET, # Internet
                          socket.SOCK_DGRAM) # UDP
socket_lr.bind((UDP_IP, UDP_PORT_FROM_LR))
socket_lr.setblocking(0) # set to non-blocking port (for the while loop)


### Initialize OpenBCI_LSL
print('')
print('Initializing OpenBCI, perform EOG calibration in popup window.')
print("Please wait for BCI Message...")
print('')
# try:
#     openbci_process = subprocess.Popen(['python', 'C:/Users/Visionlab/Desktop/ClosedEyes/OpenBCI_LSL/openbci_lsl.py', '--channels', str(eog_channels)])
#     time.sleep(20)
#     eog_conn, eog_addr = socket_eog.accept()
# except Exception as e:
#     raise

# print('')
# print('BCIMSG: You may now start streaming in OpenBCI.')
# print('')

# while 1:
#     opbci_msg = b"D"
#     try:
#         # opbci_msg = socket_eog_from.recv(2) # buffer size is 1 bytes (1 character + null)   
#         opbci_msg = eog_conn.recv(2) # buffer size is 1 bytes (1 character + null)   
#     except:
#         pass

#     if opbci_msg.startswith(b"s"):
#         eog_conn.send(b"e") # turn off LSL streaming until first trial
#         print('BCIMSG: Press ENTER to continue.')
#         break

eog_conn = None
try:
    myeog_process = subprocess.Popen(['python', 'EOG_circuit_testing/StreamSerialEOG.py'])
    time.sleep(2)
except Exception as e:
    raise

# This will ensure that the Command Prompt window will remain selected after the message is displayed
# Makes the program less annoying to use.
try:
    SetForegroundWindow(find_window(best_match='ss_run_h.py'))
except:
    pass

while 1:
    if msvcrt.kbhit():
        if ord(msvcrt.getch()) == 13: # Enter
            break

### Initialize Eyelink
if cond != 4 :
    print('')
    print("Initializing Eyelink...")
    el_outlet, elTk = InitEyelink(sr=EL_SAMPLE_RATE, fname=edfn)

if cond == 4 :
    print('')
    print("COND 4: skipping Eyelink initialization...")
    print("Select only 'EventMarker' and 'OpenBCI' streams in LabRecorder!")

### Create LSL outlet for event markers
# first create a new stream info (here we set the name to MyMarkerStream,
# the content-type to Markers, 1 channel, irregular sampling rate,
# and string-valued data) The last value would be the locally unique
# identifier for the stream as far as available, e.g.
# program-scriptname-subjectnumber (you could also omit it but interrupted
# connections wouldn't auto-recover). The important part is that the
# content-type is set to 'Markers', because then other programs will know how
#  to interpret the content
try:
    ev_info = StreamInfo('EventMarkers', 'Markers', 1, 0, 'string', 'myuidw43536')
    ev_outlet = StreamOutlet(ev_info)
    ### Pause for experimenter to connect LSL stream
    print("Established LSL event stream outlet.")
    # print("REMINDER: Set Participant (%p) to 'P" + str(pname) + "' and Session (%s) to 'S00" + str(cond) + "'")
    # print('')
    # raw_input("Start recording in LabRecorder now. Press ENTER to continue...")
except:
    print("Could not create LSL outlet.")

print('')
print("Initializing Lab Recorder...")
print("Start stream on LabRecorder to begin experiment (or press 's' to skip & begin without Lab Recorder).")
print('')
lr_exp_string = '[' + pname + "," + str(cond) + ']'
if cond != 4:
    labrecorder_process = subprocess.Popen([lab_recorder_exe, '-c', 'CE_LabRecorderSetting.cfg', '-u', lr_UDP_str, '-e', lr_exp_string])
else:
    labrecorder_process = subprocess.Popen([lab_recorder_exe, '-c', 'CE_LabRecorderSetting_noEL.cfg', '-u', lr_UDP_str, '-e', lr_exp_string])
   
while 1:
    lr_msg = b"D"
    try:
        lr_msg = socket_lr.recv(2) # buffer size is 1 bytes (1 character + null)   
    except:
        pass

    if lr_msg.startswith(b"b"):
        break

    if msvcrt.kbhit():
        if ord(msvcrt.getch()) == 115: # 's'
            break

    if labrecorder_process.poll() is not None :
        print('Lab Recorder has quit, press ENTER to continue without LSL stream recording, or ESC to exit.')
        while 1:
            if msvcrt.kbhit():
                if ord(msvcrt.getch()) == 27: # ESC
                    if openbci_process.poll() is None:
                        openbci_process.kill()
                    quit()
                if ord(msvcrt.getch()) == 13: # Enter
                    break
        break

# This will ensure that the Command Prompt window will remain selected after LabRecorder loads
# Makes the program less annoying to use.
try:
    SetForegroundWindow(find_window(best_match='ss_run2.py'))
except:
    pass

if cond == 4 :
    print('')
    raw_input("COND " +str(cond) + ": Turn off TV now! Press ENTER to continue...")

#create a window
w = visual.Window(monitor="Samsung UN55KU6290", screen=1, size=SCREEN_RES, fullscr=True, color=BG_COLOUR, units="cm")
w.mouseVisible = False

### Create text messages
st_cal_msg = createMessage(w, FG_COLOUR, 'Press any key to begin calibration.')
end_cal_msg = createMessage(w, FG_COLOUR, 'End of Calibration.\nPress any key to continue to the experiment.')

st_exp_msg = createMessage(w, FG_COLOUR, 'Please wait for experimenter instructions.\nPress any key to continue.')
end_cond1_msg = createMessage(w, FG_COLOUR, 'End of Condition.\nPress any key to continue.')
end_exp_msg = createMessage(w, FG_COLOUR, 'End of Experiment.\nPress any key to exit.')

exit_msg = createMessage(w, FG_COLOUR, 'Exiting...')



# set up a custom graphics envrionment (EyeLinkCoreGraphicsPsychopy) for calibration
genv = EyeLinkCoreGraphicsPsychoPy(elTk, w)

# play feedback beeps during calibration/validation, disabled by default
genv.enableBeep = True
pylink.setCalibrationSounds("", "", "")
# calibration target size 
genv.targetSize = 14
# Configure the calibration target, could be a 'circle', 
# a movie clip ('movie'), a 'picture', or a 'spiral', the default is a circle
genv.calTarget = 'circle'
pylink.openGraphicsEx(genv)


### Start calibration and trial procedures ##################
pylink.beginRealTimeMode(100)

# Begin Eyelink recording (continuous mode, for trial-by-trial, see trial loop functions)
# if cond != 4 :
#     elTk.startRecording(1, 1, 1, 1)
    
exp_abort = False
if run_calib_trials :
    # Run calibration procedure for Eyelink & EOG linearization
    is_calib = True
    ev_outlet.push_sample(["cal_horz_start"])
    if cond != 4 :
        elTk.sendMessage(["cal_horz_start"])
    exp_abort = ss_runcalibs(w, "H", FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_cal_msg, end_cal_msg, caliOrderHorz, numCalibTargets, numTrialsPerTarget_HorzCalib, isi, TOFFSETX_CAL_HORZ_1,TOFFSETX_CAL_HORZ_2,TOFFSETX_CAL_HORZ_3,TOFFSETX_CAL_HORZ_4,TOFFSETX_CAL_HORZ_5,TOFFSETX_CAL_HORZ_6, fix_shift, TOFFSETX_SM, TOFFSETX_MD, ev_outlet, eog_conn, el_outlet, elTk)
    ev_outlet.push_sample(["cal_horz_end"])
    if cond != 4 :
        elTk.sendMessage(["cal_horz_end"])       
    if not exp_abort :
       print("On vertical test!") 
       ev_outlet.push_sample(["cal_vert_start"])
       if cond != 4 :
           elTk.sendMessage(["cal_vert_start"])
       exp_abort = ss_runcalibs(w, "V", FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_cal_msg, end_cal_msg, caliOrderHorz, numCalibTargets, numTrialsPerTarget_HorzCalib, isi, TOFFSETX_CAL_HORZ_1,TOFFSETX_CAL_HORZ_2,TOFFSETX_CAL_HORZ_3,TOFFSETX_CAL_HORZ_4,TOFFSETX_CAL_HORZ_5,TOFFSETX_CAL_HORZ_6, fix_shift,TOFFSETX_SM, TOFFSETX_MD, ev_outlet, eog_conn, el_outlet, elTk)
       ev_outlet.push_sample(["cal_vert_end"])
       if cond != 4 :
           elTk.sendMessage(["cal_vert_end"])

if not exp_abort :
    # Run trials
    ev_outlet.push_sample(["exp_start"])
    if cond != 4 :
        elTk.sendMessage(["exp_start"])

    ev_outlet.push_sample(["cond" + str(cond) + "_start"])
    if cond != 4 :
        elTk.sendMessage(["cond" + str(cond) + "_start"])
        if cond == 1:
            ss_runtrials(w, FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_exp_msg, end_cond1_msg, trialTargetOrder, trialBlockNumber, blockNumbers, blockMarkersOn, blockTargVert, isi, TOFFSETX_SM,TOFFSETX_MD, TOFFSETX_LG, fix_shift, ev_outlet, eog_conn, el_outlet, elTk)           
        elif cond == 2:
            ss_runtrials_cond2(w, FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_exp_msg, end_cond1_msg, trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen, isi, TOFFSETX_SM,TOFFSETX_MD, TOFFSETX_LG, fix_shift, ev_outlet, eog_conn, el_outlet, elTk)
    else :
        if cond == 1:
            ss_runtrials(w, FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_exp_msg, end_cond1_msg, trialTargetOrder, trialBlockNumber, blockNumbers, blockMarkersOn, blockTargVert, isi, TOFFSETX_SM,TOFFSETX_MD, TOFFSETX_LG, fix_shift, ev_outlet, eog_conn, el_outlet, elTk)           
        elif cond == 2:
             ss_runtrials_cond2(w, FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_exp_msg, end_cond1_msg, trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen, isi, TOFFSETX_SM,TOFFSETX_MD, TOFFSETX_LG, fix_shift, ev_outlet, eog_conn, el_outlet, elTk)
    
    ev_outlet.push_sample(["cond" + str(cond) + "_end"])
    if cond != 4 :
        elTk.sendMessage(["cond" + str(cond) + "_end"])


# Ask if we would like to run Cond 2 (Non-Illuminated)
run_cond2 = False
if not exp_abort :
    if cond == 1:
                
        w.close()
        w.mouseVisible = True

        # user_prac = raw_input("End of Condition 1. Run Condition 2 now? (y/n) ")
        # if user_prac.startswith('y') or user_prac.startswith('Y') :
        #     run_cond2 = True
        #     cond = 2
        # else :
        #     user_prac = raw_input("Really skip Condition 2? (y/n)") 
        #     if user_prac.startswith('n') or user_prac.startswith('N') :
        #         run_cond2 = True

        # if run_cond2 :
        #     exit_cond2 = False
        #     print("\nPlease select the config file for Condition 2...")

        #     trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen = load_cond2_config_file()

        #     if not all(np.size(x) > 0 for x in [trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen]) :
        #         print('')
        #         print("Config Error! Could not load or parse CSV file successfully. Exiting...")
        #         exit_cond2 = True

        #     n_trials = np.size(trialTargetOrder, axis=0)
        #     num_blocks = np.size(np.unique(trialBlockNumber), axis=0)
        #     block_length = n_trials/num_blocks

        #     print("""
        #     Press ENTER to begin 15 minute period of dark adaptation, 
        #     Press S to skip and begin Cond 2 immediately,
        #     Press ESC to abort Cond 2 and exit.
        #     """)
        #     wait_cond2 = False
        #     while True :
        #         if msvcrt.kbhit() :
        #             ch = msvcrt.getch() 
        #             if ch == b'\x1b' :
        #                 run_cond2 = False
        #                 exit_cond2 = True
        #             elif ch == 's' : 
        #                 break
        #             elif ch == '\r' :
        #                 wait_cond2 = True
        #                 break
            
        #     if not exit_cond2 :
        #         if wait_cond2 :
        #             consoleTimer(15,0) # 15 minute timer coundown in console

        #         print('')
        #         raw_input("COND " +str(cond) + ": Turn off TV now! Press ENTER to continue...")

        #         block_markers_str = ' '.join(str(b) for b in blockEyesOpen)
        #         dems_row = [pname, str(cond), str(ipd), exp_date, str(block_length), block_markers_str, str(n_trials)]
        #         with open(dems_out, 'ab') as csvfile :
        #             demswriter = csv.writer(csvfile, delimiter=',')
        #             demswriter.writerow(dems_row)

        #         #create a window
        #         w = visual.Window(monitor="Samsung UN55KU6290", screen=1, size=SCREEN_RES, fullscr=True, color=BG_COLOUR, units="cm")
        #         w.mouseVisible = False

        #         # Run trials
        #         ev_outlet.push_sample(["cond" + str(cond) + "_start"])
        #         if cond != 4 :
        #             elTk.sendMessage(["cond" + str(cond) + "_start"])
        #             ss_runtrials_cond2(w, FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_exp_msg, end_exp_msg, trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen, isi, TOFFSETX_HORZ_SM, TOFFSETX_HORZ_LG, fix_shift, ev_outlet, eog_conn, el_outlet, elTk)           
        #         else :
        #             ss_runtrials_cond2(w, FG_COLOUR, BG_COLOUR, RED, GRN, PPCM, st_exp_msg, end_exp_msg, trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen, isi, TOFFSETX_HORZ_SM, TOFFSETX_HORZ_LG, fix_shift, ev_outlet, eog_conn, None, None)                
        #         ev_outlet.push_sample(["cond" + str(cond) + "_end"])
        #         if cond != 4 :
        #             elTk.sendMessage(["cond" + str(cond) + "_end"])


ev_outlet.push_sample(["exp_end"])
if cond != 4 :
    elTk.sendMessage(["exp_end"])

# End Eyelink recording after trials are complete (continuous mode only, for trial-by-trial see trial functions)
# if cond != 4 :
#     elTk.stopRecording()

# end_message.setAutoDraw(False)
# exit_message.setAutoDraw(True)
if exp_abort or (run_cond2 is not False):
    try:
        exit_msg.draw()
        w.flip()
    except :
        pass


### Terminate Eyelink stream
pylink.endRealTimeMode()

if cond != 4 :
    elTk.setOfflineMode()     

    # download EDF file to Display PC and put it in local folder ('edfData')
    if exp_abort or (run_cond2 is not False):
        try:
            msg = 'EDF data is transfering from EyeLink Host PC...'
            edfTransfer = visual.TextStim(w, text=msg, color='white')
            edfTransfer.draw()
            w.flip()
            pylink.pumpDelay(500)
        except Exception as e:
            raise
         
    # Close the file and transfer it to Display PC
    elTk.closeDataFile()
    elTk.receiveDataFile(edfn, edf_out_path+edfn)
    elTk.close()

    if exp_abort or (run_cond2 is not False):
        try:
            w.flip()
        except:
            pass

# Save lost of frame intervals
# w.saveFrameIntervals(fileName=None, clear=True)

# Terminate OpenBCI_LSL if needed
# if openbci_process.poll() is None :
#     openbci_process.kill()
# eog_conn.close()


if exp_abort or (run_cond2 is not False):
    try:
        w.close()
        w.mouseVisible = True
    except:
        pass

if labrecorder_process.poll() is None :
    time.sleep(2)
    print("")
    print('Waiting for Lab Recorder, stop recording stream to end experiment.')
    print("")
    
    while 1:
        lr_msg = b"D"
        try:
            lr_msg = socket_lr.recv(2) # buffer size is 1 bytes (1 character + null)
        except:
            pass
        
        if lr_msg.startswith(b"e"):
            break

        if labrecorder_process.poll() is not None :
            break

# Terminate Lab Recorder if needed
if labrecorder_process.poll() is None :
    labrecorder_process.kill()

print('Exiting program...')
print('')

#cleanup
core.quit()