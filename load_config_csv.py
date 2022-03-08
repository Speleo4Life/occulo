import Tkinter, tkFileDialog
import csv
import os
import numpy as np


def load_cond1_config_file() :
    
    # Init output vars (to pass to ss_trial)
    trialTargetOrder = np.array([], dtype=int)
    blockMarkersByTrial = np.array([], dtype=int)
    blockMarkersOn = np.array([], dtype=int)
    trialBlockNumber = np.array([], dtype=int)
    trialTargVert = np.array([], dtype=int)
    blockTargVert = np.array([], dtype=int)
    calOrderID = np.array([], dtype=int)

    # Open system dialog to choose config file
    root = Tkinter.Tk()
    root.withdraw()
    file_path = tkFileDialog.askopenfilename(initialdir = "Config_Files/",title = "Select file",filetypes = (("CSV files","*.csv"),("all files","*.*")))

    # Parse CSV file and store output vars
    if not file_path:
        print("error: no file selected!")
    else:
        try :
            with open(file_path, mode='r') as csv_file:
                csv_reader = csv.DictReader(csv_file)
                line_count = 0
                for row in csv_reader:
                    if line_count == 0:
                        line_count += 1
                    trialTargetOrder = np.append(trialTargetOrder, int(row["TrialTargetNum"]))
                    blockMarkersByTrial = np.append(blockMarkersByTrial, int(row["BlockMarkersOn"]))
                    trialBlockNumber = np.append(trialBlockNumber, int(row["BlockNum"]))
                    trialTargVert = np.append(trialTargVert, int(row["TargVert"]))
                    calOrderID = np.append(calOrderID, int(row["CalOrdID"]))
                    line_count += 1
                print('Got ' + os.path.basename(file_path) + ', processed ' + str(line_count-1) + ' trials total.')
                
                blockNumbers = np.unique(trialBlockNumber)
                num_blocks = np.size(blockNumbers,0)

                for bn in range(num_blocks):
                    current_block = blockNumbers[bn]
                    block_markers = blockMarkersByTrial[trialBlockNumber == current_block]
                    if np.all(block_markers == 0):
                        blockMarkersOn = np.append(blockMarkersOn, 0)
                    elif np.all(block_markers == 1):
                        blockMarkersOn = np.append(blockMarkersOn, 1)
                    elif np.all(block_markers == 2):
                        blockMarkersOn = np.append(blockMarkersOn, 2)

                for bn in range(num_blocks):
                    current_block = blockNumbers[bn]
                    targ_vert = trialTargVert[trialBlockNumber == current_block]
                    if np.all(targ_vert == 0):
                        blockTargVert = np.append(blockTargVert, 0)
                    elif np.all(targ_vert == 1):
                        blockTargVert = np.append(blockTargVert, 1)
                    elif np.all(block_markers == 2):
                        blockTargVert = np.append(blockTargVert, 2)

                calOrderID = np.unique(calOrderID) 
                

                if calOrderID == 1:
                    calOrderID = ["H", "V"]
                else:
                    calOrderID = ["V", "H"]

        except:
            print("error: could not load file!")

    return trialTargetOrder, trialBlockNumber, trialTargVert, blockNumbers, blockMarkersOn, blockTargVert, calOrderID    



def load_cond2_config_file() :
    
    # Init output vars (to pass to ss_trial)
    trialTargetOrder = np.array([], dtype=int)
    blockEyesByTrial = np.array([], dtype=int)
    blockEyesOpen = np.array([], dtype=int)
    trialBlockNumber = np.array([], dtype=int)

    # Open system dialog to choose config file
    root = Tkinter.Tk()
    root.withdraw()
    file_path = tkFileDialog.askopenfilename(initialdir = "Config_Files/",title = "Select file",filetypes = (("CSV files","*.csv"),("all files","*.*")))

    # Parse CSV file and store output vars
    if not file_path:
        print("error: no file selected!")
    else:
        try :
            with open(file_path, mode='r') as csv_file:
                csv_reader = csv.DictReader(csv_file)
                line_count = 0
                for row in csv_reader:
                    if line_count == 0:
                        line_count += 1
                    trialTargetOrder = np.append(trialTargetOrder, int(row["TrialTargetNum"]))
                    blockEyesByTrial = np.append(blockEyesByTrial, int(row["BlockEyesOpen"]))
                    trialBlockNumber = np.append(trialBlockNumber, int(row["BlockNum"]))
                    line_count += 1
                print('Got ' + os.path.basename(file_path) + ', processed ' + str(line_count-1) + ' trials total.')
                
                blockNumbers = np.unique(trialBlockNumber)
                num_blocks = np.size(blockNumbers,0)
                for bn in range(num_blocks):
                    current_block = blockNumbers[bn]
                    block_markers = blockEyesByTrial[trialBlockNumber == current_block]
                    if np.all(block_markers == 0):
                        blockEyesOpen = np.append(blockEyesOpen, 0)
                    elif np.all(block_markers == 1):
                        blockEyesOpen = np.append(blockEyesOpen, 1)
        except:
            print("error: could not load file!")

    return trialTargetOrder, trialBlockNumber, blockNumbers, blockEyesOpen