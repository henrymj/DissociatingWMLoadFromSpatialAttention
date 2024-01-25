"""Dot cloud change etection experiment.

Author - Henry Jones henryjones@uchicago.edu

Adapted experiment code from 
Colin Quirk https://github.com/colinquirk/PsychopyChangeDetection 
Will Thyer https://github.com/henrymj/Cannonball/tree/master/experiment

Note: this code relies on Colin Quirk's templateexperiments module. You can get it from
https://github.com/colinquirk/templateexperiments and either put it in the same folder as this
code or give the path to psychopy in the preferences.
"""
import os
import glob
import sys
import errno
import pickle

import json
import random
import copy

import numpy as np
from numpy.matlib import repmat
import math

import psychopy.core
import psychopy.event
import psychopy.visual
import psychopy.parallel
import psychopy.tools.monitorunittools

import template 
import eyelinker

# PARAMATERS
exp_name = 'A01'
data_directory = os.path.join(
    '..', '..', '..', 'data', 'raw')

distance_to_monitor = 80  # cm
allowed_deg_from_fix = .75

# Trial parameters
n_trials_per_block = 96  # must be divisible by 2 for load counterbalancing
n_blocks = 16
n_trials_per_drift_correction = 6

data_fields = [
    'Subject',
    'Block',
    'Trial',
    'timestamp',
    'ITI',
    'RT',
    'Response',
    'ACC',
    'args_phase',
    'args_n_clouds',
    'args_loc_c1',
    'args_loc_c2',
    'args_width_c1',
    'args_width_c2',
    'args_interrupt',
    'memArray_x',
    'memArray_y',
    'memArray_sizes',
    'memArray_dotColors',
    'correctResponse',
    'change_type',
    'probeArray_x',
    'probeArray_y',
    'probeArray_sizes',
    'probeArray_dotColors',
    'interruptArray_x',
    'interruptArray_y',
    'interruptArray_sizes',
    'interruptArray_dotColors',
    'targetCloud_color',
    'targetCloud_bins',
    'targetCloud_nDots',
    'otherCloud_color',
    'otherCloud_bins',
    'otherCloud_nDots',
    'probeCloud_color',
    'probeCloud_bins',
    'probeCloud_nDots',
    'interruptCloud_color',
    'interruptCloud_bins',
    'interruptCloud_nDots',
    'memPort',
]


# Timing parameters (in seconds)
timing = {
    'minITI': .5,
    'maxITI': .9,
    'stepITI': .01,
    'baseline': .3,
    'memArray': .25,
    'delay': 1,
    'postBreak': 1,
    'feedback': 1,
    'interrupt_delay1': .5,
    'interruption': .15,
    'interrupt_delay2': .85,
    
}

# Stimulus details
# sizes are in degrees of visual angle
stimulus = {
    'fixSize': .15,
    'fixColor': 0,
    'regionHeight': 12,
    
    'minEccentricity': 1.75,
    'maxEccentricity': 5.25,
    'minDotSize': .25,
    'maxDotSize': .34,
    
    'nBins': 8,
    'nSlices': 40,
    'nRows': 10,
    'binSizes': [1, 3], # number of bins for a cloud to take up
    'nClouds': [1, 2],
    'minDots': 12,
    'max1CloudDots': 48,
    'max2CloudDots': 24,
    
    'nTrialsPractice': 10,
    'nTrialsPerBlock': n_trials_per_block,
    'nBlocks': n_blocks,
    'nTrials': n_trials_per_block*n_blocks,
    'nTrialsInterrupt': 192,
    'nBlocksInterrupt': 2,
    
    'colors': ['blue', 'green'], # stimulus colors
    
    'keyCodes': random.sample(['right', 'left'], 2)
}

stimulus['sliceThickness'] = 2*np.pi / stimulus['nSlices']
stimulus['nDotsPerBin'] = int(stimulus['nSlices']*stimulus['nRows']/stimulus['nBins'])
stimulus['binStarts'] = np.arange(0, stimulus['nDotsPerBin']*stimulus['nBins'], stimulus['nDotsPerBin'])

colors = {
    'fix': [1,1,1],
    'blue': [-1,-1,1],
    
    # convert 0 - 255 into -1 - 1
    'green': np.array([0,52,0])/127.5 - 1,
    'grey': np.array([42,42,42])/127.5 - 1  # for interruption trials
}

#    Trial codes: 1:96
#    Block codes: 101:~118 (+ makeup blocks)
#    Condition Codes: 211:253
#        200 +
#        1X = 1 cloud
#        2X = 2 clouds
#        X1 = cloud(s) are bin width 1
#        X2 = clouds are different widths
#        X3 = cloud(s) are bin width 3
#        + 30 = Interruption Trials
#    Probe: 170 (no change) & 171 (change)
portCodes = {
    'baseline': 150,
    'delay': 160,  #delay part 1 if interruption trial
    'interruption': 161,
    'delay2': 162,  # only if there is an interruption
    # whether we've made it to the probe part of the task
    'reject': 99,
    'success': 100, 

    'response': 180,

    'taskStart': 190,
    'taskEnd': 191
}



gender_options = [
    'Other/Choose Not To Respond',
    'Male',
    'Female',
]

hispanic_options = [
    'Choose Not To Respond',
    'Yes, Hispanic or Latino/a',
    'No, not Hispanic or Latino/a',
]

race_options = [
    'Choose Not To Respond',
    'American Indian or Alaskan Native',
    'Asian',
    'Pacific Islander',
    'Black or African American',
    'White / Caucasian',
    'More Than One Race',
]

questionaire_dict = {
    'Age': 99,
    'Gender': gender_options,
    'Hispanic:': hispanic_options,
    'Race': race_options
}

# necessary for saving out all trial info as json
class NpEncoder(json.JSONEncoder):  
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

# TASK CLASS
class Acacia01(template.BaseExperiment):
    """The class that runs the  experiment.

    Parameters:
    allowed_deg_from_fix -- The maximum distance in visual degrees the stimuli can appear from
        fixation
    colors -- The list of colors (list of 3 values, -1 to 1) to be used in the experiment.
    orients -- The list of orientsations to be used in the experiment.
    stim_idx -- List of indices for colors and orientations.
    data_directory -- Where the data should be saved.
    delay_time -- The number of seconds between the stimuli display and test.
    instruct_text -- The text to be displayed to the participant at the beginning of the
        experiment.
    iti_time -- The number of seconds in between a response and the next trial.
    keys -- The keys to be used for making a response. First is used for 'same' and the second is
        used for 'different'
    max_per_quad -- The number of stimuli allowed in each quadrant. If None, displays are
        completely random.
    min_distance -- The minimum distance in visual degrees between stimuli.
    n_blocks -- The number of blocks in the experiment.
    n_trials_per_block -- The number of trials within each block.
    percent_same -- A float between 0 and 1 (inclusive) describing the likelihood of a trial being
        a "same" trial.
    questionaire_dict -- Questions to be included in the dialog.
    sample_time -- The number of seconds the stimuli are on the screen for.
    set_sizes -- A list of all the set sizes. An equal number of trials will be shown for each set
        size.
    single_probe -- If True, the test display will show only a single probe. If False, all the
        stimuli will be shown.
    stim_size -- The size of the stimuli in visual angle.

    Additional keyword arguments are sent to template.BaseExperiment().

    Methods:
    chdir -- Changes the directory to where the data will be saved.
    display_break -- Displays a screen during the break between blocks.
    display_fixation -- Displays a fixation cross.
    display_stimuli -- Displays the stimuli.
    display_test -- Displays the test array.
    generate_locations -- Helper function that generates locations for make_trial
    get_response -- Waits for a response from the participant.
    make_block -- Creates a block of trials to be run.
    make_trial -- Creates a single trial.
    run_trial -- Runs a single trial.
    run -- Runs the entire experiment.
    """
    def __init__(self,
                n_trials_per_block=stimulus['nTrialsPerBlock'],
                n_blocks=stimulus['nBlocks'],
                n_trials_per_drift_correction=n_trials_per_drift_correction,
                timing=timing,
                stimulus=stimulus,
                colors=colors,
                portCodes=portCodes,
                allowed_deg_from_fix=allowed_deg_from_fix,
                data_directory=data_directory,
                questionaire_dict=questionaire_dict,
                track_eyes=True,
                **kwargs):

        self.n_trials_per_block = n_trials_per_block
        self.n_blocks = n_blocks
        self.n_trials_per_drift_correction = n_trials_per_drift_correction
        
        self.timing = timing
        self.stim = stimulus
        self.colors = colors
        self.portCodes=portCodes

        self.allowed_deg_from_fix = allowed_deg_from_fix

        self.data_directory = data_directory
        
        self.questionaire_dict = questionaire_dict
        
        self.track_eyes = track_eyes
        
        self.instruct_text = []
        self.instruct_interruption_text = []
        
        super().__init__(**kwargs)
        
        # Initialize eye tracking rejections
        self.rejections = []


    #  EYE TRACKING FUNCTIONS
    def init_tracker(self):
        eye_file = self.experiment_name + '_' + self.experiment_info['Subject Number'] + '.edf'
        
        if (len(glob.glob(eye_file)) > 0) and (not self.overwrite_ok): # ALWAYS CREATE NEW EYE TRACKING FILE IF NOT OVERWRITING
            i = 1
            new_eye_file = self.experiment_name + '_' + str(int(self.experiment_info['Subject Number'])+i) + '.edf' 
            while len(glob.glob(new_eye_file)) > 0:
                i += 1
                new_eye_file = self.experiment_name + '_' + str(int(self.experiment_info['Subject Number'])+i) + '.edf' 
            eye_file = new_eye_file

        if self.track_eyes:
            self.tracker = eyelinker.EyeLinker(self.experiment_window,eye_file,'BOTH')
        else:
            print('NOT TRACKING EYES')
            self.tracker = eyelinker.MockEyeLinker(self.experiment_window,eye_file,'BOTH')
                

        self.tracker.initialize_graphics()
        self.tracker.open_edf()
        self.tracker.initialize_tracker()
        self.tracker.send_tracking_settings()


    def show_eyetracking_instructions(self):
        self.tracker.display_eyetracking_instructions()
        self.tracker.setup_tracker()


    def start_eyetracking(self, block_num, trial_num):
        """Send block and trial status and start eyetracking recording

        Parameters:
        block_num-- Which block participant is in
        trial_num-- Which trial in block participant is in
        """
        status = f'Block {block_num}, Trial {trial_num}'
        self.tracker.send_status(status)

        self.tracker.start_recording()


    def stop_eyetracking(self):
        """Stop eyetracking recording
        """
        self.tracker.stop_recording()


    def realtime_eyetracking(self,wait_time,trial,sampling_rate=.01):
        """Collect real time eyetracking data over a period of time

        Returns eyetracking data

        Parameters:
        wait_time-- How long in s to collect data for
        sampling_rate-- How many ms between each sample
        """
        start_time = psychopy.core.getTime()
        while psychopy.core.getTime() < (start_time + wait_time):

            realtime_data = self.tracker.gaze_data

            reject,eyes = self.check_realtime_eyetracking(realtime_data)
            
            if reject:
                psychopy.core.wait(.01)
                self.send_synced_event(self.portCodes['reject'])

                self.rejections.append(trial['args'])
                
                print(f'# of rejected trials: {len(self.rejections)}')
                    
                self.stop_eyetracking()
                self.display_eyemovement_feedback(eyes)
                return reject
            psychopy.core.wait(sampling_rate)


    def check_realtime_eyetracking(self,realtime_data):
        left_eye,right_eye = realtime_data
        if left_eye:
            lx,ly = left_eye
        if right_eye:
            rx,ry = right_eye      
        if (not left_eye) & (not right_eye):
            return False,None

        eyex = np.nanmean([lx,rx])
        eyey = np.nanmean([ly,ry])
        
        winx,winy = self.experiment_window.size/2
        
        eyex -= winx
        eyey -= winy
        eyes = np.array([eyex,eyey])

        limit_radius = psychopy.tools.monitorunittools.deg2pix(self.allowed_deg_from_fix, self.experiment_monitor)

        euclid_distance = np.linalg.norm(eyes-np.array([0,0]))

        if euclid_distance > limit_radius:
            return True,(eyex,eyey)
        else:
            return False,None


    def display_eyemovement_feedback(self,eyes):

        psychopy.visual.TextStim(win=self.experiment_window,text='Eye Movement Detected.\nPress either key to continue.',pos = [0,1.2], color = [1,-1,-1]).draw()
        self.draw_fixation()
        
        psychopy.visual.Circle(win=self.experiment_window,radius=5,pos=eyes,fillColor='red',units='pix').draw()
        
        self.experiment_window.flip()
        
        psychopy.event.waitKeys(keyList=self.stim['keyCodes'])


    def kill_tracker(self):
        """
        Turns off eyetracker and transfers EDF file
        """
        self.tracker.set_offline_mode()
        self.tracker.close_edf()
        self.tracker.transfer_edf()
        self.tracker.close_connection()


    #  EEG FUNCTIONS
    def setup_eeg(self):
        """
        Connects the parallel port for EEG port code
        """
        try:
            self.port = psychopy.parallel.ParallelPort(address=53328)
            self.send_synced_event = self.send_synced_event_wPort
            self.send_synced_event(self.portCodes['taskStart'])

        # _ is raised when tested on windows and no port is found
        # NotImplementedError is raised when tested on mac
        except (NotImplementedError) as e:
            self.port = None
            self.send_synced_event = self.send_synced_event_nullPort
            print('No parallel port connected. Port codes will not send!')
            print('Exact error: ', e)


    def send_synced_event_wPort(self, code, keyword = "SYNC"):
        """Send port code to EEG and eyetracking message for later synchronization

        Parameters:
        code-- Digits to send
        keyword-- Accompanying sync keyword (matters for later EEGLAB preprocessing)
        """
        self.tracker.send_message(keyword + ' ' + str(code))
        self.port.setData(code)
        

    def send_synced_event_nullPort(self, *args, **kwargs): # used in case of mac issues
        pass


    # TASK & STIMULUS FUNCTIONS
    def chdir(self, dir=None):
        """Changes the directory to where the data will be saved.
        """
        dir = self.data_directory if dir is None else dir
        os.makedirs(dir, exist_ok=True)
        os.chdir(dir)


    def make_practice_trials(self, n_trials_per_block):

        trials = []
        for _ in range(n_trials_per_block):
            trials.append(
                self.make_trial(
                    random.choice(np.arange(0, self.stim['nBins'])), # loc_c1
                    random.choice(np.arange(0, self.stim['nBins'])), # loc_c2
                    random.choice(self.stim['binSizes']), # width_c1
                    random.choice(self.stim['binSizes']), # width_c2
                    random.choice(self.stim['nClouds']), # n_clouds
                    phase='practice'
                )
            )
        return trials


    def make_trials(self, CB='full', **kwargs):
        """
        function for generating trials.
        CB - counterbalance style. If 'full', then balance locations and widths along with load. Used for main test trials
            if 'Load', only balance along load - used for practice and interruption trials.
        """

        if CB=='load':
            return self.make_trials_loadCB(**kwargs)
        
        # how many trials to covers a full permutation our variables of interest, how many full permutations fits into nTrials
        n_per_perm = self.stim['nBins']**2 * len(self.stim['binSizes'])**2 * len(self.stim['nClouds']) 
        n_loops = int(self.stim['nTrials'] / n_per_perm)
        
        if n_loops < 1:
            print('not enough nTrials to create a full permutation - generating randomly w/ balanced load:')
            return(self.make_trials_loadCB(self.stim['nTrials'], **kwargs))
        else:
            assert int(n_loops)==n_loops # make sure it's a nice whole number
            
        return self.make_trials_fullCB(n_loops, **kwargs)

    def make_trials_fullCB(self, n_loops, phase=None, interrupt=False):
        """
        Function to generate n_loops of fully counter-balanced trials (locations, widths, and load) 
        """
        trials = {'load1': [], 'load2': []}
        for _ in range(n_loops):
            for loc_c1 in range(self.stim['nBins']):
                for loc_c2 in range(self.stim['nBins']):
                    for width_c1 in self.stim['binSizes']:
                        for width_c2 in self.stim['binSizes']:
                            for n_clouds in self.stim['nClouds']:
                                trials[f'load{n_clouds}'].append(self.make_trial(loc_c1, loc_c2, width_c1, width_c2, n_clouds, interrupt=interrupt, phase=phase))

        # Permute the order of locations and widths
        np.random.shuffle(trials['load1'])
        np.random.shuffle(trials['load2'])
        return trials

    def make_trials_loadCB(self, n_trials, phase=None, interrupt=False):
        """
        Catch all function for making practice trials, interruption trials, and random test trials
        Always has load counterbalanced
        """
        trials = {'load1': [], 'load2': []}
        n_per_load = n_trials/len(self.stim['nClouds'])
        assert int(n_per_load)==n_per_load # make sure it's a nice whole number
        
        n_clouds_order = random.sample(self.stim['nClouds'] * int(n_per_load), n_trials)
        for n_clouds in n_clouds_order:
            trials[f'load{n_clouds}'].append(
                self.make_trial(
                    random.choice(np.arange(0, self.stim['nBins'])), # loc_c1
                    random.choice(np.arange(0, self.stim['nBins'])), # loc_c2
                    random.choice(self.stim['binSizes']), # width_c1
                    random.choice(self.stim['binSizes']), # width_c2
                    n_clouds, # n_clouds
                    interrupt=interrupt,
                    phase=phase
                )
            )
        return trials


    def make_trial(self, loc_c1, loc_c2, width_c1, width_c2, n_clouds, interrupt=False, phase=None):
        trial = {}
        trial['args'] = {
            'phase': phase,
            'n_clouds': n_clouds,
            'interrupt': interrupt,
            'loc_c1': loc_c1,
            'loc_c2': loc_c2,
            'width_c1': width_c1,
            'width_c2': width_c2
        }

        colors = random.sample(self.stim['colors'], 2)
        targetCloud = {
            'color': colors[0],
            'bins': self.get_bin_locs(loc_c1, width_c1)
        }
        otherCloud = {
            'color': colors[1],
            'bins': self.get_bin_locs(loc_c2, width_c2)
        }

        if n_clouds==1:
            targetCloud['nDots'] = np.random.choice(np.arange(self.stim['minDots'], self.stim['max1CloudDots']+1))

            [x, y, sizes, dotColors] = self.gen_dot_info(targetCloud, None);

            # set otherCloud to basically empty
            otherCloud['nDots'] = 0
            otherCloud['bins'] = []
            otherCloud['idx'] = []

            trialCodeA = 10
            if width_c1==1:
                trialCodeB = 1
            else:
                trialCodeB = 3

        elif n_clouds==2:
            targetCloud['nDots'] = np.random.choice(np.arange(self.stim['minDots'], self.stim['max2CloudDots']+1))
            otherCloud['nDots'] = np.random.choice(np.arange(self.stim['minDots'], self.stim['max2CloudDots']+1))

            [x, y, sizes, dotColors] = self.gen_dot_info(targetCloud, otherCloud);

            trialCodeA = 20
            if (width_c1==1) and (width_c2==1):
                trialCodeB = 1
            elif (width_c1==3) and (width_c2==3):
                trialCodeB = 3
            else:
                trialCodeB = 2

        trial['memArray'] = {
            'x': x,
            'y': y,
            'sizes': sizes,
            'dotColors': dotColors
        }

        # Generate Probe Display
        probeCloud = targetCloud.copy()

        change = random.choice([True, False])
        if change:
            trial['correctResponse'] = self.stim['keyCodes'][0]

            change_type = random.choice([1, 2])
            trial['change_type'] = change_type

            if change_type==1:  # shift location by 1 bin
                shift = random.choice([1, -1])
                trial['probe_bin_shift'] = shift
                angles = np.arctan2(y, x)
                angles = angles + shift*(2*np.pi)/self.stim['nBins']
                dists = np.sqrt(x**2 + y**2)
                x = np.cos(angles)*dists
                y = np.sin(angles)*dists

                # update bins, somewhat as a formality
                old_bins = probeCloud['bins']
                bins = old_bins + shift
                bins = (bins + self.stim['nBins']) % self.stim['nBins']  # -1 becomes 7
                probeCloud['bins'] = bins

            elif change_type==2:  # change width, broad <-> narrow
                old_bins = probeCloud['bins'];
                if len(old_bins)==3:
                    bins = np.array([old_bins[1]])
                else:
                    b = old_bins[0]
                    bins = np.array([b - 1, b, b+1])

                bins = (bins + self.stim['nBins']) % self.stim['nBins']  # -1 becomes 7
                probeCloud['bins'] = bins

                (x, y, sizes, dotColors) = self.gen_dot_info(probeCloud);
        else:
            trial['correctResponse'] = self.stim['keyCodes'][1]
            trial['change_type'] = 0

        trial['probeArray'] = {
            'x': x[:probeCloud['nDots']],
            'y': y[:probeCloud['nDots']],
            'sizes': sizes[:probeCloud['nDots']],
            'dotColors': dotColors[:probeCloud['nDots'], :]
        }
        
        
        if interrupt:  # make a random cloud in a random location
            interruptCloud = {
                'color': 'grey',
                'bins': self.get_bin_locs(
                    random.choice(np.arange(0, self.stim['nBins'])),
                    random.choice(self.stim['binSizes'])
                ),
                'nDots': np.random.choice(np.arange(self.stim['minDots'], self.stim['max1CloudDots']+1))
            }
            
            (x, y, sizes, dotColors) = self.gen_dot_info(interruptCloud);
            
            trial['interruptArray'] = {
                'x': x,
                'y': y,
                'sizes': sizes,
                'dotColors': dotColors
            }
        else:
            interruptCloud = {
                'color': 'grey',
                'bins': [],
                'nDots': 0
            }
            trial['interruptArray'] = {
                'x': np.array([]),
                'y': np.array([]),
                'sizes': np.array([]),
                'dotColors': np.array([])
            }
            
        trial['targetCloud'] = targetCloud
        trial['otherCloud'] = otherCloud
        trial['probeCloud'] = probeCloud
        trial['interruptCloud'] = interruptCloud
        
        trial['memPort'] = 200 + trialCodeA + trialCodeB + 30*interrupt

        return(trial)


    def gen_dot_info(self, cloud1, cloud2=None):

        n_slices = self.stim['nSlices']
        n_rows = self.stim['nRows']
        min_dist = self.stim['minEccentricity']
        max_dist = self.stim['maxEccentricity']
        row_thickness = (max_dist - min_dist)/n_rows
        slice_thickness = stimulus['sliceThickness']
        min_size = self.stim['minDotSize']
        max_size = self.stim['maxDotSize']

        row_sizes = np.linspace(min_size, max_size, n_rows)
        sizes = repmat(row_sizes, n_slices,1).T

        # set the distances from center with jitter
        initial_distances = min_dist + .5*row_thickness + np.arange(0, n_rows*row_thickness, row_thickness)  # centers of each row
        distance_jitter = (np.random.uniform(size=(n_rows, n_slices))-.5)*(row_thickness-sizes)
        distances = repmat(initial_distances, n_slices, 1).T + distance_jitter

        # set angles with jitter
        initial_angles = np.linspace(0, 2*np.pi-slice_thickness, n_slices) + slice_thickness/2 # centers of each slice
        # with jitter, adjust for the distance and size of the dot
        angle_jitter = (np.random.uniform(size=(n_rows, n_slices))-.5)*(slice_thickness - (sizes/distances)) 
        angles = repmat(initial_angles, n_rows, 1) + angle_jitter

        x = np.cos(angles)*distances
        y = np.sin(angles)*distances

        keep_idx, protected_idx = self.gen_cloud_idx(cloud1)
        dot_colors = repmat(self.colors[cloud1['color']], cloud1['nDots'], 1)

        if cloud2 is not None:
            cloud2_idx, _ = self.gen_cloud_idx(cloud2, keep_idx, protected_idx)
            keep_idx = np.concatenate((keep_idx, cloud2_idx))
            dot_colors = np.concatenate((dot_colors, repmat(self.colors[cloud2['color']], cloud2['nDots'], 1)))


        x = x.flatten(order='F')[keep_idx]
        y = y.flatten(order='F')[keep_idx]
        sizes = sizes.flatten(order='F')[keep_idx]

        return(x, y, sizes, dot_colors)


    def gen_cloud_idx(self, cloud, other_cloud_idx=None, pre_protected_idx=None):
        '''
        Take in cloud, generates indices for the appropriate bins.
        
        If the first cloud, always uses an index from the first and last column/wedge, so that it spans the full bin.
        
        If the second cloud, checks for overlap with the first cloud before deciding on index in first and last column/wedge. Avoids indices picked by the first cloud.
        '''
        pre_protected_idx = pre_protected_idx if pre_protected_idx is not None else np.array([])
        
        nBins = len(cloud['bins'])
        possible_idx = np.full(int(nBins*self.stim['nDotsPerBin']), np.nan).astype(int)

        # Get indices that the cloud can't use, so that the other cloud also always be at the edge of any bin that it's in
        protected_first_idx = np.full(nBins, np.nan).astype(int)
        protected_last_idx = np.full(nBins, np.nan).astype(int)

        for bini in range(nBins):
            bin_start_idx = self.stim['binStarts'][cloud['bins'][bini]]
            curr_idx = np.arange(bin_start_idx, bin_start_idx+self.stim['nDotsPerBin'])
            possible_idx[((bini)*self.stim['nDotsPerBin']):((bini+1)*self.stim['nDotsPerBin'])] = curr_idx

            # note an idx in the first and last slice of each bin
            # to reserve for the second cloud
            protected_first_idx[bini] = random.choice(curr_idx[:self.stim['nRows']])
            protected_last_idx[bini] = random.choice(curr_idx[-self.stim['nRows']:])

        protected_idx = np.concatenate((protected_first_idx, protected_last_idx))

        cloud_idx = np.full(cloud['nDots'], np.nan).astype(int)

        # set indices, taking into account protected indices and/or the other cloud
        if other_cloud_idx is None:
            # First cloud:
            # Must use a dot in the first and last slice of the
            # bins it spans. Those dots must not come from the protected set.
            # Remaining dots are sampled uniformly, avoiding the protected set.
            cloud_idx[0] = random.choice(np.setdiff1d(possible_idx[1:self.stim['nRows']], protected_idx))
            cloud_idx[1] = random.choice(np.setdiff1d(possible_idx[-self.stim['nRows']:], protected_idx))
            cloud_idx[2:] = np.random.choice(np.setdiff1d(possible_idx, np.concatenate((protected_idx, cloud_idx[:2]))), size=cloud['nDots']-2, replace=False)
        else: 
            # Second cloud:
            # For first and last slice, check for a reserved dot.
            # Use if it exists.
            # If if does not exist, sample from the first and last slice
            # indices, while avoiding those used by the first cloud.
            # for the remaining dots, sample uniformly while avoiding those
            # used by the first dot cloud.

            # check for a reserved dot in the first slice, use if it exists
            tmp_intersection = np.intersect1d(possible_idx[:self.stim['nRows']], pre_protected_idx)
            if len(tmp_intersection):
                assert(len(tmp_intersection)==1)
                cloud_idx[0] = tmp_intersection[0]
            else:
                cloud_idx[0] = random.choice(np.setdiff1d(possible_idx[1:self.stim['nRows']], protected_idx))

            # check for a reserved dot in the first slice, use if it exists
            tmp_intersection = np.intersect1d(possible_idx[-self.stim['nRows']:], pre_protected_idx)
            if len(tmp_intersection):
                assert(len(tmp_intersection)==1)
                cloud_idx[1] = tmp_intersection[0]
            else:
                cloud_idx[1] = random.choice(np.setdiff1d(possible_idx[-self.stim['nRows']:], protected_idx))

            cloud_idx[2:] = np.random.choice(np.setdiff1d(possible_idx, np.concatenate((other_cloud_idx, cloud_idx[:2]))), size=cloud['nDots']-2, replace=False)

        return(cloud_idx, protected_idx)


    def get_bin_locs(self, loc, width):
        bins = loc + np.arange(0, width)
        if width==3:  # center bins at loc
            bins -= 1  
        bins = (bins + self.stim['nBins']) % self.stim['nBins']  # -1 -> 7; 8 -> 0
        return bins


    # DISPLAY FUNCTIONS
    def init_fixation(self):
        # set up fixation shapes, based on Thaler et al., 2013
        self.fix_outer = psychopy.visual.Circle(
            self.experiment_window, lineColor=None, fillColor = [-1,-1,-1], 
            fillColorSpace='rgb', radius=.25, units='deg', interpolate=True
        )
        self.fix_cross = psychopy.visual.ShapeStim(
            self.experiment_window,
            vertices=((0, -0.26), (0, 0.26), (0,0), (-0.26,0), (0.26, 0)),
            lineWidth=psychopy.tools.monitorunittools.deg2pix(.16, self.experiment_monitor),
            closeShape=False,
            lineColor=[0,0,0]
        )
        self.fix_inner = psychopy.visual.Circle(
            self.experiment_window, lineColor=None, fillColor = [-1,-1,-1], 
            fillColorSpace='rgb', radius=.075, units='deg', interpolate=True
        )


    def draw_fixation(self):
        self.fix_outer.draw()
        self.fix_cross.draw()
        self.fix_inner.draw()


    def display_fixation(self, wait_time = None, text = None, keyList = None, realtime_eyetracking = False, trial = None):
        """Displays a fixation cross. A helper function for self.run_trial.

        Parameters:
        wait_time -- The amount of time the fixation should be displayed for.
        text -- Str that displays above fixation cross. 
        keyList -- If keyList is given, will wait until key press
        trial -- Trial object needed for realtime eyetracking functionality.
        real_time_eyetracking -- Bool for if you want to do realtime eyetracking or not
        """
        
        if text is not None:
            psychopy.visual.TextStim(win=self.experiment_window,text=text,pos = [0,1], color = [1,-1,-1]).draw()

        self.draw_fixation()
        
        self.experiment_window.flip()

        if realtime_eyetracking:
            reject = self.realtime_eyetracking(wait_time=wait_time, trial=trial)
            return reject    
        else:
            if keyList:
                resp = psychopy.event.waitKeys(maxWait=wait_time,keyList=keyList)
                if resp == ['p']:
                    self.display_text_screen(text='Paused',keyList = ['space'])
                    self.display_fixation(wait_time=1)
                elif resp == ['e']:
                    self.tracker.calibrate()
                    self.display_fixation(wait_time=1)
                elif resp == ['escape']:
                    resp = self.display_text_screen(text='Are you sure you want to exit?',keyList = ['y','n'])
                    if resp == ['y']:
                        psychopy.core.wait(1)
                        self.send_synced_event(self.portCodes['taskEnd'])
                        self.tracker.transfer_edf()
                        self.quit_experiment()
                    else:
                        self.display_fixation(wait_time=1)
            else:
                psychopy.core.wait(wait_time)

    def init_stimtrak(self, x=930, y=510):
        # not exactly fixation, but a stable thing we want to reuse
        self.stimtrak = psychopy.visual.Circle(
            self.experiment_window, lineColor=None, fillColor = [1,1,1], 
            fillColorSpace='rgb', radius=20, pos = [x,y], units='pix'
        )

    def draw_stimtrak(self):
        self.stimtrak.draw()


    def display_block_feedback(self,acc):

        text = (
            'Please take a short break.\n'
            f'Your accuracy this block was: {round(100*np.nanmean(acc))}%\n\n'
            'Remember:\n'
            f'If the set has changed, press the {self.stim["keyCodes"][0]} key.\n\n'
            f'If the set has not changed, press the {self.stim["keyCodes"][1]} key.\n\n'
            'When you are ready to begin the next block, press either key.'
        )
        self.display_text_screen(
                text=text,
                keyList=self.stim['keyCodes'])


    def display_memArray(self, trial, realtime_eyetracking=False):
        """Displays the stimuli. A helper function for self.run_trial.

        Parameters:
        locations -- A list of locations (list of x and y value) describing where the stimuli
            should be displayed.
        colors -- A list of colors describing what should be drawn at each coordinate.
        """
        self.draw_fixation()
        
        dots = psychopy.visual.ElementArrayStim(
            win=self.experiment_window,
            nElements=len(trial['memArray']['x']),
            xys=np.concatenate((trial['memArray']['x'][:, np.newaxis], trial['memArray']['y'][:, np.newaxis]), 1),
            units="deg",
            sizes=trial['memArray']['sizes'],
            colors=trial['memArray']['dotColors'],
            elementTex=None,
            elementMask="circle"
        )
        
        dots.draw()
        self.draw_stimtrak()
        self.experiment_window.flip()

        if realtime_eyetracking:  
            reject = self.realtime_eyetracking(wait_time=self.timing['memArray'], trial=trial)
            return reject
        else:
            psychopy.core.wait(self.timing['memArray'])


    def display_interruption(self, trial, realtime_eyetracking=False):
        """Displays the stimuli. A helper function for self.run_trial.

        Parameters:
        locations -- A list of locations (list of x and y value) describing where the stimuli
            should be displayed.
        colors -- A list of colors describing what should be drawn at each coordinate.
        """

        self.draw_fixation()
        
        dots = psychopy.visual.ElementArrayStim(
            win=self.experiment_window,
            nElements=len(trial['interruptArray']['x']),
            xys=np.concatenate((trial['interruptArray']['x'][:, np.newaxis], trial['interruptArray']['y'][:, np.newaxis]), 1),
            units="deg",
            sizes=trial['interruptArray']['sizes'],
            colors=trial['interruptArray']['dotColors'],
            elementTex=None,
            elementMask="circle"
        )
        
        dots.draw()
        self.experiment_window.flip()

        if realtime_eyetracking:  
            reject = self.realtime_eyetracking(wait_time=self.timing['interruption'], trial=trial)
            return reject
        else:
            psychopy.core.wait(self.timing['interruption'])


    def get_response_from_probe(self, trial, keyList=None):
        
        self.send_synced_event(self.portCodes['success']) #success code means no eyetracking rejection
        
        self.draw_fixation()
        
        dots = psychopy.visual.ElementArrayStim(
            win=self.experiment_window,
            nElements=len(trial['probeArray']['x']),
            xys=np.concatenate((trial['probeArray']['x'][:, np.newaxis], trial['probeArray']['y'][:, np.newaxis]), 1),
            units="deg",
            sizes=trial['probeArray']['sizes'],
            colors=trial['probeArray']['dotColors'],
            elementTex=None,
            elementMask="circle"
        )
        
        dots.draw()

        self.experiment_window.flip()
        rt_timer = psychopy.core.MonotonicClock()
        resp = psychopy.event.waitKeys(keyList=keyList)
        rt = rt_timer.getTime()*1000

        psychopy.core.wait(.01) # wait 10ms to avoid dropping port code - only helpful is participant accidentally responds _super_ quickly
        self.send_synced_event(self.portCodes['response'])
        
        # we are using their first response
        return resp[0], rt


    def send_data(self, data):
        """Updates the experiment data with the information from the last trial.

        This function is seperated from run_trial to allow additional information to be added
        afterwards.

        Parameters:
        data -- A dict where keys exist in data_fields and values are to be saved.
        """
        self.update_experiment_data([data])


    def run_trial(self, trial, block_num, trial_num, realtime_eyetracking=False):
        """Runs a single trial.

        Returns the data from the trial after getting a participant response.

        Parameters:
        trial -- The dictionary of information about a trial.
        block_num -- The number of the block in the experiment.
        trial_num -- The number of the trial within a block.
        """

        self.send_synced_event(trial_num)

        # ITI
        ITI = random.choice(np.arange(timing['minITI'], timing['maxITI'], timing['stepITI']))
        self.display_fixation(
            wait_time=ITI,
            trial=trial,
            keyList=['p','escape','e'],
            realtime_eyetracking=False # do not care about eye movements during baseline
            )
        
        self.start_eyetracking(block_num = block_num, trial_num = trial_num)
        
        # BASELINE
        self.send_synced_event(self.portCodes['baseline'])
        reject = self.display_fixation(
            wait_time=self.timing['baseline'],
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None

        # MEMORY ARRAY
        self.send_synced_event(trial['memPort'])
        reject = self.display_memArray(
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None
        
        # DELAY
        self.send_synced_event(self.portCodes['delay'])
        reject = self.display_fixation(
            self.timing['delay'],
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None

        # PROBE
        resp, rt = self.get_response_from_probe(trial, keyList=self.stim['keyCodes'])

        self.stop_eyetracking()
        
        acc = 1 if resp == trial['correctResponse'] else 0

        data = {
            'Subject': self.experiment_info['Subject Number'],
            'Block': block_num,
            'Trial': trial_num,
            'Timestamp': psychopy.core.getAbsTime(),
            'ITI': ITI,
            'RT': rt,
            'Response': resp,
            'ACC': acc,
        }
        # add in trial info, flattened
        for k in trial:
            if type(trial[k]) is dict:
                for k2 in trial[k]:
                    data[f'{k}_{k2}'] = trial[k][k2]
            else:
                data[k] = trial[k]
        
        print(f'B: {block_num}, T: {trial_num}, ACC: {acc}, RT: {rt: 2f}')
        return data


    def run_trial_wInterruption(self, trial, block_num, trial_num, realtime_eyetracking=False):
        """Runs a single trial.

        Returns the data from the trial after getting a participant response.

        Parameters:
        trial -- The dictionary of information about a trial.
        block_num -- The number of the block in the experiment.
        trial_num -- The number of the trial within a block.
        """
        self.send_synced_event(trial_num)
        

        # ITI
        ITI = random.choice(np.arange(timing['minITI'], timing['maxITI'], timing['stepITI']))
        self.display_fixation(
            wait_time=ITI,
            trial=trial,
            keyList=['p','escape','e'],
            realtime_eyetracking=False
            )
        
        self.start_eyetracking(block_num = block_num, trial_num = trial_num)
        
        # BASELINE
        self.send_synced_event(self.portCodes['baseline'])
        reject = self.display_fixation(
            wait_time=self.timing['baseline'],
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None

        # MEMORY ARRAY
        self.send_synced_event(trial['memPort'])
        reject = self.display_memArray(
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None
        
        # 1ST DELAY - 500 ms
        self.send_synced_event(self.portCodes['delay'])
        reject = self.display_fixation(
            self.timing['interrupt_delay1'],
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None
        
        # INTERRUPTION - 150ms
        self.send_synced_event(self.portCodes['interruption'])
        reject = self.display_interruption(
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None
        
        # 2ND DELAY
        self.send_synced_event(self.portCodes['delay2'])
        reject = self.display_fixation(
            self.timing['interrupt_delay2'],
            trial=trial,
            realtime_eyetracking=realtime_eyetracking
            )
        if reject:
            return None

        # PROBE
        resp, rt = self.get_response_from_probe(trial, keyList=self.stim['keyCodes'])
        
        self.stop_eyetracking()
        
        acc = 1 if resp == trial['correctResponse'] else 0

        data = {
            'Subject': self.experiment_info['Subject Number'],
            'Block': block_num,
            'Trial': trial_num,
            'Timestamp': psychopy.core.getAbsTime(),
            'ITI': ITI,
            'RT': rt,
            'Response': resp,
            'ACC': acc,
        }
        # add in trial info, flattened
        for k in trial:
            if type(trial[k]) is dict:
                for k2 in trial[k]:
                    data[f'{k}_{k2}'] = trial[k][k2]
            else:
                data[k] = trial[k]

        print(f'B: {block_num}, T: {trial_num}, ACC: {acc}, RT: {rt: 2f}')
        return data


    def retrieve_block_from_trials(self, trials, n_trials=None):
        """
        Get blocks with balanced load, return block and remaining trials
        """
        n_trials = self.n_trials_per_block if n_trials is None else n_trials
        
        assert (n_trials % 2) == 0
        
        # get block with balanced load
        block = trials['load1'][:int(n_trials/2)] + trials['load2'][:int(n_trials/2)]
        np.random.shuffle(block)
        
        # trim off this block's trials from the remaining ones
        trials['load1'] = trials['load1'][int(n_trials/2):]
        trials['load2'] = trials['load2'][int(n_trials/2):]
        
        return block, trials


    def run_makeup_blocks(self, block_num, realtime_eyetracking=False, interrupt=False, remaining_info_filename=None):
        
        if interrupt:
            run_trial = self.run_trial_wInterruption
        else:
            run_trial = self.run_trial
            
        ### LOAD BALANCED VERSION
#        # generate makeup trials
#        makeup_trials = {'load1': [], 'load2': []}
#        for t_args in self.rejections:
#            t_args['phase'] = t_args['phase'] + '_makeup'
#            makeup_trials[f'load{t_args["n_clouds"]}'].append(self.make_trial(**t_args))  # splitting by load
#            
#        n_load1 = len(makeup_trials['load1'])
#        n_load2 = len(makeup_trials['load2'])
#        
#        while (n_load1 > 0) and (n_load2 > 0):
#            if (n_load1 > self.n_trials_per_block/2) and (n_load2 > self.n_trials_per_block/2): # run full counterbalanced block
#                block, makeup_trials = self.retrieve_block_from_trials(makeup_trials)
#            else:
#                min_n = min(n_load1, n_load2)
#                if min_n > 0:
#                    block, makeup_trials = self.retrieve_block_from_trials(makeup_trials, n_trials=min_n)
#                    
#            n_load1 = len(makeup_trials['load1']) # getting updated vals for next loop
#            n_load2 = len(makeup_trials['load2'])

        rejections = self.rejections.copy()
        self.rejections = []
        np.random.shuffle(rejections)
        
        while len(rejections) > 0:
            print(f'# MAKEUP TRIALS REMAINING: {len(rejections)}')
            block_num += 1
            assert block_num > 0

            if len(rejections) > self.n_trials_per_block:
                block = []
                for t_args in rejections[:self.n_trials_per_block]:
                    block.append(self.make_trial(**t_args))
                rejections = rejections[self.n_trials_per_block:]
            else:
                block = []
                for t_args in rejections:
                    block.append(self.make_trial(**t_args))
                rejections = []
            
            acc = []
            self.tracker.calibrate()

            self.send_synced_event(100+block_num)
            psychopy.core.wait(.01)
            for trial_num, trial in enumerate(block):
                if trial_num % self.n_trials_per_drift_correction == 0:
                    self.tracker.drift_correct()
                    
                data = run_trial(trial, block_num, trial_num+1, realtime_eyetracking=realtime_eyetracking)
                if data:
                    self.send_data(data)
                    acc.append(data['ACC'])
                
            self.save_data_to_csv()
            
            if remaining_info_filename is not None:  # save out remaining trials in case of crash
                info_to_pickle = {'block_num': block_num, 'trials': {'load1': [], 'load2': []}, 'rejections': rejections}
                with open(remaining_info_filename, 'wb+') as f:
                    pickle.dump(info_to_pickle,f)
                
            self.display_block_feedback(acc)
            
        return block_num


    def run(self, realtime_eyetracking=True):
        """Runs the entire experiment.

        This function takes a number of hooks that allow you to alter behavior of the experiment
        without having to completely rewrite the run function. While large changes will still
        require you to create a subclass, small changes like adding a practice block or
        performance feedback screen can be implimented using these hooks. All hooks take in the
        experiment object as the first argument. See below for other parameters sent to hooks.

        Parameters:
        setup_hook -- takes self, executed once the window is open.
        before_first_trial_hook -- takes self, executed after instructions are displayed.
        pre_block_hook -- takes self, block list, and block num
            Executed immediately before block start.
            Can optionally return an altered block list.
        pre_trial_hook -- takes self, trial dict, block num, and trial num
            Executed immediately before trial start.
            Can optionally return an altered trial dict.
        post_trial_hook -- takes self and the trial data, executed immediately after trial end.
            Can optionally return altered trial data to be stored.
        post_block_hook -- takes self, executed at end of block before break screen (including
            last block).
        end_experiment_hook -- takes self, executed immediately before end experiment screen.
        """

        """
        Setup and Instructions
        """
        self.chdir()  # CD TO EXPERIMENT DIRECTORY TOP LEVEL

        ok = self.get_experiment_info_from_dialog(self.questionaire_dict)
        if not ok:
            print('Experiment has been terminated.')
            sys.exit(1)

        self.init_experiment_fileBase() # look for already existing files and increment if necessary
        print(self.fileBase)
        self.chdir(self.fileBase)  # CD TO SUB SPECIFIC DIRECTORY
        
        if not self.extend_ok:
            self.save_experiment_info()
            self.save_experiment_pickle(additional_fields_dict=self.stim)
        else:
            with open(self.fileBase+'.pickle', 'rb') as f:
                old_pkl = pickle.load(f)
            self.stim['keyCodes'] = old_pkl['keyCodes']
            
        self.open_csv_data_file()
            
        # init instructions with keycodes
        self.instruct_text = [(
            'Thank you for participating in this study!\n\n'
            'In this task, you will see 1 or 2 sets of colored dots - 1 blue and/or 1 green - spread around a central point.\n'
            'Afer a delay, one set will reappear.\n\n'
            f'If the set has changed (e.g. in the location or spread of dots), press the {self.stim["keyCodes"][0]} key.\n\n'
            f'If the set has not changed, press the {self.stim["keyCodes"][1]} key.\n\n'
            'We will begin with a brief practice. Press either key to begin.'
        )]

        self.instruct_interruption_text = [(
            'This task is almost identical to the previous task.\n\n'
            'However, during the delay, a set of gray dots will briefly appear. Please ignore the gray dots, and keep your eyes centered.\n'
            'Remember:\n'
            f'If the set has changed (e.g. in the location or spread of dots), press the {self.stim["keyCodes"][0]} key.\n\n'
            f'If the set has not changed, press the {self.stim["keyCodes"][1]} key.\n\n'
            'Press either key to begin.'
        )]

        self.open_window(screen=0)
        self.display_text_screen('Loading...', wait_for_input=False)

        self.init_fixation()
        self.init_stimtrak()
        self.init_tracker()

        for instruction in self.instruct_text:
            self.display_text_screen(text=instruction, keyList=self.stim['keyCodes'])
        
        self.show_eyetracking_instructions()

        """
        Practice
        """
        block_num = -1
        self.port = None
        self.send_synced_event = self.send_synced_event_nullPort
        prac = self.display_text_screen(text = f'Practice block?', keyList=['y','n'])
        
        while prac == ['y']:
            
            practice_trials = self.make_trials(CB='load', n_trials=self.stim['nTrialsPractice'], phase='practice')
            block = practice_trials['load1'] + practice_trials['load2']
            np.random.shuffle(block)
            acc = []
            
            for trial_num, trial in enumerate(block):
                data = self.run_trial(trial,block_num,trial_num+1)      
                if data:
                    self.send_data(data)
                    
                    feedback='Correct!' if data['ACC']==1 else 'Incorrect.'
                    textObject = psychopy.visual.TextStim(
                        self.experiment_window, text=feedback, color=[-1,-1,-1])
                    textObject.draw()
                    self.experiment_window.flip()
                    psychopy.core.wait(self.timing['feedback'])
                    
                    acc.append(data['ACC'])

            block_num -= 1  # blocks get more negative as practice goes on
            
            self.save_data_to_csv()
            self.display_block_feedback(acc)
        
            prac = self.display_text_screen(text = f'Practice block?', keyList=['y','n'])


        """
        Main Task
        """
        self.setup_eeg() # turns on eeg and sends start code
        psychopy.core.wait(1)
        
        # load remaining trials and rejections, if they exist, else generate new
        remaining_info_file =f'TMP_{self.fileBase}_remaining_info.pickle'

        if os.path.exists(remaining_info_file) and self.extend_ok:
            with open(remaining_info_file, 'rb+') as f:
                remaining_info = pickle.load(f)
            block_num = remaining_info['block_num'] # says how many blocks have been completed
            trials = remaining_info['trials']
            self.rejections = remaining_info['rejections']
        else:
            block_num = 0  # block_num tells you the test/interruption block you're in (when in the middle of a block), or the number of test/interruption blocks you've completed (in this case, 0)
            trials = self.make_trials(CB='full', phase='test')
            self.rejections = []
            info_to_pickle = {'block_num': block_num, 'trials': trials, 'rejections': self.rejections}
            with open(remaining_info_file, 'wb+') as f:
                pickle.dump(info_to_pickle,f)


        # run main test blocks
        while (len(trials['load1']) +  len(trials['load2'])) > 0:
            block_num += 1
            assert block_num > 0 # block_num must be positive, else conflicts w/ 100 port code

            block, trials = self.retrieve_block_from_trials(trials)
            
            acc = []
            self.tracker.calibrate()

            self.send_synced_event(100+block_num)
            psychopy.core.wait(.01)
            for trial_num, trial in enumerate(block):
                if trial_num % self.n_trials_per_drift_correction == 0:
                    self.tracker.drift_correct()

                data = self.run_trial(trial, block_num, trial_num+1, realtime_eyetracking=realtime_eyetracking)
                if data:
                    self.send_data(data)
                    acc.append(data['ACC'])

            self.save_data_to_csv()

            # save out remaining trials in case of crash
            info_to_pickle = {'block_num': block_num, 'trials': trials, 'rejections': self.rejections}
            with open(remaining_info_file, 'wb+') as f:
                pickle.dump(info_to_pickle,f)
            
            self.display_block_feedback(acc)


        """
        Makeup Blocks
        """
        block_num = self.run_makeup_blocks(block_num=block_num, realtime_eyetracking=realtime_eyetracking, remaining_info_filename=remaining_info_file)  

        """
        Interruption Blocks
        """
        do_interrupt = self.display_text_screen(text = f'G blocks?', keyList=['y','n'])
        if do_interrupt == ['y']:

            for instruction in self.instruct_interruption_text:
                self.display_text_screen(text=instruction, keyList=self.stim['keyCodes'])

            # check for remaining interruption trials, generate if they don't exist
            remaining_info_file =f'TMP_{self.fileBase}_remaining_info_wINT.pickle'
            if os.path.exists(remaining_info_file) and self.extend_ok:
                with open(remaining_info_file, 'rb+') as f:
                    remaining_info = pickle.load(f)
                block_num = remaining_info['block_num']
                interrupt_trials = remaining_info['trials']
                self.rejections = remaining_info['rejections']
            else:
                interrupt_trials = self.make_trials(CB='load', n_trials=self.stim['nTrialsInterrupt'], phase='interrupt', interrupt=True)
                self.rejections = []
                info_to_pickle = {'block_num': block_num, 'trials': interrupt_trials, 'rejections': self.rejections}
                with open(remaining_info_file, 'wb+') as f:
                    pickle.dump(info_to_pickle,f)

            
            while (len(interrupt_trials['load1']) +  len(interrupt_trials['load2'])) > 0:
                block_num+=1
                assert block_num > 0
                
                block, interrupt_trials = self.retrieve_block_from_trials(interrupt_trials)
                
                acc = []
                self.tracker.calibrate()

                self.send_synced_event(100+block_num)
                psychopy.core.wait(.01)
                for trial_num, trial in enumerate(block):
                    if trial_num % self.n_trials_per_drift_correction == 0:
                        self.tracker.drift_correct()

                    data = self.run_trial_wInterruption(trial, block_num, trial_num+1, realtime_eyetracking=realtime_eyetracking)
                    if data:
                        self.send_data(data)
                        acc.append(data['ACC'])

                self.save_data_to_csv()
                # save out remaining trials in case of crash
                info_to_pickle = {'block_num': block_num, 'trials': interrupt_trials, 'rejections': self.rejections}
                with open(remaining_info_file, 'wb+') as f:
                    pickle.dump(info_to_pickle,f)

                self.display_block_feedback(acc)
            
            """
            Interruption Makeup Blocks
            """
            block_num += 1
            block_num = self.run_makeup_blocks(block_num=block_num, realtime_eyetracking=realtime_eyetracking, interrupt=True, remaining_info_filename=remaining_info_file)

        """
        End of Experiment
        """
        psychopy.core.wait(1)
        self.send_synced_event(self.portCodes['taskEnd'])
        self.display_text_screen(
            'The experiment is now over, please get your experimenter.',
            bg_color=[0, 0, 255], text_color=[255, 255, 255])
        
        self.tracker.transfer_edf()
        self.quit_experiment()
        
if __name__ == '__main__':
    exp = Acacia01(
        # BaseExperiment parameters
        experiment_name=exp_name,
        data_fields=data_fields,
        monitor_distance=distance_to_monitor,
        
        # Custom parameters go here
        track_eyes=True # False for debugging on mac
    )
    
    try:
        exp.run()
    except Exception as e:
        exp.kill_tracker()
        raise e