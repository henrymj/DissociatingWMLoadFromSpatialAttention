import pandas as pd
import numpy as np

def clean_behavior(behavior, noArt_idx, isub):
    if isub >= 6: # exclude practice
        behavior = behavior.query('Block > 0').reset_index(drop=True) 
    if isub == 8:  # sub 8 (013) did interruption trials, but they weren't captured in the EEG
        behavior = behavior[behavior.args_phase=='test'].reset_index(drop=True)
    
    behavior = behavior.loc[noArt_idx, :].reset_index(drop=True)
    if isub < 3: # first two subs had bad column names
        behavior = behavior.rename(columns={'targetCloud_nDots': 'targetCloud_color',
                                            'targetCloud_color': 'targetCloud_nDots',
                                            'otherCloud_bins': 'otherCloud_color',
                                            'otherCloud_nDots_1': 'otherCloud_bins_1',
                                            'otherCloud_nDots_2': 'otherCloud_bins_2',
                                            'otherCloud_nDots_3': 'otherCloud_bins_3',
                                            'otherCloud_color': 'otherCloud_nDots'})
    behavior = behavior[behavior.columns.drop(list(behavior.filter(regex='idx_')))]
    
    behavior['SetSize'] = np.where(behavior.otherCloud_nDots==0, 1, 2)
    if isub >= 6:  # psychopy subjects
        behavior['port_codes'] = behavior['memPort']
        behavior['targetBinList'] = behavior.targetCloud_bins.map(lambda x: x.strip('][').split(' '))
        behavior['targetWidth'] = behavior['targetBinList'].map(len)
        behavior['targetBin'] = None
        for row in behavior.itertuples():
            behavior.loc[row.Index, 'targetBin'] = row.targetBinList[0] if len(row.targetBinList)==1 else row.targetBinList[1]
        behavior['targetBin'] = behavior['targetBin'].astype(int)
        
        behavior['otherBinList'] = behavior.otherCloud_bins.map(lambda x: x.strip('][').split(' '))
        behavior['otherWidth'] = behavior['otherBinList'].map(len)
        behavior['otherBin'] = None
        for row in behavior.itertuples():
            if len(row.otherBinList)==1:
                val = row.otherBinList[0]
                val = np.nan if val=='' else int(val)
            elif len(row.otherBinList)==3:
                val = int(row.otherBinList[1])
            else:
                print('ERROR: ', behavior.otherCloud_bins)
                
            behavior.loc[row.Index, 'otherBin'] = val
            
        # switch to 1 index to match first subjects
        behavior['targetBin'] = behavior.targetBin.values + 1
        behavior['otherBin'] = behavior.otherBin.values + 1
    else:
        behavior['ACC'] = behavior['ChoiceAcc']
        behavior['targetWidth'] = np.where(behavior.targetCloud_bins_2.isnull(), 1, 3)
        behavior['targetBin'] = np.where(behavior.targetCloud_bins_2.isnull(), behavior.targetCloud_bins_1, behavior.targetCloud_bins_2)
        
        behavior['otherWidth'] = np.where(behavior.otherCloud_bins_2.isnull(), 1, 3)
        behavior['otherBin'] = np.where(behavior.otherCloud_bins_2.isnull(), behavior.otherCloud_bins_1, behavior.otherCloud_bins_2)
        
    behavior['total_nDots'] = np.where(behavior['SetSize']==2, behavior['targetCloud_nDots']+behavior['otherCloud_nDots'], behavior['targetCloud_nDots'])
    
    behavior['sameWidth'] = behavior['targetWidth'] == behavior['otherWidth']
    
    return behavior

def setup_hyperplane_conds(behavior, separate_no_overlap=False):
    behavior['hyperplane_conditions'] = None
    behavior.loc[(behavior.SetSize==1) & (behavior.targetWidth==1), 'hyperplane_conditions'] = 0 # Narrow
    behavior.loc[(behavior.SetSize==1) & (behavior.targetWidth==3), 'hyperplane_conditions'] = 1 # broad

    # set size 2 conditions
    behavior.loc[(behavior.SetSize==2), 'hyperplane_conditions'] = 6  # set size 2
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==1) & (behavior.otherWidth==1) & (behavior.targetBin==behavior.otherBin), 'hyperplane_conditions'] = 2  # set size 2, complete overlap, narrow
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==3) & (behavior.otherWidth==3) & (behavior.targetBin==behavior.otherBin), 'hyperplane_conditions'] = 3  # set size 2, complete overlap, broad
    
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth!=behavior.otherWidth) & (((behavior.targetBin-behavior.otherBin)==0) | ((np.abs(behavior.targetBin-behavior.otherBin))%6==1)) , 'hyperplane_conditions'] = 4  # set size 2, superset overlap
    
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==3) & (behavior.otherWidth==3) & ((np.abs(behavior.targetBin-behavior.otherBin)%6)==1), 'hyperplane_conditions'] = 5  # set size 2, partial overlap, 1 shift
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==3) & (behavior.otherWidth==3) & (((np.abs(behavior.targetBin-behavior.otherBin))==2) | ((np.abs(behavior.targetBin-behavior.otherBin))==6)), 'hyperplane_conditions'] = 5  # set size 2, partial overlap, 2 shift
    
    if separate_no_overlap:
        # break apart no-overlap conditions - 6=mixed, 7=both narrow, 8=both broad
        behavior.loc[ ((behavior.hyperplane_conditions==6) & (behavior.targetWidth==1) & (behavior.otherWidth==1)), 'hyperplane_conditions'] = 7  # both narrow
        behavior.loc[ ((behavior.hyperplane_conditions==6) & (behavior.targetWidth==3) & (behavior.otherWidth==3)), 'hyperplane_conditions'] = 8  # both broad

    return behavior


def setup_attended_area(behavior):
    behavior['attended_area'] = None

    # SET SIZE 1
    behavior.loc[(behavior.SetSize==1) & (behavior.targetWidth==1), 'attended_area'] = 1 # Narrow
    behavior.loc[(behavior.SetSize==1) & (behavior.targetWidth==3), 'attended_area'] = 3 # broad

    # SET SIZE 2

    # defaults, to be overwritten below - if both 1, area = 2; if both 3, area = 6; if mixed, area = 4
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==1) & (behavior.otherWidth==1), 'attended_area'] = 2  # set size 2, complete overlap, narrow
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==3) & (behavior.otherWidth==3), 'attended_area'] = 6  # set size 2, complete overlap, broad
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth!=behavior.otherWidth), 'attended_area'] = 4  # set size 2, complete overlap, broad

    # complete overlap
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==1) & (behavior.otherWidth==1) & (behavior.targetBin==behavior.otherBin), 'attended_area'] = 1
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==3) & (behavior.otherWidth==3) & (behavior.targetBin==behavior.otherBin), 'attended_area'] = 3

    # Superset - area must be 3
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth!=behavior.otherWidth) & (((behavior.targetBin-behavior.otherBin)==0) | ((np.abs(behavior.targetBin-behavior.otherBin))%6==1)) , 'attended_area'] = 3

    # partial overlap - 4 or 5
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==3) & (behavior.otherWidth==3) & ((np.abs(behavior.targetBin-behavior.otherBin)%6)==1), 'attended_area'] = 4  # 1 shift
    behavior.loc[(behavior.SetSize==2) & (behavior.targetWidth==3) & (behavior.otherWidth==3) & (((np.abs(behavior.targetBin-behavior.otherBin))==2) | ((np.abs(behavior.targetBin-behavior.otherBin))==6)), 'attended_area'] = 5  # 2 shift
    
    return behavior