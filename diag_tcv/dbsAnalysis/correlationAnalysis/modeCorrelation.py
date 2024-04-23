#%%

def isModeCorrelation(shot):
    modeCorrelation={}
    modeCorrelation['78549'] = True
    modeCorrelation['79797'] = False
    if str(shot) not in modeCorrelation:
        print(' --- modeCorrelation not filled: putting False as default --- ')
        return False
    return modeCorrelation[str(shot)]

