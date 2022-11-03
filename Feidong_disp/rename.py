## convert the EGFAnalysisTimeFreq output format to DisperPicker input format.

import os
import glob

os.system('mkdir -p TestData/group_image')
os.system('mkdir -p TestData/phase_image')
os.chdir('./disp')
for f in glob.glob('C.*'):
    sta = f.split('.')[1].split('_')
    new_name = sta[0] + '.' + sta[1] + '.dat'
    os.system('cp %s ../TestData/phase_image/%s'%(f, new_name))
for f in glob.glob('G.*'):
    sta = f.split('.')[1].split('_')
    new_name = sta[0] + '.' + sta[1] + '.dat'
    os.system('cp %s ../TestData/group_image/%s'%(f, new_name))
