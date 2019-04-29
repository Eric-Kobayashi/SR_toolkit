import os
import os.path as op
import shutil as sh
import json
from subprocess import call

if __name__ == '__main__':
    # Here is where the script folder will be
    setup_path = r"C:\Users\Public\Documents\Python Scripts\SR_toolkit"
    
    # Put imageJ path here
    imageJ_GDSCSMLM1_path = r"C:\Users\Eric\fiji-win64\Fiji.app\ImageJ-win64.exe"
    imageJ_GDSCSMLM2_path = r"C:\Users\Public\Documents\fiji-win64-GDSCSMLM2\Fiji.app\ImageJ-win64.exe"
    
    # ======= setup code, no need to change =======
    if op.isdir(setup_path) and op.basename(setup_path) == 'SR_toolkit':
        sh.rmtree(setup_path)
    if op.basename(setup_path) != 'SR_toolkit':
        setup_path = op.join(setup_path, 'SR_toolkit')
    script_path = os.path.dirname(os.path.abspath(__file__))
    sh.copytree(op.join(script_path, 'SR_toolkit'), setup_path)
    if setup_path not in os.getenv('PYTHONPATH'):
        call('setx PYTHONPATH "{}"'.format(setup_path))
    
    config_file = {}
    config_file['imageJ_GDSCSMLM1_path'] = imageJ_GDSCSMLM1_path
    config_file['imageJ_GDSCSMLM2_path'] = imageJ_GDSCSMLM2_path
    config_file['GDSCSMLM1_script_path'] = op.join(setup_path, 'Peakfit_GDSC_SMLM.py')
    config_file['GDSCSMLM2_script_path'] = op.join(setup_path, 'Peakfit_GDSC_SMLM2.py')
    config_file['Rendering_script_path1'] = op.join(setup_path, 'Rendering_SR.py')
    config_file['Rendering_script_path2'] = op.join(setup_path, 'Rendering_SR2.py')
    with open(op.join(setup_path, 'config.txt'), 'w') as f:
        json.dump(config_file, f)
    
    