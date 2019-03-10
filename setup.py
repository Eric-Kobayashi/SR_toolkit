import os
import os.path as op
import shutil as sh
import json
import site

if __name__ == '__main__':
    # Here is where the script folder will be
    setup_path = r"C:\Users\Public\Documents\Python Scripts\SR_toolkit"
    
    # Put imageJ path here
    imageJ_GDSCSMLM1_path = r"C:\Users\yz520\Downloads\fiji-win64-20170530\Fiji.app\ImageJ-win64.exe"
    imageJ_GDSCSMLM2_path = r"C:\Users\yz520\Desktop\fiji-win64-GDSCSMLM2\Fiji.app\ImageJ-win64.exe"
    
    # ======= setup code, no need to change =======
    if op.isdir(setup_path):
        sh.rmtree(setup_path)
    script_path = os.path.dirname(os.path.abspath(__file__))
    sh.copytree(op.join(script_path, 'SR_toolkit'), setup_path)
    site.addsitedir(setup_path)
    
    config_file = {}
    config_file['json_file'] = op.join(setup_path, 'to_ImageJ.json').replace('\\','//')
    config_file['imageJ_GDSCSMLM1_path'] = imageJ_GDSCSMLM1_path
    config_file['imageJ_GDSCSMLM2_path'] = imageJ_GDSCSMLM2_path
    config_file['GDSCSMLM1_script_path'] = op.join(setup_path, 'Peakfit_GDSC_SMLM.py')
    config_file['GDSCSMLM2_script_path'] = op.join(setup_path, 'Peakfit_GDSC_SMLM2.py')
    config_file['Rendering_script_path1'] = op.join(setup_path, 'Rendering_SR.py')
    config_file['Rendering_script_path2'] = op.join(setup_path, 'Rendering_SR2.py')
    with open(op.join(setup_path, 'config.txt'), 'w') as f:
        json.dump(config_file, f)
    
    