# -*- coding: utf-8 -*- 

###########################################################################
## Python code generated with wxFormBuilder (version Jun 17 2015)
## http://www.wxformbuilder.org/
###########################################################################

import wx
import wx.xrc
import os
import os.path as op
import shutil as sh
import json
from subprocess import call

class installer ( wx.Frame ):
    
    def __init__( self ):
        super().__init__ (parent=None, id = wx.ID_ANY, title = u"SR toolkit v3.6 installer", pos = wx.DefaultPosition, size = wx.Size( 450,250 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
        panel = wx.Panel(self)

        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
        
        bSizer1 = wx.BoxSizer( wx.VERTICAL )
        
        self.m_staticText1 = wx.StaticText( panel, wx.ID_ANY, u"Install SR toolkit to...", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText1.Wrap( -1 )
        bSizer1.Add( self.m_staticText1, 0, wx.ALL, 5 )
        
        self.m_dirPicker1 = wx.DirPickerCtrl( panel, wx.ID_ANY, wx.EmptyString, u"Select a folder", wx.DefaultPosition, wx.DefaultSize, wx.DIRP_DEFAULT_STYLE )
        bSizer1.Add( self.m_dirPicker1, 0, wx.ALL, 5 )
        
        self.m_staticText2 = wx.StaticText( panel, wx.ID_ANY, u"ImageJ app path where GDSC SMLM 1 is installed...", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText2.Wrap( -1 )
        bSizer1.Add( self.m_staticText2, 0, wx.ALL, 5 )
        
        self.m_filePicker2 = wx.FilePickerCtrl( panel, wx.ID_ANY, wx.EmptyString, u"Select a file", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
        bSizer1.Add( self.m_filePicker2, 0, wx.ALL, 5 )
        
        self.m_staticText3 = wx.StaticText( panel, wx.ID_ANY, u"ImageJ app path where GDSC SMLM 2 is installed (leave blank if not needed)...", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText3.Wrap( -1 )
        bSizer1.Add( self.m_staticText3, 0, wx.ALL, 5 )
        
        self.m_filePicker3 = wx.FilePickerCtrl( panel, wx.ID_ANY, wx.EmptyString, u"Select a file", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
        bSizer1.Add( self.m_filePicker3, 0, wx.ALL, 5 )
        
        self.m_button1 = wx.Button( panel, wx.ID_ANY, u"Install", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer1.Add( self.m_button1, 0, wx.ALL, 5 )
        
        self.m_button1.Bind(wx.EVT_BUTTON, self.on_press)
        
        panel.SetSizer( bSizer1 )
        panel.Layout()
        
        panel.Centre( wx.BOTH )
        self.Show()

    def on_press(self, event):
        install_dir = self.m_dirPicker1.GetPath()
        imageJ1_path = self.m_filePicker2.GetPath()
        imageJ2_path = self.m_filePicker3.GetPath()
    
        if not install_dir:
            wx.MessageBox("Please enter the install directory for SR toolkit.", "Error", wx.OK | wx.ICON_ERROR)
            return
        if not imageJ1_path:
            wx.MessageBox("Please enter the path to ImageJ.", "Error", wx.OK | wx.ICON_ERROR)
            return
        if not imageJ2_path:
            dlg = wx.MessageDialog(None, "GDSC SMLM2 path left blank, SR toolkit will run with GDSC SMLM1 only. Do you want to proceed?", "GDSC SMLM2", wx.YES_NO | wx.ICON_QUESTION)
            result = dlg.ShowModal()

            if result == wx.ID_YES:
                imageJ2_path = imageJ1_path
            else:
                return
        exitcode = self.install_dependency()
        if exitcode == 1:
            wx.MessageBox('Installation of dependencies failed, please contact Eric (feedbackeric@gmail.com).', 'Failed', wx.OK | wx.ICON_INFORMATION)

        exitcode = self.install_toolkit(install_dir, imageJ1_path, imageJ2_path)
        if exitcode == 0:
            wx.MessageBox('Installation completed.', 'Completed', wx.OK | wx.ICON_INFORMATION)
        else:
            wx.MessageBox('Installation of SR toolkit failed, please contact Eric (feedbackeric@gmail.com).\n Error code: {}'.format(exitcode), 'Failed', wx.OK | wx.ICON_INFORMATION)
        exit()

    def install_dependency(self):
        try:
            call('pip install -r requirements.txt')
        except:
            return 1
        return 0

    def install_toolkit(self, setup_path, imageJ_GDSCSMLM1_path, imageJ_GDSCSMLM2_path):
        try:
            setup_path = op.join(setup_path, 'SR_toolkit')
            if op.isdir(setup_path):
                sh.rmtree(setup_path)
            sh.copytree('SR_toolkit', setup_path)
            if os.getenv('PYTHONPATH') is None or setup_path not in os.getenv('PYTHONPATH'):
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
        except Exception as e:
            return e
        else:
            return 0
        
if __name__ == '__main__':
    app = wx.App()
    frame = installer()
    app.MainLoop()
    