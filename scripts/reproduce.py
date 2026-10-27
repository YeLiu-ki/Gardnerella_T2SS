"""Run from any directory using the repository's installed Python environment."""
from pathlib import Path
import subprocess, sys, os
root=Path(__file__).resolve().parents[1]
env=os.environ.copy()
env['MPLBACKEND']='Agg'
env.setdefault('MPLCONFIGDIR',str(root/'generated/.matplotlib'))
for script in ['export_source_data.py','plot_figure_E.py','plot_figure_F.py']:
    subprocess.run([sys.executable,str(root/'scripts'/script)],cwd=root,env=env,check=True)
