import mdtraj as md
import itertools
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import os, urllib, subprocess, glob
from tqdm import tqdm

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)
        
### pull data down from AWS and post-process trajectories

# download xtc and (gro/top)

url_prefix = 'https://fah-public-data-covid19-absolute-free-energy.s3.us-east-2.amazonaws.com'
project = 14823
runs = range(1) # just RUN0
clones = range(1) # just CLONE0
for run in runs:
    for clone in clones:
        if not os.path.exists(f'data/P{project}_R{run}_C{clone}'):
            os.makedirs(f'data/P{project}_R{run}_C{clone}')
        download(f'{url_prefix}/setup_files/p{project}/RUN0/npt.gro',
          f'data/P{project}_R{run}_C{clone}/npt.gro')
        download(f'{url_prefix}/setup_files/p{project}/RUN0/topol.top',
          f'data/P{project}_R{run}_C{clone}/topol.top')
        gen = 0
        while True:
            try:
                print(f'\nProcessing P{project}_R{run}_C{clone}_G{gen}')
                download(f'{url_prefix}/PROJ{project}/RUN{run}/CLONE{clone}/results{gen}/traj_comp.xtc',
                  f'data/P{project}_R{run}_C{clone}/traj_comp.xtc')
            except Exception as e:
                print(e)
                break
            
            path = f'data/P{project}_R{run}_C{clone}'
            
            write_mdp_cmd = f'echo "integrator          = steep" > {path}/xtc.mdp'
            make_index_cmd = f'echo "5|2\nq\n" | gmx make_ndx -f {path}/npt.gro -o {path}/index.ndx'
            make_xtctop_cmd = f'head -n-3 {path}/topol.top > {path}/xtc.top'
            make_xtcgro_cmd = f'echo "24\n" | gmx editconf -f {path}/npt.gro -n {path}/index.ndx -o {path}/xtc.gro'
            make_xtcndx_cmd = f'echo "3|2\nq\n" | gmx make_ndx -f {path}/xtc.gro -o {path}/xtc.ndx'
            make_xtctpr_cmd = f'gmx grompp -f {path}/xtc.mdp -c {path}/xtc.gro -p {path}/xtc.top -o {path}/xtc.tpr'
            pbc_correct_cmd = f'echo "3\n14\n" | gmx trjconv -f {path}/traj_comp.xtc -s {path}/xtc.tpr -n {path}/xtc.ndx -pbc mol -center -o {path}/traj_{str(gen).zfill(4)}.xtc'
            
            for cmd in [write_mdp_cmd, make_index_cmd, make_xtctop_cmd, make_xtcgro_cmd,
              make_xtcndx_cmd, make_xtctpr_cmd, pbc_correct_cmd]:
                subprocess.check_output(cmd, stderr=subprocess.STDOUT,shell=True).decode().split('\n')
            
            
            traj = md.load(f'{path}/traj_{str(gen).zfill(4)}.xtc',top = f'{path}/xtc.gro')
            PHE140_indices = [a.index for a in traj.topology.atoms if a.residue.index in [141] and a.name in ['CG','CD1','CD2','CE1','CE2','CZ']]
            HIS163_indices = [a.index for a in traj.topology.atoms if a.residue.index in [164]and a.name in ['CG','ND1','CD2','CE1','NE2']]
            traj_PHE140_indices = traj.atom_slice(PHE140_indices)
            traj_HIS163_indices = traj.atom_slice(HIS163_indices)
            coords_PHE140_com = md.compute_center_of_mass(traj_PHE140_indices)
            coords_HIS163_com = md.compute_center_of_mass(traj_HIS163_indices)
            hacked_traj = traj

            ## creating hacked traj 0 and 1
            hacked_traj.xyz[:,0,:] = coords_PHE140_com # PHE140 trajectory
            hacked_traj.xyz[:,1,:] = coords_HIS163_com # HIS163 trajectory


            ## computing the distance between the center of mass of the PHE140 and HIS163 ring
            PHE140_HIS163_distances = md.compute_distances(hacked_traj, [[0,1]])[:,0]
            np.save(f'{path}/PHE140_HIS163_distnces_G{str(gen).zfill(4)}', PHE140_HIS163_distances)

            gen += 1
        
            file_list = ['traj_comp.xtc','xtc.mdp','xtc.top','xtc.ndx','xtc.gro','xtc.tpr','index.ndx']
            for file in glob.glob(f'{path}/*'):
                if any(substring in file for substring in file_list):
                    os.remove(file)
                    
           #gen += 1
