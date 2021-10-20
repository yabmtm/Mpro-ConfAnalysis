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

crystal_runs = [[14365, 770], [14365, 771], [14367, 147], [14367, 148], [14367, 149], [14367, 150], [14367, 171], [14367, 172], [14367, 173], [14367, 174], [14367, 195], [14724, 204], [14724, 205], [14724, 206], [14724, 211], [14724, 212], [14724, 213], [14724, 214], [14724, 292], [14724, 298], [14724, 965], [14724, 973], [14724, 975], [14724, 983], [14823, 2194], [14723, 2527], [14723, 770], [14723, 771], [14724, 147], [14724, 148], [14724, 149], [14724, 150], [14724, 1586], [14724, 1587], [14824, 1703], [14824, 1704], [14824, 1705], [14824, 1706], [14824, 171], [14824, 172], [14824, 1725], [14824, 1726], [14824, 1727], [14367, 196], [14367, 213], [14369, 771], [14723, 2526], [14724, 172], [14724, 203], [14823, 2495], [14824, 1702], [14824, 1728], [14824, 1748], [14371, 147], [14371, 148], [14371, 149], [14371, 150], [14371, 171], [14371, 172], [14371, 173], [14371, 174], [14371, 195], [14371, 196], [14371, 197], [14371, 198], [14371, 203], [14371, 204], [14371, 205], [14371, 206], [14371, 211], [14371, 212], [14371, 213], [14371, 214], [14371, 292], [14371, 298], [14371, 965], [14371, 973], [14371, 975], [14371, 983], [14723, 2194], [14723, 2495], [14723, 2496], [14723, 2497], [14723, 2525], [14724, 1725], [14724, 1726], [14724, 1727], [14724, 1728], [14724, 1729], [14724, 173], [14724, 1730], [14724, 174], [14724, 1743], [14724, 1744], [14724, 1745], [14724, 1746], [14724, 1747], [14724, 1748], [14724, 195], [14724, 196], [14724, 197], [14724, 198], [14824, 195], [14824, 196], [14824, 197], [14824, 198], [14824, 203], [14824, 204], [14824, 205], [14824, 206], [14824, 211], [14824, 212], [14824, 213], [14824, 214], [14824, 292], [14824, 298], [14824, 965], [14824, 973], [14824, 975], [14824, 983], [14367, 197], [14367, 198], [14367, 203], [14367, 204], [14367, 205], [14367, 206], [14367, 211], [14367, 212], [14367, 214], [14367, 292], [14367, 298], [14367, 965], [14367, 973], [14367, 975], [14367, 983], [14369, 770], [14823, 2496], [14823, 2497], [14823, 2525], [14823, 2526], [14823, 2527], [14823, 770], [14823, 771], [14824, 147], [14824, 148], [14824, 149], [14824, 150], [14824, 1586], [14824, 1587], [14824, 1588], [14824, 1589], [14824, 1590], [14824, 1591], [14824, 1701], [14724, 1588], [14724, 1589], [14724, 1590], [14724, 1591], [14724, 1701], [14724, 1702], [14724, 1703], [14724, 1704], [14724, 1705], [14724, 1706], [14724, 171], [14824, 1729], [14824, 173], [14824, 1730], [14824, 174], [14824, 1743], [14824, 1744], [14824, 1745], [14824, 1746], [14824, 1747]]

url_prefix = 'https://fah-public-data-covid19-absolute-free-energy.s3.us-east-2.amazonaws.com'
clones = range(1) # just CLONE0
for [project,run] in crystal_runs:
    for clone in clones:
        if not os.path.exists(f'data/P{project}_R{run}_C{clone}'):
            os.makedirs(f'data/P{project}_R{run}_C{clone}')
        for file in ['npt.gro', 'topol.top', 'prod.mdp']:
            download(f'{url_prefix}/setup_files/p{project}/RUN0/{file}',
              f'data/P{project}_R{run}_C{clone}/{file}')
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
            
            try:
              for cmd in [write_mdp_cmd, make_index_cmd, make_xtctop_cmd, make_xtcgro_cmd,
                make_xtcndx_cmd, make_xtctpr_cmd, pbc_correct_cmd]:
                  subprocess.check_output(cmd, stderr=subprocess.STDOUT,shell=True).decode().split('\n')
            except Exception as e:
                print(f'Error processing P{project}_R{run}_C{clone}: {e}')
                continue
      
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
