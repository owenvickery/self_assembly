import os, sys
import numpy as np
from shutil import copyfile
import time
import re
import multiprocessing as mp
import subprocess 
from shutil import copyfile
from distutils.dir_util import copy_tree
from time import gmtime, strftime
from pathlib import Path
import argparse
import distutils.spawn
import math
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/database')

parser = argparse.ArgumentParser(description='Converts CG representation into an atomistic representation', prog='CG2AT', epilog='Enjoy the program and best of luck!\n', allow_abbrev=True)

parser.add_argument('-c', help='atomistic coordinates (Optional)',metavar='pdb/gro/tpr',type=str, required=True)
parser.add_argument('-l', help='residues to add number:residue',type=str, nargs='*')
parser.add_argument('-m', help='martini version to use',metavar='2.2 / 3.0',type=float, required=True)
parser.add_argument('-loc', help='output folder name, (default = SA_timestamp)',metavar='CG2AT',type=str)
parser.add_argument('-dssp', help='output folder name, (default = SA_timestamp)',metavar='CG2AT',type=str)
parser.add_argument('-gromacs', help='gromacs executable name (Optional)',metavar='gmx_avx',type=str)
parser.add_argument('-v', action="count", default=0, help="increase output verbosity (eg -vv, 3 levels) (Optional)")
parser.add_argument('-conc', help='salt concentration (nMol, default = 0.15)',metavar='0.15',default=0.15, type=float)
parser.add_argument('-run', help='run final SA step', action='store_true')
# parser.add_argument('-ions', help='Ions to be added (default = Na+ Cl-)',default=['Na+','Cl-'], nargs=2)
args = parser.parse_args()
options = vars(args)

def file_copy_and_check(file_in,file_out):
    if not os.path.exists(file_out):
        copyfile(file_in, file_out)

def folder_copy_and_check(folder_in,folder_out):
    if not os.path.exists(folder_out):
        copy_tree(folder_in, folder_out)

def gromacs(gro):
    cmd,output = gro[0], gro[1]
    if os.path.exists(output):
        pass
    else:
    #### if the flag gromacs is used every gromacs command will be printed to the terminal 
        if args.v >= 3:
            print('\nrunning gromacs: \n '+cmd+'\n')
        output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        err, out = output.communicate()
        exitcode = output.returncode
        out=out.decode("utf-8")
    #### all gromacs outputs will be saved into gromacs_outputs within the folder it is run
        with open('gromacs_outputs', 'a') as checks:
            checks.write(out)
    #### standard catch for failed gromacs commands
            if 'File input/output error:' in out:
                sys.exit('\n'+out)
            elif 'Error in user input:' in out:
                sys.exit('\n'+out)
            elif 'did not converge to Fmax ' in out:
                sys.exit('\n'+out)
            elif 'Segmentation fault' in out:
                sys.exit('\n'+out)
            elif 'but did not reach the requested Fmax' in out:
                sys.exit('\n'+out)
            elif 'Fatal error:' in out:
                sys.exit('\n'+out)
    if len(gro) == 4: 
        gro[3].put(gro[2])
        return gro[2]

def pdbatom(line):
### get information from pdb file
### atom number, atom name, residue name,chain, resid,  x, y, z, backbone (for fragment), connect(for fragment)
    try:
        return dict([('atom_number',int(line[7:11].replace(" ", ""))),('atom_name',str(line[12:16]).replace(" ", "")),('residue_name',str(line[17:21]).replace(" ", "")),\
            ('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54]))])
    except:
        sys.exit('\npdb line is wrong:\t'+line) 

def read_in_merged_pdbs(location):
    if os.path.exists(location):
        merge=[]
    #### opens pdb files and writes straight to merged_cg2at pdb
        with open(location, 'r') as pdb_input:
            for line in pdb_input.readlines():
                if line.startswith('ATOM'):
                    line_sep = pdbatom(line)
                    merge.append(line_sep)
        return merge
    else:
        sys.exit('cannot find minimised residue: \n'+ location) 

def fetch_gromacs():
    if args.gromacs != None:
        gmx=distutils.spawn.find_executable(args.gromacs)
    else:
        gmx=distutils.spawn.find_executable('gmx')
    if gmx==None or type(gmx) != str:
        for root, dirs, files in os.walk(os.environ.get("GMXBIN")):
            for file_name in files:
                if file_name.startswith('gmx') and file_name.islower() and '.' not in file_name:
                    gmx=distutils.spawn.find_executable(file_name)
                    if type(gmx) == str and gmx != None :
                        break
                    else:
                        gmx=None
            break
        if gmx == None:
            sys.exit('Cannot find gromacs installation')    
    return gmx

def fetch_dssp():
    if args.dssp != None:
        dssp=distutils.spawn.find_executable(args.dssp)
    else:
        dssp=distutils.spawn.find_executable('dssp')
    if dssp == None:
        sys.exit('Cannot find gromacs installation')    
    return dssp

def mkdir_directory(directory):
#### checks if folder exists, if not makes folder
    if not os.path.exists(directory):
        os.mkdir(directory)

def martinise(database_dir,final_dir, input_directory, dssp):
    if float(args.m) < 3:
        gromacs(['python '+database_dir+'/martinise_scripts/martinise.py -f '+input_directory+'AT_conversion.pdb -merge all -ff martini22 -elastic -ef 1000 -el 0.5 '+
                    '-eu 1 -ea 0 -ep 0 -x protein-cg.pdb -o protein-cg.top -dssp '+dssp, 'protein-cg.pdb']) 
        folder_copy_and_check(database_dir+'forcefields/martini_2-2.ff', final_dir+'martini_2-2.ff')
        return 'martini_2-2.ff'

    # else:
    #     os.system('python '+database_dir+'/martinise_scripts/martinise.py -f '+input_directory+'AT_conversion.pdb -merge all -ff martini22 -elastic -ef 1000 -el 0.5 '+
    #                 '-eu 1 -ea 0 -ep 0 -x protein-cg.pdb -o protein-cg.top -dssp '+dssp)

def create_directories(gmx):
    timestamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    if args.loc != None:
        working_dir_name = args.loc
    else:
        working_dir_name =  'SA_'+timestamp
    if not os.path.exists(args.c):
        sys.exit('\ncannot find input file: '+args.c)
    start_dir        = os.getcwd()+'/'  ### initial working directory
    working_dir      = os.getcwd()+'/'+working_dir_name+'/'
    input_directory  = os.getcwd()+'/'+working_dir_name+'/INPUT/' ### contains input run files
    merged_directory = os.getcwd()+'/'+working_dir_name+'/MERGED/'
    final_dir        = os.getcwd()+'/'+working_dir_name+'/FINAL/'  ### final directory for run files
    SA               = os.getcwd()+'/'+working_dir_name+'/SA/'
    database_dir     = os.path.dirname(os.path.realpath(__file__))+'/database/' ### contains script files
    mkdir_directory(working_dir)
    mkdir_directory(input_directory)
    mkdir_directory(merged_directory)
    mkdir_directory(final_dir)
    mkdir_directory(SA)
    os.chdir(input_directory)
    gromacs([gmx+' editconf -f '+start_dir+args.c+' -resnr 0 -o '+input_directory+'AT_conversion.pdb', input_directory+'AT_conversion.pdb'])

    return start_dir, working_dir, input_directory, merged_directory, final_dir, database_dir, SA 

def find_alignment(cg_protein):
    vector_coord=[]
    prev_resid=0
    vec_temp=[]
    reset=False
    COM=[]
    for bead in cg_protein:
        COM.append(np.array([bead['x'], bead['y'], bead['z']]))
        if bead['atom_name'] == 'BB':
            if bead['residue_name'] in ['VAL', 'LEU', 'ILE', 'PHE', 'ALA', 'MET', 'PRO', 'GLY', 'TRP','TYR','MET']:
                if bead['residue_id'] == prev_resid+1 or prev_resid == 0 or reset==True: 
                    prev_resid=bead['residue_id'] 
                    vec_temp.append(np.array([bead['x'], bead['y'], bead['z']]))
                    reset=False
            else:
                reset=True
                if len(vec_temp) > 4:
                    vector_coord.append(vec_temp)
                vec_temp=[]
    vectors=[]
    for vec in vector_coord:
        for coord_val, coord in enumerate(vec):
            if coord_val >= 4:
                vectors.append(convert_to_pos(coord-vec[coord_val-4]))
    normalised = np.sum(vectors, axis=0)/np.max(np.sum(vectors, axis=0), axis=0)
    return normalised, np.mean(COM, axis=0)

def convert_to_pos(coord):
    new=[]
    for number in coord:
        new.append(abs(number))
    return np.array(new)

def align_to_vector(v1, v2):
#### returns the rotation matrix to rotate v1 to v2
    v = np.cross(v1,v2)
    c = np.dot(v1,v2)
    s = np.linalg.norm(v)

    rotation=np.array([[   0, -v[2],  v[1]],
                       [v[2],     0, -v[0]],
                       [-v[1], v[0],    0], 
                       ])

    r = np.identity(3) - rotation + np.matmul(rotation,rotation) * ((1 - c)/(s**2))
    return r

def rotate_protein(cg_protein, rotation_matrix, COM):
    box = []
    for bead in cg_protein:
        coord_out = (np.array([bead['x'], bead['y'], bead['z']]) - COM).dot(rotation_matrix) + COM
        box.append(coord_out)
        bead['x']=coord_out[0]
        bead['y']=coord_out[1]
        bead['z']=coord_out[2]
    box_vec = return_box(np.array(box))
    return cg_protein, box_vec

def return_box(box):
    box_xy=[(np.max(box[:,i])-np.min(box[:,i]))+80 for i in range(2)]
    if abs(box_xy[0]-box_xy[1]) > 30:
        box_vec = box_xy
    else:
        box_vec = np.array([np.max(box_xy), np.max(box_xy)])
    box_vec = np.append(box_vec, (np.max(box[:,2])-np.min(box[:,2])))
    return box_vec

def write_pdb(merge, merged_directory, box_vec):
    
    pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"
    pdb_output=create_pdb(merged_directory+'protein_aligned.pdb', box_vec) 
    for line_val, line in enumerate(merge):
        pdb_output.write(pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
            line['x'],line['y'],line['z'],1,0))+'\n')

def create_pdb(file_name, box_vec_values):
    box_line="CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n"
    box_vec = box_line%(box_vec_values[0], box_vec_values[1], box_vec_values[2])
    pdb_output = open(file_name, 'w')
    pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n'+
                    box_vec+'MODEL        1\n')
    return pdb_output

def sort_lipids(res,res_num,loc, box):
    to_add={}
    box_area = box[0]*box[1]
    area_to_fit = box_area*(res_num/100)
    area_per_tail = {'MONO':30,'DI':60,'TRI':100,'TETRA':120}
    for tail in area_per_tail:
        if tail in loc:
            return str(int(area_to_fit/area_per_tail[tail]))

def add_molecules(merged_directory, box):
    res_to_add={}
    for res_val, res in enumerate(args.l): 
        res_info = args.l[res_val].split(':')
        res, res_num = res_info[0], int(res_info[1])
        loc = fragment_location(database_dir, res)
        res_to_add[res] = sort_lipids(res,res_num,loc, box)   
        if res_val == 0:
            print('Inserting '+res_to_add[res]+' '+res+' to the system') 
            gromacs([gmx+' insert-molecules -f '+merged_directory+'protein_centered.pdb -o '+merged_directory+'protein_'+res+'.pdb -ci '+
                    loc+' -radius 0.25 -nmol '+res_to_add[res], merged_directory+'protein_'+res+'.pdb'])
            new_file = 'protein_'+res
        else:
            print('Inserting '+res_to_add[res]+' '+res+' to the system') 
            gromacs([gmx+' insert-molecules -f '+merged_directory+new_file+'.pdb -o '+merged_directory+new_file+'_'+res+'.pdb -ci '+
                    loc+' -radius 0.25 -nmol '+res_to_add[res], merged_directory+new_file+'_'+res+'.pdb'])
            new_file+='_'+res
    file_copy_and_check(merged_directory+new_file+'.pdb', merged_directory+'protein_lipids.pdb')
    return res_to_add


def solvate():
    print('inserting water beads')
    if args.m < 3:
        loc = fragment_location(database_dir, 'W')
    else:
        loc = fragment_location(database_dir, 'WN')
    gromacs([gmx+' solvate -cp '+merged_directory+'protein_lipids.pdb -o '+merged_directory+'protein_lipids_water.pdb  -cs '+
          loc+' -radius 0.25 -box '+str(box_vec[0]/10)+' '+str(box_vec[1]/10)+' '+str((box_vec[2]+60)/10), merged_directory+'protein_lipids_water.pdb'])

def fetch_water_number():
    cg_solvated = read_in_merged_pdbs(merged_directory+'protein_lipids_water.pdb')
    sol=0
    for bead in cg_solvated:
        if bead['residue_name'] =='W':
            sol+=1
    return sol

def write_topol(final_dir, residues, forcefield):
    if not os.path.exists('topol_sol.top'):
        with open('topol_sol.top', 'w') as topol_write:       
            topol_write.write('; Include forcefield parameters\n')
            for file in os.listdir(final_dir+forcefield):
                if file.endswith('.itp'):
                    topol_write.write('#include \"'+final_dir+'/'+forcefield+'/'+file+'\"\n')
            topol_write.write('#define RUBBER_BANDS\n')
            topol_write.write('#include \"Protein_A.itp\"\n')
            topol_write.write('[ system ]\n; Name\nSomething clever....\n\n[ molecules ]\n; Compound        #mols\n')
            topol_write.write('{0:20}{1:^10}\n'.format('Protein_A', '1'))
            for residue_type in residues:
                topol_write.write('{0:20}{1:^10}\n'.format(residue_type, residues[residue_type]))
   
def write_min_mdp():
#### makes em.mdp file for each residue
    if not os.path.exists('em.mdp'):
        with open('em.mdp','w') as em:
            em.write('define = \n integrator = steep\nnsteps = 10000\nemtol = 750\nemstep = 0.001\ncutoff-scheme = Verlet\n')

def write_mdp(loc):
    if not os.path.exists(loc):
        with open(loc, 'w') as steered_md:
            steered_md.write('integrator = md\nnsteps = 10000000\nnstxout-compressed   = 10000\ndt = 0.02\ncontinuation   = no\nconstraint_algorithm = lincs\n')
            steered_md.write('constraints = none\nns_type = grid\nnstlist = 20\nrcoulomb = 1.1\nrvdw = 1.1\ncoulombtype  = Reaction_field\n')
            steered_md.write('tcoupl = v-rescale\ntc-grps = system\ntau_t = 1\nref_t = 350\nvdw-modifier = Potential-shift-verlet\nepsilon_rf = 0\n')
            steered_md.write('pcoupl = Berendsen\npcoupltype = semiisotropic\ntau_p = 12.0\nref_p = 1.0 1.0\ncompressibility = 3e-4 3e-4\n')
            steered_md.write('pbc = xyz\ngen_vel = yes\ngen_temp = 350\ngen_seed = -1\nrefcoord_scaling = com\ncutoff-scheme = Verlet\n')
            steered_md.write('lincs_order = 8\nlincs-iter = 2\nlincs_warnangle = 30\nepsilon_r = 15\nvdw_type = cutoff')

def minimise_protein_chain(chain, input):
#### grompps each protein chain
    gromacs([g_var.gmx+' grompp '
                '-f em.mdp '+
                '-p topol_sol.top '+
                '-c '+merged_directory+'protein_lipids_water.pdb'+
                '-o ions.tpr'
                '-maxwarn 1 ', 'ions.tpr'])
#### minimises chain
    os.chdir('min')
    gromacs([g_var.gmx+' mdrun -v -deffnm PROTEIN_'+input+str(chain)+' -c PROTEIN_'+input+str(chain)+'.pdb', 'PROTEIN_'+input+str(chain)+'.pdb'])
    os.chdir('..')

def genion(merged_directory):
    gromacs([gmx+' grompp'
            ' -f em.mdp'+
            ' -p topol_sol.top'+
            ' -c '+merged_directory+'protein_lipids_water.pdb'+
            ' -o ions.tpr'
            ' -maxwarn 1 ', 'ions.tpr'])

    gromacs([gmx+' genion -s ions.tpr -o complete_system.pdb -neutral -p topol_sol.top -pname NA+ -nname CL- -conc '+str(args.conc)+' << EOF\nW\nWN\nEOF\n', 'complete_system.pdb'])

def minimise_system(merged_directory):
    os.chdir(merged_directory)
    mkdir_directory('min')
    gromacs([gmx+' grompp'
            ' -f em.mdp'+
            ' -p topol_sol.top'+
            ' -c '+merged_directory+'complete_system.pdb'+
            ' -o min/system_min.tpr'
            ' -maxwarn 1 ', 'min/system_min.tpr'])    
    os.chdir('min')
    gromacs([gmx+' mdrun -v -deffnm system_min -c system_min.pdb', 'system_min.pdb'])
    os.chdir('..')  

def move_to_final(merged_directory, final_dir):
    file_copy_and_check(merged_directory+'topol_sol.top', final_dir+'topol.top')
    gromacs([gmx+' make_ndx -f '+merged_directory+'min/system_min.pdb -o '+merged_directory+'index.ndx << EOF\nq\nEOF\n', merged_directory+'index.ndx'])
    file_copy_and_check(merged_directory+'index.ndx', final_dir+'index.ndx')
    gromacs([gmx+' editconf -f '+merged_directory+'min/system_min.pdb -c -o '+final_dir+'final_system.pdb -n '+merged_directory+'index.ndx << EOF\nProtein\nSystem\nq\nEOF\n', final_dir+'final_system.pdb'])
    file_copy_and_check(merged_directory+'Protein_A.itp', final_dir+'Protein_A.itp')

def grompp_final(final_dir, SA):
    os.chdir(final_dir)
    
    gromacs([gmx+' grompp'
            ' -f self_assemble.mdp'+
            ' -p topol.top'+
            ' -c final_system.pdb'+
            ' -o ../SA/SA.tpr'
            ' -maxwarn 1 ', '../SA/SA.tpr'])     
    if args.run:
        print('Running 100 ns of self assembly')
        os.chdir('../SA')
        gromacs([gmx+' mdrun -v -pin on -nt 8 -deffnm SA -c SA.pdb', 'SA.pdb'])

def fragment_location(database_dir, residue):  
#### runs through dirctories looking for the atomistic fragments returns the correct location
    for root, dirs, files in os.walk(database_dir+'residues'):
        for directory in dirs:
            for root, dirs, files in os.walk(database_dir+'residues/'+directory):
                if residue+'.gro' in files:
                    # print(database_dir+'residues/'+directory+'/'+residue+'.gro')
                    return database_dir+'residues/'+directory+'/'+residue+'.gro' 
                break    
        break
    sys.exit('cannot find fragment: '+residue)

def check_percentages():
    lip_per=0
    for lipid in args.l:
        lip_per+= float(lipid.split(':')[1])
    if lip_per != 100:
        sys.exit('lipid percentages does not equal 100%')


check_percentages()
gmx = fetch_gromacs()
dssp = fetch_dssp()
start_dir, working_dir, input_directory, merged_directory, final_dir, database_dir, SA  = create_directories(gmx)
os.chdir(merged_directory)
forcefield = martinise(database_dir,final_dir, input_directory, dssp)


#### align protein roughly to the z-axis
cg_protein = read_in_merged_pdbs('protein-cg.pdb')
align, COM= find_alignment(cg_protein)
rotation_matrix = align_to_vector(align, [0,0,1])
cg_protein_aligned, box_vec  = rotate_protein(cg_protein, rotation_matrix, COM)
write_pdb(cg_protein_aligned, merged_directory, box_vec)
gromacs([gmx+' editconf -f '+merged_directory+'protein_aligned.pdb -c -o '+merged_directory+'protein_centered.pdb', merged_directory+'protein_centered.pdb'])


#### add lipid molecules
if args.l != None:
    res_to_add = add_molecules(merged_directory, box_vec)
else:
    file_copy_and_check(merged_directory+'protein_centered.pdb', merged_directory+'protein_lipids.pdb')
### solvate
solvate()
res_to_add['w'] = fetch_water_number()

### write setup files
write_topol(final_dir,res_to_add, forcefield)
write_min_mdp()
write_mdp(final_dir+'self_assemble.mdp')

### add_ions
genion(merged_directory)

minimise_system(merged_directory)
move_to_final(merged_directory, final_dir)
grompp_final(final_dir, SA)



