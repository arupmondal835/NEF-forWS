meld_NMR_script = '''#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, leader
import meld.system.montecarlo as mc
from meld import system
from meld.system import patchers
from meld import comm, vault
from meld import parse
from meld import remd
from meld.system import param_sampling
from openmm import unit as u
import mdtraj as md
import glob as glob
import sys
from Bio import SeqIO


N_REPLICAS = 10
N_STEPS = 1000
BLOCK_SIZE = 50

def gen_state(s, index):
    state = s.get_state_template()
    state.alpha = index / (N_REPLICAS - 1.0)
    return state


def get_dist_restraints(filename, s, scaler, ramp, seq):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]
            dist = float(cols[4])/10.0
            #
            # note: manually over riding distances
            #
            #dist = 0.45
        
            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=dist*u.nanometer, r4=(dist+0.2)*u.nanometer, 
                                                 k=350*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)
    return dists


def get_torsion_restraints(filename, s, scaler,ramp,seq):
    torsion_rests = []
    rotamer_group = []
    lines = open(filename).read().splitlines()
    #for line in open(filename,'r'):
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            torsion_rests.append(s.restraints.create_restraint_group(rotamer_group, 1))
            rotamer_group = []
        else:
            (res1, at1, res2, at2, res3, at3, res4, at4, rotamer_min, rotamer_max) = line.split()
            rotamer_max = float(rotamer_max)
            rotamer_min = float(rotamer_min)
            rotamer_avg = (rotamer_max+rotamer_min)/2.
            rotamer_sd = abs((rotamer_max - rotamer_min)/2.)
            rotamer_rest = s.restraints.create_restraint('torsion', scaler, ramp,
                                                 phi=rotamer_avg*u.degree, delta_phi=rotamer_sd*u.degree, k=2.5*u.kilojoule_per_mole/u.degree **2,
                                                 atom1=s.index.atom(int(res1)-1,at1, expected_resname=seq[int(res1)-1][-3:]),
                                                 atom2=s.index.atom(int(res2)-1,at2, expected_resname=seq[int(res2)-1][-3:]),
                                                 atom3=s.index.atom(int(res3)-1,at3, expected_resname=seq[int(res3)-1][-3:]),
                                                 atom4=s.index.atom(int(res4)-1,at4, expected_resname=seq[int(res4)-1][-3:]))
                                        
            rotamer_group.append(rotamer_rest)
    return torsion_rests

def get_seq(temp):
    f=open('sequence.dat','w')
    with open(temp, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            #print('>' + record.id)
            f.write('{}'.format(record.seq))


def setup_system():
    # load the sequence

    sequence = parse.get_sequence_from_AA3(filename='sequence_A.dat')
    n_res = len(sequence.split())

    # build the system
    p = system.subsystem.SubSystemFromSequence(sequence)
    b = system.builder.SystemBuilder()

    #rdc_patcher = patchers.RdcAlignmentPatcher(n_tensors=1)
    #s = b.build_system([p], patchers=[rdc_patcher])
    s = b.build_system([p])


    s.temperature_scaler = system.temperature.GeometricTemperatureScaler(0, 0.4, 300.*u.kelvin, 500.*u.kelvin)

    ramp = s.restraints.create_scaler('nonlinear_ramp', start_time=1, end_time=200,
                                      start_weight=1e-3, end_weight=1, factor=4.0)
    seq = sequence.split()
    for i in range(len(seq)):
        if seq[i][-3:] =='HIE': seq[i]='HIS'
    print(seq)


    #
    # Setup Scaler
    #
    scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=0.8, factor=4.0)
    scaler_short = s.restraints.create_scaler('nonlinear', alpha_min=0.8, alpha_max=1.0, factor=4.0)
    scaler_rot=s.restraints.create_scaler('constant')
    # Torsion Restraints
    #

    #
    # Distance Restraints

    if "--PS" in sys.argv:
        print('using parameter sampling')
        for noe in glob.glob('NOE_*.dat'):
            print('loading {} ...'.format(noe))
            NOESY = get_dist_restraints(noe,s,scaler,ramp,seq)
            prior_c13 = param_sampling.ScaledExponentialDiscretePrior(u0=1.0, temperature_scaler=s.temperature_scaler, scaler=scaler)
            sampler_c13 = param_sampling.DiscreteSampler(int(1), int(1.00 * len(NOESY)), 1)
            param_c13 = s.param_sampler.add_discrete_parameter("param_c13", int(1.00 * len(NOESY)), prior_c13, sampler_c13)
            s.restraints.add_selectively_active_collection(NOESY,param_c13)
    

        if len(glob.glob('NOE_*.dat')) < 1:
            print('WARNING: no highContactOrder NOE data loaded')

        for noe in glob.glob('local_NOE_*.dat'):
            print('loading {} ...'.format(noe))
            NOESY2 = get_dist_restraints(noe,s,scaler_short,ramp,seq)
            prior_n15 = param_sampling.ScaledExponentialDiscretePrior(u0=1.0, temperature_scaler=s.temperature_scaler, scaler=scaler_short)
            sampler_n15 = param_sampling.DiscreteSampler(int(1), int(1.00 * len(NOESY2)), 1)
            param_n15 = s.param_sampler.add_discrete_parameter("param_n15", int(1.00 * len(NOESY2)), prior_n15, sampler_n15)
            s.restraints.add_selectively_active_collection(NOESY2, param_n15)        

        if len(glob.glob('local_NOE_*.dat')) < 1:
            print('WARNING: no shortContactOrder NOE data loaded')

    if "--PS" not in sys.argv:
        print('NOT using parameter sampling')
        for noe in glob.glob('NOE_*.dat'):
            print('loading {} ...'.format(noe))
            NOESY = get_dist_restraints(noe,s,scaler,ramp,seq)
            s.restraints.add_selectively_active_collection(NOESY, int( len(NOESY)*1.00 ) )

        if len(glob.glob('NOE_*.dat')) < 1:
            print('WARNING: no highContactOrder NOE data loaded')

        for noe in glob.glob('local_NOE_*.dat'):
            print('loading {} ...'.format(noe))
            NOESY2 = get_dist_restraints(noe,s,scaler_short,ramp,seq)
            s.restraints.add_selectively_active_collection(NOESY2, int( len(NOESY2)*1.00 ) )

        # Torsion Restraints
    #
    for rotamer in glob.glob('rotamers_*.dat'):
        print('loading {} ...'.format(rotamer))
        TALOS = get_torsion_restraints(rotamer, s, scaler_rot,ramp, seq)
        s.restraints.add_selectively_active_collection(TALOS, int( len(TALOS) * 1.00) )

    if len(glob.glob('rotamers_*.dat')) < 1:
        print('WARNING: no rotamer data found')



    # setup mcmc at startup
    movers = []
    n_atoms = s.n_atoms
    for i in range(0, n_res):
        n = s.index.atom(i, 'N', expected_resname=seq[i][-3:])
        ca = s.index.atom(i, 'CA', expected_resname=seq[i][-3:])
        c = s.index.atom(i, 'C', expected_resname=seq[i][-3:])
 
        atom_indxs = list(system.indexing.AtomIndex(j) for j in range(ca,n_atoms))
        mover = mc.DoubleTorsionMover(index1a=n, index1b=ca, atom_indices1=list(system.indexing.AtomIndex(i) for i in range(ca, n_atoms)),
                                      index2a=ca, index2b=c, atom_indices2=list(system.indexing.AtomIndex(j) for j in range(c, n_atoms)))

        movers.append((mover, 1))

    sched = mc.MonteCarloScheduler(movers, n_res * 60)

    # create the options
    options = system.options.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.use_big_timestep = True
    options.cutoff = 1.8*u.nanometers
    options.remove_com = False
    options.use_amap = False
    options.amap_beta_bias = 1.0
    options.timesteps = 14286
    options.minimize_steps = 20000
    options.min_mc = sched
    options.param_mcmc_steps=200

    # create a store
    store = vault.DataStore(s.get_state_template(),N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1, min_acc_prob=0.02)

    remd_runner = remd.leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS,
                                                            ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS, timeout=60000)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms


setup_system()
'''

meld_gpu_job = '''#!/bin/bash
#SBATCH --job-name={0}     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1      
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-gpu=1
#SBATCH --mem-per-cpu=2000mb            
#SBATCH --partition=gpu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2000mb          # Memory per processor
#SBATCH --time=100:00:00              # Time limit hrs:min:sec
#SBATCH --output=meld_{0}_%j.log     # Standard output and error log
#SBATCH --constraint=2080ti

ml cuda/11.0.207 gcc/9.3.0 mkl/2020.0.166 openmpi/4.0.4 openmm/7.6.0 python-core/3.8

export LD_PRELOAD=/opt/pmix/3.1.2/lib64/libpmix.so
export LD_PRELOAD=/home/arup.mondal/apps/meld_0.5/5.0/lib/libMeldPlugin.so:$LD_PRELOAD

export PATH=/apps/cuda/11.0.207/gcc/9.3.0/openmpi/4.0.4/amber/20/bin:$PATH
export PATH=/home/arup.mondal/apps/meld_0.5/5.0/bin:$PATH

export PYTHONPATH=/home/arup.mondal/apps/meld_0.5/5.0/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=/apps/cuda/11.0.207/gcc/9.3.0/openmpi/4.0.4/openmm/7.6.0/lib64/python3.8/site-packages:$PYTHONPATH

export OPENMM_PLUGIN_DIR=$HPC_OPENMM_LIB/plugins
export OPENMM_CUDA_COMPILER=$HPC_CUDA_BIN/nvcc
export LD_LIBRARY_PATH=/home/arup.mondal/apps/meld_0.5/5.0/lib:$LD_LIBRARY_PATH

for i in 1 2 3 4 5 6
do
if [ -e Logs/remd_000.log ]; then             #If there is a remd.log we are conitnuing a killed simulation
    prepare_restart --prepare-run  #so we need to prepare_restart
      fi
launch_remd_multiplex --debug
done
'''
