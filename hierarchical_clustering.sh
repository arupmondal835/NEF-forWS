#! /bin/bash
#ASK FOR 1 node with 8 CPUs
#SBATCH --job-name=cluster      # Job name
#SBATCH --ntasks=1                # Number of MPI ranks
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --mem-per-cpu=14gb          # Memory per processor
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=cluster%j.log     # Standard output and error log
#SBATCH  --account=alberto.perezant
#SBATCH  --qos=alberto.perezant
pwd; hostname; date
 





for r in 1.5
do

cat<<EOF>get_res.py
#! /usr/bin/env python

f = open('ss.dat','r')

l = f.readlines()[0].rstrip()

ok = False
res = []
for i,s in enumerate(l):
    if ok and s not in 'HE':
        ok = False
        end = i
        res.append("{}-{}".format(ini,end))

    if (s in 'HE'):
        if not ok:
            ok = True
            ini = i + 1
        if (ok):
            pass
all=",".join(res)
print(":{}".format(all))
EOF

res=`python ./get_res.py`
echo $res

# Select residues on which to do the clustering
# If sse < 50%, cluster on all CA,CB
# If sse > 50%, cluster on sse CA,CB
# seqlen
seqlen=$(awk '{print $1,gsub(/./,"")}' ss.dat | awk '{print $2}')
H=$(awk '{print $1,gsub(/H/,"")}' ss.dat | awk '{print $2}')
B=$(awk '{print $1,gsub(/E/,"")}' ss.dat | awk '{print $2}')
totsse=$(echo "$H + $B" | bc -l)
fracss=$(echo "($B + $H) / $seqlen" | bc -l)

if (( $(echo "$fracss < 0.50" | bc -l) )); then
    RMSDMASK=":1-$seqlen@CA,CB"
else
    RMSDMASK="$res@CA,CB"
fi

#generate topology
cat<<EOF>get_top.py
#! /usr/bin/env python
import pickle
x = pickle.load(open('Data/system.dat','rb'))
f = open('topol.top', 'w')
f.write(x.top_string)
EOF
python get_top.py

mkdir Cpptraj_linkage_sieve_eps_$r
cd Cpptraj_linkage_sieve_eps_$r

cpptraj ../topol.top<<EOF>&ptraj_err
trajin  ../trajectory.00.dcd 1 1 
go
EOF

Frames=`grep reading ptraj_err|awk '{print $NF}'|sed 's/)//'`
#Let's use half of the trajectory for equilibration... read the last half of the trajectory
start=$(echo "scale=0;$Frames/2"|bc -l)

# cluster
cat<<EFO>cpptraj.in
trajin ../trajectory.00.dcd $start $Frames
trajin ../trajectory.01.dcd $start $Frames
#trajin ../trajectory.02.dcd $start $Frames
#trajin ../trajectory.03.dcd $start $Frames
#trajin ../trajectory.04.dcd $start $Frames

rms PredSSE first $RMSDMASK out trajrmsd.dat
cluster hieragglo epsilon $r linkage rms $RMSDMASK sieve 10 summary summary singlerepout representative repout unique repfmt pdb clusterout clusttraj avgout avg avgfmt pdb out frame_vs_cluster.txt
go
EFO

for j in 0 1 2 3 4 5 6 7 8 9 
do
    cat<<EOF>>cpptraj.in
   clear trajin
   trajin  clusttraj.c$j
   trajout clusttraj.c$j.pdb model
   average average.c$j.pdb pdb
   reference unique.c$j.pdb [uniq$j]
   rms centr$j ref [uniq$j] @CA 
   atomicfluct out back.$j.apf @C,CA,N,O byres
   go
EOF
done


cpptraj -p ../topol.top -i cpptraj.in
cd ..
done


