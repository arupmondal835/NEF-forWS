#! /bin/env python
import os

#pdb = ['1pqx','2jr2','2juw','2k2e','2kcu','2kko','2ko1','2k07','2kpu','2kw5','2kzn','2loy','2luz','2png','6nbn']
pdb = ['2loy']
#pdb = ['2jr2']
for i in pdb:
    os.mkdir(i)
    os.chdir(i)
    os.system('python /orange/alberto.perezant/arup.mondal/NEF/test_for_web_server/MELD_NEF_Classes.py -nef /orange/alberto.perezant/arup.mondal/NEF/test_for_web_server/nef_files/{0}/{0}.nef -dir . -out {0}_MELD_it1.nef'.format(i))
    os.system('python /orange/alberto.perezant/arup.mondal/NEF/test_for_web_server/processNMR.py -nef {0}_MELD_it1.nef.pkl -name {0}'.format(i))
    os.system('python MELD_NMR_setup.py --PS')
    os.system('bash /orange/alberto.perezant/arup.mondal/NEF/test_for_web_server/chain.sh')
    os.chdir('..')


