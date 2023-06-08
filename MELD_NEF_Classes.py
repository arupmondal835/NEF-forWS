#!/usr/bin/env python
# coding: utf-8

import numpy as np 
import pandas as pd
import tempfile
import argparse
import os
from tabulate import tabulate
import pickle



header = '''data_nef_MELD_{}
   
  save_nef_nmr_meta_data
      _nef_nmr_meta_data.sf_category          nef_nmr_meta_data
      _nef_nmr_meta_data.sf_framecode         nef_nmr_meta_data
      _nef_nmr_meta_data.format_name          nmr_exchange_format
      _nef_nmr_meta_data.format_version       1.1
      _nef_nmr_meta_data.program_name         MELD-NEF 
      _nef_nmr_meta_data.program_version      1.0
      _nef_nmr_meta_data.creation_date        2020-03-31
      _nef_nmr_meta_data.uuid                 MELD-NEF_2020-03-31
      _nef_nmr_meta_data.coordinate_file_name . 
      
      loop_
          _nef_program_script.program_name
          _nef_program_script.script_name
          
          MELD MELD-NEF
       stop_
  save_
  
  
  '''

def parse_args():                              #in line argument parser with help 
    '''
    Parse arguments for the program from command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-nef', type=str, help='NEF input file')
    parser.add_argument('-out', type=str, help='NEF output file name',default='MELD.nef')
    parser.add_argument('-directory', type=str, help='Where to place the files',default='.')
    return(parser.parse_args())
 




class NEF_system:
    '''This will be an object keeping track of all the blocks that compose the NEF file.
    We will have a list telling us block orders and a dictionary to map to each block'''
    def __init__(self, filename, directory):
        '''Start the system keeping track of the original file and processing all blocks inside'''
        self.block_order = []
        self.block_content = {}
        self.block_types = {}
        self.active = []
        self.original_file = filename
        self.directory = directory
        self.get_blocks(filename)
        
    def add_block(self,block):
        '''Given a block object, we insert it into the global NEF system'''
        if block.name:
            block_naming = "_".join([block.type,block.name])
        else:
            block_naming = block.type
            
        self.block_order.append(block_naming)
        self.block_content[block_naming] = block
        try: 
            self.block_types[block.type].append(block)
        except:
            self.block_types[block.type] = [block]
    
    def get_blocks(self,file_name):
        '''open a nef file and get the blocks from the file.
        Input: a string as file name.
        Return: a list of all the blocks(also a list). 
        In one block list, lines are elements in the list.
        '''
        all_blocks = []
        block = []
        logic = False
        with open(file_name, 'r') as fin:
            for line in fin:
                line_strip = '{}\n'.format(line.strip())
                if 'save_\n' in line_strip:
                    logic = False
                    self.add_block(NEF_block(block))
                    block = []
                if 'save_nef_' in line_strip:
                    logic = True
                if logic:
                    block.append(line_strip)
                    
    def write(self,file_name='MELD_test1.nef'):
        '''write the whole NEF pipeline. The first block will change from whatever we were given to know that
        it has been hadnled by our MELD pipeline. For all other blocks, we will only print the ones that are active.
        If no block has been set as active, it should print all'''
        
        txt = ""
        txt += header.format(file_name)
        if len(self.active) < 1:
            self.active = self.block_order[1:]

        
        for block_name in self.active:
            #update the output object
            self.block_content[block_name].write()
            #and set it to our txt 
            txt += self.block_content[block_name].output
            txt += '   save_\n'
        
        with open(file_name,'w') as fo:
            fo.write(txt)
            
    def sequence(self):
        '''NEF.block_content['molecular_system'].loop_type_data['_nef_sequence']
        This has a pandas kind of structure. index/chain/sequence/residue/linking/residue/cis
        '''
        data = self.block_content['molecular_system'].loop_type_data['_nef_sequence']
        chains = data.chain_code.unique()
        self.chains = chains
        self.sequence = {}
        self.sequence_names = []
        for chain in chains:
            residues = data.residue_name[data.chain_code == chain].values
            # Add caps. This will fail if not a protein. TODO: make more general
            allowed_residues_in_meld = ["ASH","GLH","HIE","HID","HIP","LYN","ACE","OHE","NME","NHE","HIS","ALA","CYS","ASP","GLU","PHE","GLY","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
            
            #residues[0] = '{}'.format(residues[0])
            #residues[-1] = '{}'.format(residues[-1])
            protein_residues=[]
            for AA3 in residues:
                if AA3 not in allowed_residues_in_meld:
                    print('Warning: non protein residues')
                else:
                    protein_residues.append(AA3)
            with open('{}/sequence_{}.dat'.format(self.directory,chain),'w') as fo:
                fo.write(" ".join(protein_residues))
                fo.write("\n")
            self.sequence['sequence_{}'.format(chain)] = " ".join(protein_residues)
            self.sequence_names.append('sequence_{}'.format(chain))
            


class NEF_block:
    '''Each block is refered to by the save_nef_[type]_[name] keyworkd. Inside the block
    each _nef_[type] has several attributes. Then some of those attributes have a corresponding loop_
    structure, where first the headers of that attribute name are given, and then the information as a table.
    We want to create the same structure for all blocks'''
    def __init__(self,block):
        '''Each block should have the original data before processing. 
        Attributes:
        _raw: contains the original data
        type: the type of block we are dealing with
        name: if present will help distinguish between different blocks with same type
        attributes: the list of keywords this particular block has. These will effectively be
            the keys to dictionaries containing such data
        '''
        self._raw = block
        self.attributes = {}
        self.attributes_ordered = []
        self.loop_type = []
        self.loop_type_headers = {}
        self.loop_type_data = {}
        self.type = None
        self.name = None
        self.header = None
        
        
        self.process_type(block[0].rstrip())
        self.block_from_raw()
        
    def process_type(self,header):
        '''Gets the type of data in this block as well as the name for this type if any
        Will return None otherwise. Block header has this format: save_nef_[type]_[name]'''
        if ('list' in header):
            tmp = header.split('list_')
            self.type = "_".join(tmp[0].split('_')[2:]) + 'list'
            self.name = tmp[1]
        elif ('spectrum' in header):
            tmp = header.split('spectrum_')
            self.type = "_".join(tmp[0].split('_')[2:]) + 'spectrum'
            self.name = tmp[1]
        else:
            self.type = "_".join(header.split('_')[2:])
            self.name = None
        self.header = header
        print(self.name,self.type,self.header)
        
    def process_attributes(self,data):
        '''First part of the block. Tells us about the attributes in this block'''

        for line in data.split('\n'):
            if self.type in line and not 'save_' in line:
                result = line.split('{}.'.format(self.type))
                result = result[1].split()
                self.attributes_ordered.append(result[0])
                self.attributes[result[0]] = result[1].rstrip()
            else:
                # Empty lines
                pass
        print(self.type,self.attributes_ordered)
    
    def process_loop(self,data):
        '''Each loop has headers and then values. We will likely create a pandas structure'''
        headers = []
        txt = ""
        loop_type = None
        header = None
        #print('processing loop')
        for line in data.split('\n'):
            if '_nef_' in line:
                [loop_type, header] = line.split('.')
                headers.append(header)
            elif line.strip():
                #this is the table
                if ('stop_' or 'save_') in line:
                    pass
                else:
                    txt += '{}\n'.format(line)
            else:
                pass
            
        #print(headers)
        #print(txt)
            
        #write a temporary file and read into a pandas data frame
        with open('tmp.txt','w') as fo:
            fo.write("{}\n".format(" ".join(headers)))
            fo.write(txt)
        self.loop_type.append(loop_type)
        self.loop_type_headers[loop_type] = headers
        self.loop_type_data[loop_type] = pd.read_csv('tmp.txt',sep=" ",skipinitialspace=True)
        
    def block_from_raw(self):
        '''Give the initial block structure as specified in the nef file. 
        All blocks start with the header, then the properties in that header, then the tables of data
        for each specific property (keyword loop_)
        Will generate one big text and then split into loop regions.
        The first one has the attributes for this block
        The next ones (if any) have the table information -- all blocks have at least one loop_ keyword
        '''
        block = "".join(self._raw)
        data = block.split('loop_')
        
        for i,loop in enumerate(data):
            #print(i)
            #print(loop)
            if i < 1:
                self.process_attributes(loop)
            else:
                self.process_loop(loop)

    def write(self):
        '''Each block object should be able to write itself. First the attributes, then each loop'''
        self.output = '{}\n\n'.format(self.header)
        for attribute in self.attributes_ordered:
            self.output += '      _nef_{}.{} {}'.format(self.type, attribute, self.attributes[attribute])
            self.output += '\n'
        self.output += '\n'
        
        #Now process the loops
        for loop_type in self.loop_type:
            self.output += '      loop_\n'
            #Handle all the headers
            for loop_header in self.loop_type_headers[loop_type]:
                self.output += '         {}.{}\n'.format(loop_type,loop_header)
            self.output += '\n'
            #Use pandas to write out a formatted table of the data in the loop
            #Note that we have the formatters option to tells us how to write! formatterslist
            #If each loop type has a formatting option we just have to define all of them somewhere at the 
            #beginning and then the output will always be exact!
            #For now i'm just going to do standard formatting
            self.loop_type_data[loop_type]
            #prepend 7 white spaces to the data frame
            self.output += '\n'
            df = self.loop_type_data[loop_type] 
            df[df.columns[0]] = '      ' + df[df.columns[0]].astype(str)
            #self.output += df.to_string(header=False, justify='match-parent',index=False)
            self.output += tabulate(df, showindex=False, tablefmt='plain')
            #Finish the loop
            self.output += '\n    stop_\n\n'
        
        #Finish the block
        self.output += "\n\n"


'''
Example:
#This is an example of how to use the block object
all_blocks = get_blocks('CCPN_Commented_Example.nef')
a = NEF_block(all_blocks[0])
print('loop type: ',a.loop_type)
print('loop type headerse',a.loop_type_headers['         _nef_program_script'])
print('loop type data', a.loop_type_data['         _nef_program_script'])




#This is an example of how to use the NEF_system class, which reads and processes all blocks
NEF = NEF_system('CCPN_Commented_Example.nef')

#we need the block order to write out file in same order
print(NEF.block_order)
#This will group all types together in case we want all the noesy for example, and will return a list of the 
#block objects
print(NEF.block_types.keys())
print(NEF.block_content['molecular_system'].loop_type_headers)
print(NEF.block_types['distance_restraint_list'])

print(NEF.block_content['molecular_system'].loop_type)
print(NEF.block_content['molecular_system'].loop_type_headers)
print(NEF.block_content['molecular_system'].loop_type_data['         _nef_sequence'].to_string(header=False))
NEF.block_content['molecular_system'].loop_type_data['         _nef_sequence']

NEF.write()
for i in NEF.block_content.keys():
    print(i,NEF.block_content[i].attributes_ordered)
'''





def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def main():
    args = parse_args()
    #Work in a temporary directory
    if args.directory == '.':
        args.directory = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        NEF = NEF_system(args.nef,args.directory)
        print(NEF.block_content['molecular_system'].loop_type_data['_nef_sequence'])
        NEF.sequence()
        #Now save the object in a way we can use for MELD
        save_object(NEF, '{}/{}.pkl'.format(args.directory,args.out))
        #And save a NEF file
        NEF.write(file_name='{}/{}'.format(args.directory,args.out))


if __name__ == '__main__': #Python way to execute main()
    main()
