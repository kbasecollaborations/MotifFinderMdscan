import sys                                                                                                      
import os                                                                                                       
import json                            
import re
import numpy as np
from Bio import motifs
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from io import StringIO
#from installed_clients.DataFileUtilClient import DataFileUtil
import shutil
import subprocess

class MdscanUtil:
  def __init__(self):
      pass

  def build_mdscan_motif_command(self, inputFilePath, motiflen, prb):
      cmd1 = 'cp -r /kb/module/deps/kb_mdscan/MDscan.linux /kb/module/work/tmp/mdscan_out'
      subprocess.call('mkdir /kb/module/work/tmp/mdscan_out', shell=True)
      subprocess.call(cmd1, shell=True)
      os.chdir("/kb/module/work/tmp/mdscan_out")
      #command = '/kb/module/work/tmp/mdscan_out/MDscan.linux -f /kb/module/work/tmp/SeqSet.fa -w '+ str(motiflen)
      command = '/kb/module/work/tmp/mdscan_out/MDscan.linux -i /kb/module/work/tmp/SeqSet.fa -w '+ str(motiflen) + ' -o /kb/module/work/tmp/mdscan_out/SeqSet.out'
      return command

  def run_mdscan_command(self, command):
      print (command)
      os.system(command)

  def write_obj_ref(self, path, obj_ref):
      file = open(path+"/mdscan_obj.txt","w")
      file.write(obj_ref)
      file.close() 

  def parse_mdscan_output(self, path):
      outputFileList = []
 
      seqflag=False
      motifList={}
      motifDict={}
      locList=[]
      alphabet=['A','C','G','T']
      print(path)
  
      
      motifList['Condition']='temp'
      motifList['SequenceSet_ref']='123'

      background={}
      background['A']=0.0
      background['C']=0.0
      background['G']=0.0
      background['T']=0.0

      
      

      
      pwmList=[]
      pwmDict={}
      pfmDict={}
      rowList = []
      rowDict={}
      sequence=''
      idflag=0
      for filename in os.listdir(path):
          outputFileList.append(path + '/' + filename)
          if(filename=="SeqSet.out"):
             outputFilePath=path+'/'+filename
             #print(outputFilePath)
             samplerFile = open(outputFilePath,'r')
             motifSet=[]
             for line in samplerFile:
                 line=line.rstrip()
                 
                 if(line.startswith("Motif")):
                   
                   if(len(motifDict)!=0):
                      #print("+++")
                      motifSet.append(motifDict)
                      #print(motifDict)
                      motifDict={}
                      #print("+++")
                   motifDict['Motif_Locations'] = []
                   motifDict['PWM'] = []
                   motifDict['PFM'] = []
                   locList=[]
                   seqflag=True
                   seq=line.split(";")
                   consensus=(seq[3]).split(" ")[2]
                   pwmDict['A']=[]
                   pwmDict['C']=[]
                   pwmDict['G']=[]
                   pwmDict['T']=[]

                   pfmDict['A']=[]
                   pfmDict['C']=[]
                   pfmDict['G']=[]
                   pfmDict['T']=[]
                   motifDict['PWM']=pwmDict
                   motifDict['PFM']=pfmDict
                   motifDict['Iupac_sequence']=consensus
                   #print(consensus)
                   locDict={}
                 if(seqflag):
                     
                     if(line.startswith(">")):
                        out=line.split("\t")
                        idflag=True
                        #print(out)
                        rec=(out[3]).split(" ")

                        seqid=(out[0]).replace(">", "")
                        #print(seqid)
                        orientation=rec[0]
                        if(orientation == "f"):
                           seq_start=int(rec[1])
                           seq_end=int(rec[1])+8
                           orientation='+'
                        else:
                           seq_start=int(rec[1])-8
                           seq_end=int(rec[1])
                           orientation='-'

                        '''sequence=(out[3]).replace(";", "").replace("\"", "")'''
                        #locDict={}
                        locDict['sequence_id']=seqid;
                        locDict['start']=seq_start;
                        locDict['end']=seq_end;
                        #locDict['sequence']=sequence;
                        locDict['orientation']=orientation;
                        #motifDict['Motif_Locations'].append(locDict)
                     else:
                        if(idflag):
                           sequence=line
                           #print(line)
                           locDict['sequence']=sequence;
                           #print(locDict)
                           motifDict['Motif_Locations'].append(locDict)
                           locDict={}
                           idflag=False
             else:
                 #print("+++")
                 motifSet.append(motifDict)
                 #print(motifDict)
                 motifDict={}
                 #print("+++")
                 
             #print(motifSet)       
    
      motifList['Motifs']=motifSet
      motifList['Background']=background
      motifList['Alphabet']=alphabet   
      return motifList

  def UploadFromMdscan(self, callback_url, params):
          """
          :param params: instance of type "UploadmfmdInParams" -> structure:
             parameter "path" of String, parameter "ws_name" of String,
             parameter "obj_name" of String
          :returns: instance of type "UploadOutput" -> structure: parameter
             "obj_ref" of String
          """
          # ctx is the context object
          # return variables are: output
          #BEGIN UploadFromMdscan
          print('Extracting motifs')
          '''motifList = self.parse_sampler_output(params['path'])
          print(motifList)
       
          MSO = {}
          MSO=motifList
        
          dfu = DataFileUtil(callback_url)
          save_objects_params = {}
          save_objects_params['id'] = dfu.ws_name_to_id(params['ws_name'])
          save_objects_params['objects'] = [{'type': 'KBaseGeneRegulation.MotifSet' , 'data' : MSO , 'name' : params['obj_name']}]
          
          info = dfu.save_objects(save_objects_params)[0]
          print('SAVED OBJECT')
          print(info)
          motif_set_ref = "%s/%s/%s" % (info[6], info[0], info[4])
          print(motif_set_ref)
          output = {'obj_ref' : motif_set_ref}
          print(output)'''

        
          #exit("test")
          #END UploadFromMdscan

          # At some point might do deeper type checking...
          if not isinstance(output, dict):
              raise ValueError('Method UploadFrommfmd return value ' +
                             'output is not type dict as required.')


          # return the results
          return [output]

MDU=MdscanUtil()
output=MDU.parse_mdscan_output("/home/manish/Desktop/reorganization/reorg/MotifFinderMdscan/test_local/workdir/tmp/mdscan_out")
print(output)


