#!/usr/bin/python2

import sys
import os
import commands
from commands import getstatusoutput
from commands import getoutput
import datetime
import argparse

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

class customAction(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        setattr(args, self.dest, values.split(','))

def get_comma_separated_args(self, arg_line):
    return arg_line.split(',')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def lxbatchSubmitJob (run, path, cfg, outdir, queue, job_dir, dryrun):
    jobname = job_dir+'/H4Reco_'+queue+'_'+run+'.sh'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('export X509_USER_PROXY=/cmshome/marzocchi/myProxy/x509up_u28021 \n\n')
#    f.write ('git clone --recursive https://github.com/simonepigazzini/H4Analysis.git \n')
#    f.write ('cd H4Analysis/ \n')
#    f.write ('cp /cmshome/marzocchi/H4Analysis/cfg/*.cfg cfg/ \n')
    f.write ('cd ../.. \n') 
#    f.write ('source scripts/setup.sh \n')
#    f.write ('make -j 2 \n')
#    f.write ('cp '+path+cfg+' job.cfg \n\n')
#    f.write ('cp '+path+'/ntuples/Template*.root ./ntuples/ \n\n')
    f.write ('bin/H4Reco '+path+cfg+' '+run+'\n\n')
    f.write ('cd /cmshome/marzocchi/CMSSW_8_0_8_patch1/src/ \n')
    f.write ('eval `scramv1 runtime -sh` \n')
    f.write ('cd /cmshome/marzocchi/TestBeam_latest/H4Analysis/ntuples/ \n')
#    f.write ('ls \n')
    if "store/" in outdir:
        f.write ('xrdcp -f *'+run+'.root '+outdir+'\n')
    else:
        f.write ('cp ntuples/*'+run+'.root '+outdir+'\n')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    if not dryrun:
        getstatusoutput ('cd '+job_dir+'; bsub -q ' + queue + ' ' + '-u badder.marzocchi@cern.ch ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def herculesSubmitJob (run, path, cfg, outdir, queue, job_dir, dryrun):
    jobname = job_dir+'/H4Reco_'+queue+'_'+run+'.sh'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('mkdir -p /gwpool/users/$USER/pool/'+run+'/ \n')
    f.write ('cd /gwpool/users/$USER/pool/'+run+'/ \n')
    f.write ('wget https://github.com/simonepigazzini/H4Analysis/archive/master.zip \n')
    f.write ('unzip master.zip \n')
    f.write ('cd H4Analysis-master/ \n\n')
    f.write ('cp '+path+'/ntuples/Template*.root ./ntuples/ \n')
    f.write ('cp '+path+cfg+' job.cfg \n')
    f.write ('source scripts/setup.sh \n')
    f.write ('make -j \n\n')
    f.write ('bin/H4Reco job.cfg '+run+'\n\n')
    f.write ('cp ntuples/*'+run+'.root '+outdir+'\n')
    f.write ('cd /gwpool/users/$USER/pool/ \n')
    f.write ('rm -r '+run+' \n')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    if not dryrun:
        getstatusoutput ('cd '+job_dir+'; qsub -q ' + queue + ' ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'submit H4Reco step to lsf')
    parser.add_argument('-r', '--runs' , action=customAction, help='run to be processed, either list or file')
    parser.add_argument('-q', '--queue' , default = '1nh', help='batch queue (1nh)')
    parser.add_argument('-s', '--storage' , default = '/store/group/dpg_ecal/alca_ecalcalib/ECALTB_H4_Fall2015/', help='storage path')
    parser.add_argument('-v', '--version' , default = 'v1', help='production version')
    parser.add_argument('-c', '--cfg' , default = '../cfg/H4DAQ_base.cfg', help='production version')
    parser.add_argument('--dryrun' , action="store_true", default=False, help='do not submit the jobs, just create them')
    parser.add_argument('--batch' , default='lxbatch', help='batch system to use')
    
    args = parser.parse_args ()

    ## check ntuple version
    stageOutDir = args.storage+'ntuples_'+args.version+'/'
    
    if args.batch == 'lxbatch':
        if getoutput('gfal-ls root://eoscms/'+stageOutDir) == "":
            print "ntuples version "+args.version+" directory on eos already exist! no jobs created."
            #exit(0)
        getstatusoutput('cmsMkdir '+stageOutDir)    
    else:
        getstatusoutput('mkdir -p '+stageOutDir)    
    
    ## job setup
    local_path = getoutput('pwd')
    date = datetime.datetime.now().strftime("%d-%m-%Y")
    job_dir = local_path+"/"+date+"_ntuples_"+args.version
    getstatusoutput('mkdir -p '+job_dir)
    if local_path.endswith('scripts'):
        local_path = local_path[:-len('scripts')]

    if args.cfg.startswith('../'):
        args.cfg = args.cfg[len('../'):]
    
    if len(args.runs) == 1 and os.path.isfile(args.runs[0]):
        runs_file = open(args.runs[0], 'r')
        args.runs  = []
        if runs_file:
            for run in runs_file:
                args.runs.append(run.rstrip())

    ## create jobs
    print 'submitting jobs to queue', args.queue
    for run in args.runs:
        print 'submitting run: ', run
        if args.batch == 'lxbatch':
            lxbatchSubmitJob(run, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun) 
        if args.batch == 'hercules':
            herculesSubmitJob(run, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun) 
