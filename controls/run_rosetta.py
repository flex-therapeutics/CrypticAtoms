"""Two apo-only, sequence-preserving Rosetta controls, scored separately.
Run one sample per task from the same prepared apo input. No holo or ligand input.
"""
from pathlib import Path
import os,sys,json,time,argparse,hashlib,subprocess,concurrent.futures
ROOT=Path(__file__).resolve().parent;REVIEW=ROOT.parent;REPO=REVIEW
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
def prepare():
 import numpy as np
 import biotite.structure.io.pdb as pdb
 sys.path.insert(0,str(REPO))
 from benchmark.metrics import load_protein_chain,parse_contig_residues
 targets=json.loads((REPO/'benchmark/targets.json').read_text());manifest=[]
 for t in targets:
  a=load_protein_chain(REPO/t['apo']['path'],t['apo']['chain'])
  # Use the same observed apo contig as the benchmark generation inputs.
  ids=parse_contig_residues(t['apo']['contig'])[t['apo']['chain']] if t['apo'].get('contig') else set(a.res_id.tolist())
  a=a[np.isin(a.res_id,list(ids)) & (a.element!='H')]
  assert np.sum(a.atom_name=='CA')>0,t['name']
  key=t['apo']['pdb_code']+'_'+t['apo']['chain']
  if t['name']=='HIV_RT':key+='_original' # distinct archived apo preparation
  out=ROOT/'inputs'/f'{key}.pdb';out.parent.mkdir(exist_ok=True)
  pp=pdb.PDBFile();pp.set_structure(a);tmp=out.with_suffix(f'.{os.getpid()}.tmp');pp.write(tmp);os.replace(tmp,out)
  manifest.append({'target':t['name'],'input_key':key,'path':str(out),'n_res':int(np.sum(a.atom_name=='CA'))})
 (ROOT/'inputs.json').write_text(json.dumps(manifest,indent=2))
 return manifest
def worker(args):
 import fcntl
 lockpath=ROOT/'locks'/args.arm/args.key;lockpath.mkdir(parents=True,exist_ok=True)
 lock=(lockpath/f'{args.sample:03d}.lock').open('w')
 try:fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
 except BlockingIOError:
  print('Already running in another worker',args.arm,args.key,args.sample,flush=True);return
 affected={'1HMV_A_original', '7T9I_R', '7W80_B', '6N4B_R', '1SWX_A', '4L6R_A', '3U84_A', '2SHP_A'}
 while args.key in affected and (ROOT/'PAUSE_AFFECTED').exists():time.sleep(.2)
 done=ROOT/'outputs'/args.arm/args.key/f'sample_{args.sample:03d}.json'
 if done.exists() and (args.key not in affected or json.loads(done.read_text()).get('topology_revision')=='density-jumps-v1'):return
 import pyrosetta
 from pyrosetta.rosetta.protocols.relax import FastRelax
 from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
 from pyrosetta.rosetta.core.pack.task import TaskFactory
 from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline,RestrictToRepacking,IncludeCurrent
 from pyrosetta.rosetta.core.kinematics import MoveMap
 seed=20260915+int(hashlib.sha256(f'{args.key}:{args.arm}:{args.sample}'.encode()).hexdigest()[:7],16)
 flags=f'-mute all -constant_seed -jran {seed} -ex1 -ex2 -use_input_sc -ignore_zero_occupancy false -pack_missing_sidechains true -missing_density_to_jump true'
 pyrosetta.init(flags,silent=True)
 pose=pyrosetta.pose_from_pdb(str(ROOT/'inputs'/f'{args.key}.pdb'))
 seq=pose.sequence();sf=pyrosetta.create_score_function('ref2015')
 tf=TaskFactory();tf.push_back(InitializeFromCommandline());tf.push_back(RestrictToRepacking());tf.push_back(IncludeCurrent())
 before={(i,n):[pose.residue(i).xyz(n)[j] for j in range(3)] for i in range(1,pose.total_residue()+1) for n in ['N','CA','C','O'] if pose.residue(i).has(n)}
 t0=time.time();initial=float(sf(pose))
 if args.arm=='repack':
  task=tf.create_task_and_apply_taskoperations(pose);mover=PackRotamersMover(sf,task)
 else:
  mover=FastRelax(sf,5);mover.set_task_factory(tf)
  mm=MoveMap();mm.set_bb(True);mm.set_chi(True);mm.set_jump(False);mover.set_movemap(mm)
 mover.apply(pose)
 assert pose.sequence()==seq
 import numpy as np
 displacement=np.array([[pose.residue(i).xyz(n)[j]-x[j] for j in range(3)] for (i,n),x in before.items()])
 bbmax=float(np.linalg.norm(displacement,axis=1).max())
 if args.arm=='repack':assert bbmax<1e-6,bbmax
 out=ROOT/'outputs'/args.arm/args.key;out.mkdir(parents=True,exist_ok=True)
 target=out/f'sample_{args.sample:03d}.pdb';pose.dump_pdb(str(target))
 data={'topology_revision':'density-jumps-v1','arm':args.arm,'input_key':args.key,'sample_idx':args.sample,'seed':seed,'version':pyrosetta.version(),'flags':flags,'initial_score':initial,'final_score':float(sf(pose)),'seconds':time.time()-t0,'max_backbone_displacement':bbmax,'n_res':pose.total_residue(),'sequence_unchanged':True}
 checkpoint=target.with_suffix('.json.tmp');checkpoint.write_text(json.dumps(data,indent=2));os.replace(checkpoint,target.with_suffix('.json'));print(json.dumps(data),flush=True)
def launch_one(task):
 arm,key,i=task;out=ROOT/'outputs'/arm/key/f'sample_{i:03d}.json'
 if out.exists():return {'arm':arm,'key':key,'sample':i,'status':'cached'}
 logs=ROOT/'logs'/arm/key;logs.mkdir(parents=True,exist_ok=True)
 cmd=[sys.executable,str(Path(__file__).resolve()),'--worker','--arm',arm,'--key',key,'--sample',str(i)]
 with (logs/f'{i:03d}_{os.getpid()}.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 return {'arm':arm,'key':key,'sample':i,'returncode':r.returncode}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--worker',action='store_true');p.add_argument('--arm',choices=['repack','fastrelax']);p.add_argument('--key');p.add_argument('--sample',type=int);p.add_argument('--samples',type=int,default=30);p.add_argument('--start-sample',type=int,default=0);p.add_argument('--workers',type=int,default=24);args=p.parse_args()
 if args.worker:worker(args)
 else:
  rows=prepare();keys=list(dict.fromkeys(x['input_key'] for x in rows))
  tasks=[(arm,key,i) for i in range(args.start_sample,args.samples) for key in keys for arm in ['repack','fastrelax']]
  print('Tasks',len(tasks),'unique inputs',len(keys),flush=True)
  with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
   for result in pool.map(launch_one,tasks):print(json.dumps(result),flush=True)
