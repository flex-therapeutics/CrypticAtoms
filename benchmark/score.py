"""Score a conformation already mapped to the reference residue numbering."""
import argparse,json
from pathlib import Path
import numpy as np
from .metrics import load_protein_chain,superimpose_on_pocket_ca,all_atom_rmsd_over_residues
from .diagnostics import stability,pocket_clashes,backbone_sidechain_sums,residue
ROOT=Path(__file__).resolve().parents[1]
def score(target, sample, chain, residue_map=None):
    targets={t['name']:t for t in json.loads((ROOT/'benchmark/targets.json').read_text())};t=targets[target]
    a=load_protein_chain(ROOT/t['apo']['path'],t['apo']['chain']);h=load_protein_chain(ROOT/t['holo']['path'],t['holo']['chain']);x=load_protein_chain(sample,chain)
    if residue_map is not None:
        assert set(map(str,x.res_id))<=set(residue_map),'residue map must cover the selected chain'
        x.res_id=np.array([int(residue_map[str(i)]) for i in x.res_id])
    assert len(set(zip(x.res_id,x.atom_name)))==len(x),'duplicate residue/atom keys'
    pocket=t['pocket']['residues_5A'];moving=t['pocket']['moving_residues_2A'];D=t['baselines']['align_pocket']['aa_rmsd_moving']
    xa,_,na=superimpose_on_pocket_ca(x,a,pocket);xh,_,nh=superimpose_on_pocket_ca(x,h,pocket);ha,_,_=superimpose_on_pocket_ca(h,a,pocket)
    da=all_atom_rmsd_over_residues(xa,a,moving);dh=all_atom_rmsd_over_residues(xh,h,moving)
    if not np.isfinite([da,dh,D]).all():raise ValueError('No finite score for the moving-residue set')
    tier1=bool(dh<D and dh<da)
    out=dict(target=target,Dsa=da,Dsh=dh,D=D,tier1=tier1,tier2=bool(tier1 and dh<1.5),tier3=bool(tier1 and dh<1),n_align_apo=na,n_align_holo=nh,missing_pocket_residues=sorted(set(pocket)-set(x.res_id)),moving_residues_present=sorted(set(moving)&set(x.res_id)))
    out.update(stability(xa,xh,a,h,ha,pocket));out.update(pocket_clashes(x,a,h,pocket))
    complete_ids=[r for r in moving if set(residue(a,r).atom_name)==set(residue(h,r).atom_name) and set(residue(a,r).res_name)==set(residue(h,r).res_name)]
    diagnostic_ids=[r for r in complete_ids if set(residue(xh,r).atom_name)==set(residue(h,r).atom_name) and set(residue(xh,r).res_name)==set(residue(h,r).res_name)]
    out['diagnostic_residues']=diagnostic_ids
    vals=backbone_sidechain_sums(xh,h,diagnostic_ids)
    for prefix in ['bb','sc']:
        n=sum(z[prefix+'_n'] for z in vals);out[prefix+'_rmsd']=float(np.sqrt(sum(z[prefix+'_sum'] for z in vals)/n)) if n else None
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--target',required=True);p.add_argument('--sample',type=Path,required=True);p.add_argument('--chain',required=True);p.add_argument('--residue-map',type=Path);p.add_argument('--output',type=Path);args=p.parse_args()
    result=score(args.target,args.sample,args.chain,json.loads(args.residue_map.read_text()) if args.residue_map else None)
    text=json.dumps(result,indent=2,default=lambda x:x.item())+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')
