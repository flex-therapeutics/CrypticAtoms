from pathlib import Path
import json,numpy as np,pymol
pymol.finish_launching(['pymol','-qc']);from pymol import cmd
r=Path(__file__).parent;w=r.parent;rows=json.loads((r/'gallery_selection.json').read_text());targets={t['name']:t for t in json.loads((w/'fig2/targets.json').read_text())};aud=[]
for target in dict.fromkeys(q['target'] for q in rows):
 cmd.delete('all');cmd.load(str(w/'fig2/sessions'/f'{target}.pse'));view=cmd.get_view();cmd.disable('all');ext=[]
 for q in [q for q in rows if q['target']==target]:
  name=q['method'];cmd.load(q['protein'],name);cmd.load(q['ligand'],name+'_lig');cmd.remove(name+' and (hydro or not polymer.protein)');cmd.remove(name+'_lig and hydro')
  # Fit the protein and apply exactly the same rigid transform to its generated ligand.
  ref={int(a.resi):a.index for a in cmd.get_model('holo and name CA').atom if a.resi.lstrip('-').isdigit()};mob={int(a.resi):a.index for a in cmd.get_model(name+' and name CA').atom if a.resi.lstrip('-').isdigit()};common=[i for i in targets[target]['pocket']['residues_5A'] if i in ref and i in mob];assert len(common)>=3,(target,name)
  before=cmd.get_coords(name).copy();args=[]
  for i in common:args += [f'{name} and index {mob[i]}',f'holo and index {ref[i]}']
  fit=cmd.pair_fit(*args);after=cmd.get_coords(name);bc=before.mean(0);ac=after.mean(0);u,sv,vt=np.linalg.svd((before-bc).T@(after-ac));rot=u@vt
  assert np.sqrt(np.mean(np.sum(((before-bc)@rot+ac-after)**2,axis=1)))<.01
  xyz=cmd.get_coords(name+'_lig');cmd.load_coords((xyz-bc)@rot+ac,name+'_lig');lig=cmd.get_coords(name+'_lig');prot=cmd.get_coords(name);dist=np.sqrt(((lig[:,None,:]-prot[None,:,:])**2).sum(-1)).min();assert dist<6,(target,name,dist)
  ext.extend(lig.tolist());cmd.hide('everything',name+' or '+name+'_lig');cmd.disable(name);cmd.disable(name+'_lig');aud.append(dict(target=target,method=name,sample_idx=q['sample_idx'],fit_rmsd=fit,min_ligand_protein_distance=float(dist)))
 cmd.set_view(view);cmd.zoom(' or '.join(q['method']+'_lig' for q in rows if q['target']==target),5.5);camera=cmd.get_view()
 for q in [q for q in rows if q['target']==target]:
  name=q['method'];cmd.disable('all');cmd.enable(name);cmd.enable(name+'_lig');cmd.hide('everything');cmd.show('cartoon',name+' and byres ('+name+' within 8 of '+name+'_lig)');cmd.show('sticks',name+' and byres ('+name+' within 4.5 of '+name+'_lig)');cmd.show('sticks',name+'_lig');cmd.set_color('samplepink',[1,.6,.6]);cmd.set_color('ligandnavy',[.059,.278,.380]);cmd.color('samplepink',name);cmd.color('ligandnavy',name+'_lig');cmd.util.cnc(name+'_lig');cmd.set('cartoon_transparency',.35,name);cmd.set('stick_radius',.12,name);cmd.set('stick_radius',.24,name+'_lig');cmd.bg_color('white');cmd.set('ray_opaque_background',0);cmd.set('max_threads',8);cmd.set_view(camera);cmd.deselect();cmd.ray(1000,900);cmd.png(str(r/'renders'/f'{target}_{name}.png'),dpi=400);cmd.save(str(r/'sessions'/f'{target}_{name}.pse'));print('DONE',target,name,flush=True)
(r/'gallery_audit.json').write_text(json.dumps(aud,indent=2));cmd.quit()
