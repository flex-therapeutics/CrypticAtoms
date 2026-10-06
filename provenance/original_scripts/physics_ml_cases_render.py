from pathlib import Path
import json,csv,pymol,sys
pymol.finish_launching(['pymol','-qc']);from pymol import cmd
r=Path(__file__).parent;w=r.parent;t={x['name']:x for x in json.loads((w/'fig2/targets.json').read_text())};(r/'renders').mkdir(exist_ok=True);(r/'sessions').mkdir(exist_ok=True)
for target in (sys.argv[1:] or ['cAbl','HIV_RT__3V81','GLTP']):
 cmd.delete('all');cmd.load(str(w/'fig2/sessions'/f'{target}.pse'));view=cmd.get_view();cmd.delete('apo');cmd.delete('holo');cmd.delete('apo_front');cmd.delete('holo_front')
 for n in ['apo','holo','physics','learned']:cmd.load(str(r/'structures'/f'{target}_{n}.pdb'),n);cmd.dss(n)
 cmd.set_color('apogray',[.54,.54,.54]);cmd.set_color('holoteal',[.102,.6,.6]);cmd.set_color('samplepink',[1,.6,.6]);cmd.set_color('ligandnavy',[.059,.278,.380]);mov='+'.join(map(str,t[target]['pocket']['moving_residues_2A']));cmd.set_view(view);cmd.zoom('(holo and resi '+mov+') or ligand',3.0);camera=cmd.get_view()
 for mobile in ['apo','physics','learned']:
  cmd.disable('all');cmd.hide('everything');cmd.enable('holo');cmd.enable(mobile);cmd.enable('ligand');cmd.color('holoteal','holo');cmd.color('apogray' if mobile=='apo' else 'samplepink',mobile)
  for n in ['holo',mobile]:
   cmd.show('cartoon',n+' and byres ('+n+' within 8 of ligand)');cmd.show('sticks',n+' and resi '+mov);cmd.set('stick_radius',.19,n);cmd.set('cartoon_transparency',.65,n)
  cmd.show('sticks','ligand');cmd.color('ligandnavy','ligand');cmd.util.cnc('ligand');cmd.set('stick_radius',.22,'ligand');cmd.set('stick_transparency',.15,'holo');cmd.set('ray_opaque_background',0);cmd.set('max_threads',8);cmd.bg_color('white');cmd.set_view(camera);cmd.deselect();cmd.ray(1100,1100);cmd.png(str(r/'renders'/f'{target}_{mobile}.png'),dpi=400);cmd.save(str(r/'sessions'/f'{target}_{mobile}.pse'));print('DONE',target,mobile,flush=True)
cmd.quit()
