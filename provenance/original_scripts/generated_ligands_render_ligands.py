from pathlib import Path
import json,pymol
pymol.finish_launching(['pymol','-qc']);from pymol import cmd
r=Path(__file__).parent
for q in json.loads((r/'gallery_selection.json').read_text()):
 t=q['target'];m=q['method'];cmd.delete('all');cmd.load(str(r/'sessions'/f'{t}_{m}.pse'));cmd.disable('all');cmd.enable(m+'_lig');cmd.hide('everything');cmd.show('sticks',m+'_lig');cmd.zoom(m+'_lig',.7);cmd.set('max_threads',8);cmd.ray(1000,700);cmd.png(str(r/'renders'/f'{t}_{m}_ligand.png'),dpi=400);cmd.save(str(r/'sessions'/f'{t}_{m}_ligand.pse'));print('DONE',t,m,flush=True)
cmd.quit()
