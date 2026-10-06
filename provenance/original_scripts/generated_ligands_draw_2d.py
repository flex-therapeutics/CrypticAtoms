"""Depict the existing gallery ligand graphs without repairing or selecting new samples."""
from pathlib import Path
import json
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
root=Path(__file__).resolve().parent
out=root/'renders';out.mkdir(exist_ok=True)
audit=[]
selection=root/'gallery_selection.json'
if not selection.exists(): selection=root/'generated_ligands_gallery_selection.json'
for row in json.loads(selection.read_text()):
    source=root/'ligand_2d_sources'/f"{row['target']}_{row['method']}.sdf"
    mol=Chem.SDMolSupplier(str(source),sanitize=True,removeHs=True)[0]
    assert mol is not None,row
    original=Chem.MolToSmiles(mol)
    drawing=Chem.Mol(mol)
    rdDepictor.Compute2DCoords(drawing,canonOrient=True,clearConfs=True)
    prepared=rdMolDraw2D.PrepareMolForDrawing(drawing,addChiralHs=False)
    stem=out/f"{row['target']}_{row['method']}_ligand_2d"
    for kind in ('png','svg'):
        drawer=rdMolDraw2D.MolDraw2DCairo(900,900) if kind=='png' else rdMolDraw2D.MolDraw2DSVG(900,900)
        opt=drawer.drawOptions();opt.clearBackground=False;opt.padding=.04
        opt.bondLineWidth=7;opt.minFontSize=50;opt.maxFontSize=68
        opt.updateAtomPalette({6:(.07,.27,.33)})
        drawer.DrawMolecule(prepared);drawer.FinishDrawing()
        data=drawer.GetDrawingText()
        if kind=='png':Path(str(stem)+'.png').write_bytes(data)
        else:Path(str(stem)+'.svg').write_text(data)
    assert Chem.MolToSmiles(mol)==original
    audit.append({'target':row['target'],'method':row['method'],'sample_idx':row['sample_idx'],'smiles':original,'heavy_atoms':mol.GetNumHeavyAtoms(),'fragments':len(Chem.GetMolFrags(mol))})
(root/'ligand_2d_audit.json').write_text(json.dumps({'rdkit_version':rdBase.rdkitVersion,'source':'Same SDF graphs and gallery samples as the existing 3D views; no graph repairs.','examples':audit},indent=2))
print('Rendered',len(audit),'unchanged ligand graphs in 2D.')
