import csv,json,unittest
from pathlib import Path
import numpy as np
from benchmark.metrics import load_protein_chain,superimpose_on_pocket_ca,all_atom_rmsd_over_residues
from benchmark.score import score
ROOT=Path(__file__).resolve().parents[1]
class MetricRegression(unittest.TestCase):
 def test_all_reference_baselines(self):
  targets=json.loads((ROOT/'benchmark/targets.json').read_text())
  self.assertEqual(len(targets),21)
  for t in targets:
   with self.subTest(target=t['name']):
    a=load_protein_chain(ROOT/t['apo']['path'],t['apo']['chain']);h=load_protein_chain(ROOT/t['holo']['path'],t['holo']['chain'])
    ah,_,_=superimpose_on_pocket_ca(a,h,t['pocket']['residues_5A']);d=all_atom_rmsd_over_residues(ah,h,t['pocket']['moving_residues_2A'])
    self.assertAlmostEqual(d,t['baselines']['align_pocket']['aa_rmsd_moving'],delta=.02)
    # A rigid-body transformation cannot change the pocket-aligned metric.
    shifted=a.copy();shifted.coord=a.coord@np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])+[10.,-4.,7.]
    shifted,_,_=superimpose_on_pocket_ca(shifted,h,t['pocket']['residues_5A']);v=all_atom_rmsd_over_residues(shifted,h,t['pocket']['moving_residues_2A'])
    self.assertAlmostEqual(d,v,delta=1e-4)
 def test_figure4_examples(self):
  for q in csv.DictReader((ROOT/'examples/manifest.csv').read_text().splitlines()):
   with self.subTest(target=q['target'],group=q['group']):
    result=score(q['target'],ROOT/q['path'],q['chain'])
    for key in ['Dsa','Dsh','bb_rmsd','sc_rmsd']:
     self.assertAlmostEqual(result[key],float(q[key]),delta=.02)
    self.assertEqual(result['tier1'],q['tier1'].lower()=='true')
if __name__=='__main__':unittest.main()
