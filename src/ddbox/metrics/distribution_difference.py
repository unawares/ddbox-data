from typing import List

from ddbox.metrics.SA_Score import sascorer
from ddbox.registries import metrics
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import qed
from scipy.stats import wasserstein_distance


@metrics.register('DistributionDifferenceLogP')
def distribution_difference_logp(generated_molecules: List[Chem.rdchem.Mol], reference_molecules: List[Chem.rdchem.Mol]):
    generated_molecules_values = [Chem.Crippen.MolLogP(molecule) for molecule in generated_molecules]
    reference_molecules_values = [Chem.Crippen.MolLogP(molecule) for molecule in reference_molecules]
    return wasserstein_distance(generated_molecules_values, reference_molecules_values)


@metrics.register('DistributionDifferenceSA')
def distribution_difference_sa(generated_molecules: List[Chem.rdchem.Mol], reference_molecules: List[Chem.rdchem.Mol]):
    generated_molecules_values = [sascorer.calculateScore(molecule) for molecule in generated_molecules]
    reference_molecules_values = [sascorer.calculateScore(molecule) for molecule in reference_molecules]
    return wasserstein_distance(generated_molecules_values, reference_molecules_values)


@metrics.register('DistributionDifferenceQED')
def distribution_difference_sa(generated_molecules: List[Chem.rdchem.Mol], reference_molecules: List[Chem.rdchem.Mol]):
    generated_molecules_values = [qed(molecule) for molecule in generated_molecules]
    reference_molecules_values = [qed(molecule) for molecule in reference_molecules]
    return wasserstein_distance(generated_molecules_values, reference_molecules_values)


@metrics.register('DistributionDifferenceWeight')
def distribution_difference_sa(generated_molecules: List[Chem.rdchem.Mol], reference_molecules: List[Chem.rdchem.Mol]):
    generated_molecules_values = [Descriptors.MolWt(molecule) for molecule in generated_molecules]
    reference_molecules_values = [Descriptors.MolWt(molecule) for molecule in reference_molecules]
    return wasserstein_distance(generated_molecules_values, reference_molecules_values)
