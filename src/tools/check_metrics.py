from models import MoleculeModel
from app.db import Session
import logging

from ddbox.metrics import metrics
from ddbox.metrics.utils import get_molecule_from_smiles_if_valid_or_none

logger = logging.getLogger(__name__)


logger = logging.getLogger(__name__)


session = Session()

query = session.query(MoleculeModel)

molecules = query \
    .order_by(MoleculeModel.id) \
    .offset(0) \
    .limit(1000) \
    .all()

smiles_a = [molecule.smiles for molecule in molecules]

molecules = query \
    .order_by(MoleculeModel.id) \
    .offset(100) \
    .limit(1000) \
    .all()

smiles_b = [molecule.smiles for molecule in molecules]


session.close()


def check():
    molecules_a = [get_molecule_from_smiles_if_valid_or_none(smiles) for smiles in smiles_a]
    molecules_b = [get_molecule_from_smiles_if_valid_or_none(smiles) for smiles in smiles_b]

    molecules_a = [molecule for molecule in molecules_a if molecule is not None]
    molecules_b = [molecule for molecule in molecules_b if molecule is not None]

    logger.info("FCD: %s" % (metrics.eval('FCD', molecules_a, molecules_b)))
    logger.info("SNN: %s" % (metrics.eval('SNN', molecules_a, molecules_b)))
    logger.info("Frag: %s" % (metrics.eval('Frag', molecules_a, molecules_b)))
    logger.info("Scaf: %s" % (metrics.eval('Scaf', molecules_a, molecules_b)))
    logger.info("IntDiv: %s" % (metrics.eval('IntDiv', molecules_a)))
    logger.info("IntDiv2: %s" % (metrics.eval('IntDiv2', molecules_a)))
    logger.info("fraction_valid: %s" % (metrics.eval('fraction_valid', smiles_a)))
    logger.info("fraction_unique: %s" % (metrics.eval('fraction_unique', smiles_a)))
    logger.info("fraction_passes: %s" % (metrics.eval('fraction_passes', molecules_a)))
    logger.info("DistributionDifferenceLogP: %s" % (metrics.eval('DistributionDifferenceLogP', molecules_a, molecules_b)))
    logger.info("DistributionDifferenceSA: %s" % (metrics.eval('DistributionDifferenceSA', molecules_a, molecules_b)))
    logger.info("DistributionDifferenceQED: %s" % (metrics.eval('DistributionDifferenceQED', molecules_a, molecules_b)))
    logger.info("DistributionDifferenceWeight: %s" % (metrics.eval('DistributionDifferenceWeight', molecules_a, molecules_b)))
