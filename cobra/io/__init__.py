from warnings import warn

from .sbml3 import read_sbml_model, write_sbml_model
from .json import load_json_model, save_json_model, to_json

# These functions have other dependencies
try:
    import libsbml
    from .sbml import read_legacy_sbml
    from .sbml import write_cobra_model_to_sbml_file as write_legacy_sbml
except ImportError:
    warn("cobra.io.sbml requires libsbml")
    libsbml = None
    read_legacy_sbml = None
    write_legacy_sbml = None

try:
    import scipy
    from .mat import load_matlab_model, save_matlab_model
except ImportError:
    warn("cobra.io.mat requires scipy")
    scipy = None
    load_matlab_model = None
    save_matlab_model = None

try:
    import requests
except ImportError:
    warn("cobra.io.modelseed requires requests")
    requests = None
else:
    from .modelseed import reconstruct_modelseed_model, gapfill_modelseed_model, optimize_modelseed_model
    from .modelseed import delete_modelseed_model, get_modelseed_model_data, get_modelseed_model_stats
    from .modelseed import get_modelseed_gapfill_solutions, get_modelseed_fba_solutions, list_modelseed_models
    from .modelseed import get_workspace_object_meta, get_workspace_object_data, list_workspace_objects
    from .modelseed import create_cobra_model_from_modelseed_model
    from .PatricClient import get_patric_token

del libsbml, scipy, requests, warn
