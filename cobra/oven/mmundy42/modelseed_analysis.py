from six import print_
import pandas as pd
from tabulate import tabulate
import re
import json
from warnings import warn

from cobra.io.modelseed import modelseed_suffix_re, get_workspace_object_data


def modelseed_metabolites_consumed(fba_solution, threshold=1E-8, floatfmt='.3f'):
    """ Print the metabolites consumed by the organism with IDs from a ModelSEED FBA solution.

    Parameters
    ----------
    fba_solution : dict
        Dictionary of ModelSEED FBA solution
    threshold : float, optional
        Tolerance for determining if a flux is zero
    floatfmt : str, optional
        Format for floats when printed by tabulate
    """

    # Build a dictionary of metabolite consumption from the exchanges in the solution.
    metabolite_fluxes = {}
    for met_id in fba_solution['exchanges']:
        metabolite = fba_solution['exchanges'][met_id]
        if metabolite['x'] >= threshold:
            metabolite_fluxes[met_id] = {
                'id': met_id,
                'flux': metabolite['x']}

    # Build a data frame from the dictionary.
    consumed = pd.DataFrame(metabolite_fluxes).T

    # Sort by flux value.
    consumed = consumed.sort_values(by=['flux', 'id'], ascending=[False, True])

    # Print a nice table.
    print_('{0} metabolites consumed\n'.format(len(consumed.index)))
    print_(tabulate(consumed.loc[:, ['id', 'flux']].values,
                    floatfmt=floatfmt, tablefmt='simple', headers=['ID', 'FLUX']))

    return


def modelseed_metabolites_produced(fba_solution, threshold=1E-8, floatfmt='.3f'):
    """ Print the metabolites produced by the organism with IDs from a ModelSEED FBA solution.

    Parameters
    ----------
    fba_solution : dict
        Dictionary of ModelSEED FBA solution
    threshold : float, optional
        Tolerance for determining if a flux is zero
    floatfmt : str, optional
        Format for floats when printed by tabulate
    """

    # Build a dictionary of metabolite consumption from the exchanges in the solution.
    metabolite_fluxes = {}
    for met_id in fba_solution['exchanges']:
        metabolite = fba_solution['exchanges'][met_id]
        if metabolite['x'] <= -threshold:
            metabolite_fluxes[met_id] = {
                'id': met_id,
                'flux': abs(metabolite['x'])}

    # Build a data frame from the dictionary.
    produced = pd.DataFrame(metabolite_fluxes).T

    # Sort by flux value.
    produced = produced.sort_values(by=['flux', 'id'], ascending=[False, True])

    # Print a nice table.
    print_('{0} metabolites produced\n'.format(len(produced.index)))
    print_(tabulate(produced.loc[:, ['id', 'flux']].values,
                    floatfmt=floatfmt, tablefmt='simple', headers=['ID', 'FLUX']))

    return


def modelseed_remove_compartment_index(model):
    """ Remove the compartment index number from reaction and metabolite IDs.

    Parameters
    ----------
    model : cobra.Model
        Model created from a ModelSEED model (usually by importing SBML)

    """

    # Remove the compartment index from all of the reaction IDs.
    for r in model.reactions:
        r.id = re.sub(modelseed_suffix_re, r'_\g<1>', r.id)
    model.repair()

    # Remove the compartment index from all of the metabolite IDs.
    for m in model.metabolites:
        m.id = re.sub(modelseed_suffix_re, r'_\g<1>', m.id)
    model.repair()

    return


def modelseed_convert_media(media_reference, media_filename):

    # Get the ModelSEED media object from the workspace.
    source_media = get_workspace_object_data(media_reference, json_data=False)
    lines = source_media.split('\n')
    header = lines[0].split('\t')
    if header[0] != 'id' or header[3] != 'minflux' or header[4] != 'maxflux':
        warn('Media must have ID in first column, minflux in fourth column, and maxflux in fifth column')
        return

    # Cobrapy media is expressed in terms of exchange reaction IDs.
    cobra_media = dict()
    for index in range(1, len(lines)):
        line = lines[index]
        if len(line) > 0:  # Skip empty lines
            fields = line.split('\t')
            if len(fields) < 5:
                warn('Skipped line {0} because fields are missing: {1}'.format(index, line))
                continue
            exchange_id = 'EX_{0}_e'.format(fields[0])
            cobra_media[exchange_id] = (float(fields[3]), float(fields[4]))

    # Write JSON file.
    json.dump(cobra_media, open(media_filename, 'w'))
    return
