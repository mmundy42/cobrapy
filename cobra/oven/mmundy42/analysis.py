from six import iteritems, print_
import pandas as pd
from tabulate import tabulate


# Notes:
#   For ModelSEED models, I don't understand the boundary reactions rxn13782_c, rxn13783_c,
#   and rxn13784_c which I think are sink reactions for protein biosynthesis, DNA replication,
#   and RNA transcrption.  The compounds show up as being consumed but are produced by the
#   biomass reaction.

def format_long_string(string, max_length):
    """ Format a string so it fits in column of a specific width.

    Args:
        string: String to format
        max_length: Maximum length of returned string

    Returns:
        Formated string
    """

    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string


def metabolites_consumed(model, threshold=1E-8, floatfmt='.3f'):
    """ Print the metabolites consumed by the organism with IDs and names.

    Args:
        model: Model to analyze (must be optimized)
        threshold (float): Tolerance for determining if a flux is zero
        floatfmt (str): Format for floats when printed by tabulate

    Returns:
        Nothing
    """

    # Build a dictionary of metabolite consumption and production from the boundary reactions.
    boundary_reactions = model.reactions.query(lambda x: x, 'boundary')
    metabolite_fluxes = {}
    for rxn in boundary_reactions:
        for met, stoich in iteritems(rxn.metabolites):
            metabolite_fluxes[met] = {
                'id': format_long_string(met.id, 15),
                'name': met.name,
                'flux': stoich * rxn.x}

    # Build a data frame to find metabolites consumed.
    metabolite_fluxes = pd.DataFrame(metabolite_fluxes).T

    # Drop unused boundary fluxes.
    metabolite_fluxes = metabolite_fluxes[metabolite_fluxes.flux.abs() > threshold].copy()

    # Keep positive flux values which mean metabolite is consumed.
    # @todo This does not support fva.
    consumed = metabolite_fluxes[metabolite_fluxes.flux > 0.].copy()

    # Sort by flux value.
    consumed = consumed.sort_values(by=['flux', 'id'], ascending=[False, True])

    # Print a nice table.
    print_(tabulate(consumed.loc[:, ['id', 'flux', 'name']].values,
                    floatfmt=floatfmt, tablefmt='simple', headers=['ID', 'FLUX', 'NAME']))

    return


def metabolites_produced(model, threshold=1E-8, floatfmt='.3f'):
    """ Print the metabolites produced by the organism with IDs and names.

    Args:
        model: Model to analyze (must be optimized)
        threshold (float): Tolerance for determining if a flux is zero
        floatfmt (str): Format for floats when printed by tabulate

    Returns:
        Nothing
    """

    # Build a dictionary of metabolite consumption and production from the boundary reactions.
    boundary_reactions = model.reactions.query(lambda x: x, 'boundary')
    metabolite_fluxes = {}
    for rxn in boundary_reactions:
        for met, stoich in iteritems(rxn.metabolites):
            metabolite_fluxes[met] = {
                'id': format_long_string(met.id, 15),
                'name': met.name,
                'flux': stoich * rxn.x}

    # Build a data frame to find metabolites produced.
    metabolite_fluxes = pd.DataFrame(metabolite_fluxes).T

    # Drop unused boundary fluxes.
    metabolite_fluxes = metabolite_fluxes[metabolite_fluxes.flux.abs() > threshold].copy()

    # Keep negative flux values which mean metabolite is produced.
    # @todo This does not support fva.
    produced = metabolite_fluxes[metabolite_fluxes.flux < 0.].copy()

    # Make flux value positive and sort by flux value.
    produced.flux = produced.flux.abs()
    produced = produced.sort_values(by=['flux', 'id'], ascending=[False, True])

    # Print a nice table.
    print_(tabulate(produced.loc[:, ['id', 'flux', 'name']].values,
                    floatfmt=floatfmt, tablefmt='simple', headers=['ID', 'FLUX', 'NAME']))

    return


def exchange_reactions(model):
    """ Print the exchange reactions in a model with id, name, and definition.

    Args:
        model: Model to analyze

    Returns:
        Nothing
    """

    boundary_reactions = model.reactions.query(lambda x: x, 'boundary')
    output = [[rxn.id, rxn.name, rxn.reaction] for rxn in boundary_reactions]
    print_(tabulate(output, tablefmt='simple', headers=['ID', 'NAME', 'DEFINITION']))

    return
