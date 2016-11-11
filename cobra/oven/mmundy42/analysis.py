from six import iteritems, print_
import pandas as pd
from tabulate import tabulate
from cobra.core import DictList

"""
Notes
-----
    For ModelSEED models, I don't understand the boundary reactions rxn13782_c, rxn13783_c,
    and rxn13784_c which I think are sink reactions for protein biosynthesis, DNA replication,
    and RNA transcrption.  The compounds show up as being consumed but are produced by the
    biomass reaction.
"""


def format_long_string(string, max_length):
    """ Format a string so it fits in column of a specific width.

    Parameters
    ----------
    string : str
        String to format
    max_length : int
        Maximum length of returned string

    Returns
    -------
    str
        Formated string
    """

    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string


def metabolites_consumed(model, threshold=1E-8, floatfmt='.3f'):
    """ Print the metabolites consumed by the organism with IDs and names.

    Parameters
    ----------
    model : cobra.Model
        Model to analyze (must be optimized)
    threshold : float, optional
        Tolerance for determining if a flux is zero
    floatfmt : str, optional
        Format for floats when printed by tabulate
    """

    # Build a dictionary of metabolite consumption and production from the boundary reactions.
    reactions = model.reactions.query(lambda x: x, 'boundary')
    metabolite_fluxes = {}
    for rxn in reactions:
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
    print_('{0} metabolites consumed by {1}\n'.format(len(consumed.index), model.id))
    print_(tabulate(consumed.loc[:, ['id', 'flux', 'name']].values,
                    floatfmt=floatfmt, tablefmt='simple', headers=['ID', 'FLUX', 'NAME']))

    return


def metabolites_produced(model, threshold=1E-8, floatfmt='.3f'):
    """ Print the metabolites produced by the organism with IDs and names.

    Parameters
    ----------
    model : cobra.Model
        Model to analyze (must be optimized)
    threshold : float, optional
        Tolerance for determining if a flux is zero
    floatfmt : str, optional
        Format for floats when printed by tabulate

    """

    # Build a dictionary of metabolite consumption and production from the boundary reactions.
    reactions = model.reactions.query(lambda x: x, 'boundary')
    metabolite_fluxes = {}
    for rxn in reactions:
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


def boundary_reactions(model):
    """ Print the boundary reactions in a model with id, name, and definition.

    Parameters
    ----------
    model : cobra.Model
        Model to analyze
    """

    reactions = model.reactions.query(lambda x: x, 'boundary')
    reactions.sort(key=lambda x: x.id)
    print_('{0} boundary reactions in {1}\n'.format(len(reactions), model.id))
    output = [[rxn.id, rxn.name, rxn.reaction] for rxn in reactions]
    print_(tabulate(output, tablefmt='simple', headers=['ID', 'NAME', 'DEFINITION']))

    return


def exchange_reactions(model):
    """ Print the exchange reactions in a model with id, name, and definition.

        An exchange reaction is defined as a reaction with an ID that starts with 'EX_'.

    Parameters
    ----------
    model : cobra.Model
        Model to analyze
    """

    reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
    reactions.sort(key=lambda x: x.id)
    print_('{0} exchange reactions in {1}\n'.format(len(reactions), model.id))
    output = [[rxn.id, rxn.name, rxn.reaction] for rxn in reactions]
    print_(tabulate(output, tablefmt='simple', headers=['ID', 'NAME', 'DEFINITION']))

    return


def compare_models(model1, model2, detail=False, boundary=False):
    """ Compare two models and report differences.

    Parameters
    ----------
    model1 : cobra.Model
        First model to analyze
    model2 : cobra.Model
        Second model to analyze
    detail : bool, optional
        When true, print details on differences
    boundary : bool, optional
        When true, print info about boundary reactions
    """

    header = ['ID', 'NAME']

    # Compare reactions.
    print_('REACTIONS\n' + '---------')
    print_('{0} reactions in {1}'.format(len(model1.reactions), model1.id))
    print_('{0} reactions in {1}\n'.format(len(model2.reactions), model2.id))

    # See if reactions from first model are in the second model.
    num_matched = 0
    reaction_only_in_one = DictList()
    for r1 in model1.reactions:
        if model2.reactions.has_id(r1.id):
            num_matched += 1
        else:
            reaction_only_in_one.append(r1)
    print_('{0} reactions in {1} and {2}'.format(num_matched, model1.id, model2.id))
    print_('{0} reactions only in {1}\n'.format(len(reaction_only_in_one), model1.id))
    if detail:
        output = [[rxn.id, rxn.name] for rxn in reaction_only_in_one]
        print_(tabulate(output, tablefmt='simple', headers=header) + '\n')

    # See if reactions from second model are in the first model.
    num_matched = 0
    reaction_only_in_two = DictList()
    for r2 in model2.reactions:
        if model1.reactions.has_id(r2.id):
            num_matched += 1
        else:
            reaction_only_in_two.append(r2)
    print_('{0} reactions in both {1} and {2}'.format(num_matched, model2.id, model1.id))
    print_('{0} reactions only in {1}\n'.format(len(reaction_only_in_two), model2.id))
    if detail:
        output = [[rxn.id, rxn.name] for rxn in reaction_only_in_two]
        print_(tabulate(output, tablefmt='simple', headers=header) + '\n')

    # Compare metabolites.
    print_('METABOLITES\n' + '-----------')
    print_('{0} metabolites in {1}'.format(len(model1.metabolites), model1.id))
    print_('{0} metabolites in {1}\n'.format(len(model2.metabolites), model2.id))

    # See if metabolites from first model are in the second model.
    num_matched = 0
    metabolite_only_in_one = DictList()
    for m1 in model1.metabolites:
        if model2.metabolites.has_id(m1.id):
            num_matched += 1
        else:
            metabolite_only_in_one.append(m1)
    print_('{0} metabolites in both {1} and {2}'.format(num_matched, model1.id, model2.id))
    print_('{0} metabolites only in {1}\n'.format(len(metabolite_only_in_one), model1.id))
    if detail:
        output = [[met.id, met.name] for met in metabolite_only_in_one]
        print_(tabulate(output, tablefmt='simple', headers=header) + '\n')

    # See if metabolites from second model are in the first model.
    num_matched = 0
    metabolite_only_in_two = DictList()
    for m2 in model2.metabolites:
        if model1.metabolites.has_id(m2.id):
            num_matched += 1
        else:
            metabolite_only_in_two.append(m2)
    print_('{0} metabolites in both {1} and {2}'.format(num_matched, model2.id, model1.id))
    print_('{0} metabolites only in {1}\n'.format(len(metabolite_only_in_two), model2.id))
    if detail:
        output = [[met.id, met.name] for met in metabolite_only_in_two]
        print_(tabulate(output, tablefmt='simple', headers=header) + '\n')

    # See about system boundary reactions.
    if boundary:
        # Get the list of system boundary reactions from first model.
        model1_boundary = model1.reactions.query(lambda x: x, 'boundary')
        print_('{0} reactions are system boundary reactions in {1}'.format(len(model1_boundary), model1.id))

        # Get the list of system boundary reactions from second model.
        model2_boundary = model2.reactions.query(lambda x: x, 'boundary')
        print_('{0} reactions are system boundary reactions in {1}'.format(len(model2_boundary), model2.id))

    return
