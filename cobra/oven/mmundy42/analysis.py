from six import iteritems, print_
import pandas as pd
from tabulate import tabulate


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
        model
            Model to analyze (must be optimized)
        threshold : float)
            Tolerance for determining if a flux is zero
        floatfmt : str
            Format for floats when printed by tabulate
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

    Parameters
    ----------
        model
            Model to analyze (must be optimized)
        threshold : float, optional
            Tolerance for determining if a flux is zero
        floatfmt : str, optional
            Format for floats when printed by tabulate

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

    Parameters
    ----------
        model
            Model to analyze
    """

    boundary_reactions = model.reactions.query(lambda x: x, 'boundary')
    output = [[rxn.id, rxn.name, rxn.reaction] for rxn in boundary_reactions]
    print_(tabulate(output, tablefmt='simple', headers=['ID', 'NAME', 'DEFINITION']))

    return

def compare_models():
    parser = argparse.ArgumentParser(prog='compare-models')
    parser.add_argument('model1', help='path to sbml file of first model')
    parser.add_argument('model2', help='path to sbml file of second model')
    parser.add_argument('--source1', help='source system that generated first model', default='modelbuilder')
    parser.add_argument('--source2', help='source system that generated second model', default='modelseed')
    parser.add_argument('--detail', help='show details on reaction and metabolite differences', action='store_true', default=False)
    parser.add_argument('--boundary', help='compare boundary reactions', action='store_true', default=False)
    args = parser.parse_args()

    # Read the SBML files for the two models to compare.
    print 'Reading first model from {0}'.format(args.model1)
    model1 = cobra.io.read_sbml_model(args.model1)
    print 'Reading second model from {0}'.format(args.model2)
    model2 = cobra.io.read_sbml_model(args.model2)

    # A model created by the ModelSEED service includes a compartment index number at the end of
    # every reaction ID and metabolite ID. Remove the compartment index number so comparisons by
    # ID do not adjustment later.
    if args.source1 == 'modelseed':
        trim_index(model1)

    if args.source2 == 'modelseed':
        trim_index(model2)
        # trim_index = re.compile(r'[\d]+$')
        # for r2 in model2.reactions:
        #     r2.id = re.sub(trim_index, '', r2.id)
        # model2.repair()

    # Compare reactions.
    print '{0} reactions in first model'.format(len(model1.reactions))
    print '{0} reactions in second model'.format(len(model2.reactions))

    # See if reactions from first model are in the second model.
    num_matched = 0
    num_only_in_one = 0
    for r1 in model1.reactions:
        try:
            lookup = r1.id #+'0' # Need to adjust based on source options
            r2 = model2.reactions.get_by_id(lookup)
            num_matched += 1
        except KeyError:
            num_only_in_one += 1
            if args.detail:
                print 'Reaction {0} is not in second model'.format(lookup)
    print '{0} reactions in both models'.format(num_matched)
    print '{0} reactions only in first'.format(num_only_in_one)

    # See if reactions from second model are in the first model.
    num_matched = 0
    num_only_in_two = 0
    for r2 in model2.reactions:
        try:
            lookup = r2.id[:-1] # Need to adjust based on source options
            r1 = model1.reactions.get_by_id(lookup)
            num_matched += 1
        except KeyError:
            num_only_in_two += 1
            if args.detail:
                print 'Reaction {0} is not in first model'.format(lookup)
    print '{0} reactions in both models'.format(num_matched)
    print '{0} reactions only in second'.format(num_only_in_two)

    # Compare metabolites.
    print '{0} metabolites in first model'.format(len(model1.metabolites))
    print '{0} metabolites in second model'.format(len(model2.metabolites))

    # See if metabolites from first model are in the second model.
    num_matched = 0
    num_only_in_one = 0
    for m1 in model1.metabolites:
        try:
            if m1.id.endswith('_b'): # Need to adjust based on source options
                lookup = m1.id
            else:
                lookup = m1.id+'0'
            m2 = model2.metabolites.get_by_id(lookup)
            num_matched += 1
        except KeyError:
            num_only_in_one += 1
            if args.detail:
                print 'Metabolite {0} is not in second model'.format(lookup)
    print '{0} metabolites in both models'.format(num_matched)
    print '{0} metabolites only in first'.format(num_only_in_one)

    # See if metabolites from second model are in the first model.
    num_matched = 0
    num_only_in_two = 0
    for m2 in model2.metabolites:
        try:
            if m2.id.endswith('_b'): # Need to adjust based on source options
                lookup = m2.id
            else:
                lookup = m2.id[:-1]
            m1 = model1.metabolites.get_by_id(lookup)
            num_matched += 1
        except KeyError:
            num_only_in_two += 1
            if args.detail:
                print 'Metabolite {0} is not in first model'.format(lookup)
    print '{0} metabolites in both models'.format(num_matched)
    print '{0} metabolites only in second'.format(num_only_in_two)

    # See about system boundary reactions.
    if args.boundary:
        # Get the list of system boundary reactions from first model.
        model1_boundary = list()
        for r1 in model1.reactions:
            if r1.boundary == 'system_boundary':
                model1_boundary.append(r1)
        print '{0} reactions are system boundary reactions in first model'.format(len(model1_boundary))

        # Get the list of system boundary reactions from second model.
        model2_boundary = list()
        for r2 in model2.reactions:
            if r2.boundary == 'system_boundary':
                model2_boundary.append(r1)
        print '{0} reactions are system boundary reactions in second model'.format(len(model2_boundary))

    return
