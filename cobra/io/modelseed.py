from __future__ import absolute_import
from time import sleep
from operator import itemgetter
from warnings import warn
from os.path import join
import re
import json

from .. import Model, Reaction, Metabolite, Gene

from . import PatricClient

# ModelSEED service endpoint
modelseed_url = 'https://p3.theseed.org/services/ProbModelSEED'

# Workspace service endpoint
workspace_url = 'https://p3.theseed.org/services/Workspace'

# Regular expression for compartment suffix on ModelSEED IDs
modelseed_suffix_re = re.compile(r'_([a-z])\d*$')

# Regular expression for prefix on PATRIC gene IDs
patric_gene_prefix_re = re.compile(r'^fig\|')


def delete_modelseed_model(reference):
    """ Delete a ModelSEED model from the workspace.

        Parameters
        ----------
            reference : str
                Workspace reference to model
    """

    try:
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        ms_client.call('delete_model', {'model': reference})
    except PatricClient.ServerError as e:
        handle_server_error(e, [reference])

    return


def gapfill_modelseed_model(reference, media_reference=None, likelihood=False, comprehensive=False, solver=None):
    """ Run gap fill on a ModelSEED model.

        @param reference: Workspace reference to model
        @param media_reference: Workspace reference to media to gap fill on (default is complete media)
        @param likelihood: True to use likelihood-based gap fill
        @param comprehensive: True to run a comprehensive gap fill
        @param solver: Name of LP solver
        @return: Model statistics dictionary
    """

    params = dict()
    params['model'] = reference
    params['integrate_solution'] = 1
    if media_reference is not None:
        params['media'] = media_reference
    if likelihood:
        params['probanno'] = 1
    else:
        params['probanno'] = 0
    if comprehensive:
        params['comprehensive_gapfill'] = 1
    if solver is not None:
        params['solver'] = solver

    try:
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        job_id = ms_client.call('GapfillModel', params)
        _wait_for_job(job_id)
    except PatricClient.ServerError as e:
        references = [reference]
        if media_reference is not None:
            references.append(media_reference)
        handle_server_error(e, references)

    return get_modelseed_model_stats(reference)


def get_modelseed_fba_solutions(reference):
    """ Get the list of fba solutions available for a ModelSEED model.

        @param reference: Workspace reference to model
        @return: List of fba solution data structures
    """

    solutions = None
    try:
        get_modelseed_model_stats(reference)  # Confirm model exists
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        solutions = ms_client.call('list_fba_studies', {'model': reference})
    except PatricClient.ServerError as e:
        handle_server_error(e, [reference])

    # For each solution in the list, get the referenced fba object, and add the
    # details on the flux values to the solution. Note ModelSEED stores the
    # results of each flux balance analysis separately.
    for sol in solutions:
        try:
            solution_data = get_workspace_object_data(sol['ref'])
        except PatricClient.ServerError as e:
            handle_server_error(e, sol['ref'])

        # A ModelSEED model does not have exchange reactions so instead the results of a flux
        # balance analysis reports flux values on metabolites in the extracellular compartment.
        # For ModelSEED, a positive flux means the metabolite is consumed and a negative flux
        # means the metabolite is produced.
        # @todo Should the fluxes be flipped to match COBRA models?
        sol['exchanges'] = dict()
        for flux in solution_data['FBACompoundVariables']:
            exchange_id = flux['modelcompound_ref'].split('/')[-1]
            sol['exchanges'][exchange_id] = {
                'x': flux['value'],
                'lower_bound': flux['lowerBound'],
                'upper_bound': flux['upperBound']}

        # Flux values for all of the reactions are reported separately.
        sol['reactions'] = dict()
        for flux in solution_data['FBAReactionVariables']:
            reaction_id = flux['modelreaction_ref'].split('/')[-1]
            sol['reactions'][reaction_id] = {
                'x': flux['value'],
                'lower_bound': flux['lowerBound'],
                'upper_bound': flux['upperBound']}

    solutions.sort(key=itemgetter('rundate'), reverse=True)  # Sort so last completed fba is first in list
    return solutions


def get_modelseed_gapfill_solutions(reference):
    """ Get the list of gap fill solutions for a ModelSEED model.

        @param reference: Workspace reference to model
        @return: List of gap fill solution data structures
    """

    solutions = None
    try:
        get_modelseed_model_stats(reference)  # Confirm model exists
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        solutions = ms_client.call('list_gapfill_solutions', {'model': reference})
    except PatricClient.ServerError as e:
        handle_server_error(e, [reference])

    # Convert the data about the gapfilled reactions from a list to a dictionary
    # keyed by reaction ID.
    for sol in solutions:
        if len(sol['solution_reactions']) > 1:
            warn('Gapfill solution {0} has {1} items in solution_reactions list'
                 .format(sol['id'], len(sol['solution_reactions'])))
        sol['reactions'] = dict()
        if len(sol['solution_reactions']) > 0:  # A gap fill solution can have no reactions
            for reaction in sol['solution_reactions'][0]:
                reaction_id = '{0}_{1}'.format(re.sub(modelseed_suffix_re, '', reaction['reaction'].split('/')[-1]),
                                               reaction['compartment'])
                sol['reactions'][reaction_id] = reaction
        del sol['solution_reactions']

    # Sort so last completed gap fill is first in list.
    solutions.sort(key=itemgetter('rundate'), reverse=True)
    return solutions


def get_modelseed_model_data(reference):
    """ Get the model data for a ModelSEED model.

        @param reference: Workspace reference to model
        @return: Model data dictionary
    """

    try:
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        return ms_client.call('get_model', {'model': reference, 'to': 1})
    except PatricClient.ServerError as e:
        handle_server_error(e, [reference])


def get_modelseed_model_stats(reference):
    """ Get the model statistics for a ModelSEED model.

        @param reference: Workspace reference to model
        @return: Model statistics dictionary
    """

    # The metadata for the model object has the data needed for the dictionary.
    metadata = get_workspace_object_meta(reference)

    # Build the model statistics dictionary.
    stats = dict()
    stats['fba_count'] = int(metadata[7]['fba_count'])
    stats['gapfilled_reactions'] = int(metadata[7]['gapfilled_reactions'])
    stats['gene_associated_reactions'] = int(metadata[7]['gene_associated_reactions'])
    stats['genome_ref'] = metadata[7]['genome_ref']
    stats['id'] = metadata[0]
    stats['integrated_gapfills'] = int(metadata[7]['integrated_gapfills'])
    stats['name'] = metadata[7]['name']
    stats['num_biomass_compounds'] = int(metadata[7]['num_biomass_compounds'])
    stats['num_biomasses'] = int(metadata[7]['num_biomasses'])
    stats['num_compartments'] = int(metadata[7]['num_compartments'])
    stats['num_compounds'] = int(metadata[7]['num_compounds'])
    stats['num_genes'] = int(metadata[7]['num_genes'])
    stats['num_reactions'] = int(metadata[7]['num_reactions'])
    stats['ref'] = metadata[7]['ref']
    stats['rundate'] = metadata[3]
    stats['source'] = metadata[7]['source']
    stats['source_id'] = metadata[7]['source_id']
    stats['template_ref'] = metadata[7]['template_ref']
    stats['type'] = metadata[7]['type']
    stats['unintegrated_gapfills'] = int(metadata[7]['unintegrated_gapfills'])

    return stats


def list_modelseed_models(root_path=None, print_output=False):
    """ List the ModelSEED models for the user.

        @param root_path: Root path to search for models
        @return: List of model what @todo Need details on structure
    """

    params = dict()
    if root_path is not None:
        params['path'] = root_path

    try:
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        output = ms_client.call('list_models', params)
    except PatricClient.ServerError as e:
        handle_server_error(e)
    if not print_output:
        return output
    for model in output:
        print('Model {0} for organism {1} with {2} reactions and {3} metabolites'
              .format(model['ref'], model['name'], model['num_reactions'], model['num_compounds']))
    return


def _convert_compartment(modelseed_id, format_type):
    """ Convert a compartment ID in ModelSEED source format to another format.

        @param modelseed_id: Compartment ID in ModelSEED source format
        @param format_type: Type of format to convert to
        @return: ID in specified format
    """

    if format_type == 'modelseed' or format_type == 'bigg':
        return modelseed_id[0]

    return modelseed_id


def _convert_suffix(modelseed_id, format_type):
    """ Convert a string with a compartment suffix from ModelSEED source format to another format.

        @param modelseed_id: String with compartment suffix in ModelSEED source format
        @param format_type: Type of format to convert to
        @return: ID in specified format
    """

    # Remove compartment index number for ModelSEED format. For example, rxn00001_c0 becomes rxn00001_c.
    # ModelSEED always uses compartment index 0 anyway.
    if format_type == 'modelseed':
        compartment = re.search(modelseed_suffix_re, modelseed_id).group(1)
        return re.sub(modelseed_suffix_re, '', modelseed_id) + '_{0}'.format(compartment)

    # Convert to BiGG type format. For example, rxn00001_c0 becomes rxn00001[c].
    elif format_type == 'bigg':
        match = re.search(modelseed_suffix_re, modelseed_id)
        compartment = match.group(1)
        return re.sub(modelseed_suffix_re, '', modelseed_id) + '[{0}]'.format(compartment)

    # No conversion is needed or format_type is unknown.
    return modelseed_id


def _add_metabolite(modelseed_compound, model, id_type):
    """ Create a COBRApy Metabolite object from a ModelSEED compound dictionary and add it to COBRApy model.

        @param modelseed_compound: Compound dictionary from ModelSEED model
        @param model: COBRApy Model object to add metabolite to
        @param id_type: Type of metabolite ID
        @return: True when metabolite is a duplicate
    """

    # Convert from ModelSEED format to COBRApy format.
    cobra_id = _convert_suffix(modelseed_compound['id'], id_type)
    cobra_compartment = _convert_compartment(modelseed_compound['modelcompartment_ref'].split('/')[-1], id_type)
    cobra_name = _convert_suffix(modelseed_compound['name'], id_type)

    # A ModelSEED model usually contains duplicate compounds. Confirm that the duplicate
    # compound is an exact duplicate and ignore it.
    if model.metabolites.has_id(cobra_id):
        metabolite = model.metabolites.get_by_id(cobra_id)
        if metabolite.name != cobra_name:
            warn('Duplicate ModelSEED compound ID {0} has different name, {1} != {2}'
                 .format(cobra_id, metabolite.name, cobra_name))
        if metabolite.formula != modelseed_compound['formula']:
            warn('Duplicate ModelSEED compound ID {0} has different formula, {1} != {2}'
                 .format(cobra_id, metabolite.formula, modelseed_compound['formula']))
        if metabolite.charge != modelseed_compound['charge']:
            warn('Duplicate ModelSEED compound ID {0} has different charge {1} != {2}'
                 .format(cobra_id, metabolite.charge, modelseed_compound['charge']))
        if metabolite.compartment != cobra_compartment:
            warn('Duplicate ModelSEED compound ID {0} has different compartment {1} != {2}'
                 .format(cobra_id, metabolite.compartment, cobra_compartment))
        return True

    # Create the Metabolite object and add it to the model.
    metabolite = Metabolite(id=cobra_id,
                            formula=modelseed_compound['formula'],
                            name=cobra_name,
                            charge=modelseed_compound['charge'],
                            compartment=cobra_compartment)
    model.add_metabolites([metabolite])

    return False


def _add_reaction(modelseed_reaction, model, id_type, likelihoods):
    """ Create a COBRApy Reaction object from a ModelSEED reaction dictionary and add it to COBRApy model.

        @param modelseed_reaction: Reaction dictionary from ModelSEED model
        @param model: COBRApy Model object to add reaction to
        @param id_type: Type of reaction ID
        @param likelihoods: Dictionary with reaction likelihoods from ModelSEED model
        @return: Nothing
    """

    # Set upper and lower bounds based directionality. Switch reverse reactions to forward reactions (ModelSEED
    # does this when exporting to SBML).
    reverse = 1.0
    if modelseed_reaction['direction'] == '=':
        lower_bound = -1000.0
        upper_bound = 1000.0
    elif modelseed_reaction['direction'] == '>':
        lower_bound = 0.0
        upper_bound = 1000.0
    elif modelseed_reaction['direction'] == '<':
        lower_bound = 0.0
        upper_bound = 1000.0
        reverse = -1.0
    else:
        warn('Reaction direction {0} assumed to be reversible for reaction {1}'
             .format(modelseed_reaction['direction'], modelseed_reaction['id']))
        lower_bound = -1000.0
        upper_bound = 1000.0

    # Create the Reaction object.
    reaction = Reaction(id=_convert_suffix(modelseed_reaction['id'], id_type),
                        name=re.sub(modelseed_suffix_re, '', modelseed_reaction['name']),
                        lower_bound=lower_bound,
                        upper_bound=upper_bound)

    # Create dictionary of metabolites and add them to the reaction.
    metabolites = dict()
    for reagent in modelseed_reaction['modelReactionReagents']:
        cobra_metabolite_id = _convert_suffix(reagent['modelcompound_ref'].split('/')[-1], id_type)
        metabolite = model.metabolites.get_by_id(cobra_metabolite_id)
        metabolites[metabolite] = float(reagent['coefficient']) * reverse
    reaction.add_metabolites(metabolites)

    # If there are proteins associated with the reaction, build the gene reaction rule and
    # add corresponding Gene objects to the model.
    if len(modelseed_reaction['modelReactionProteins']) > 0:
        # Build a list of proteins associated with the reaction.
        protein_list = list()
        for protein in modelseed_reaction['modelReactionProteins']:
            # Spontaneous and universal reactions can have an entry in the list of proteins
            # that does not have any protein subunits.
            if len(protein['modelReactionProteinSubunits']) == 0:
                continue

            # Build a list of protein subunits in this protein.
            subunit_list = list()
            for subunit in protein['modelReactionProteinSubunits']:
                # A protein with multiple subunits can have a subunit that is not linked
                # to a feature in the genome of the organism.
                if len(subunit['feature_refs']) == 0:
                    continue

                # Build a list of features in this protein subunit.
                feature_list = list()
                for feature in subunit['feature_refs']:
                    # Extract the gene ID from the reference to the feature in the genome and
                    # remove the "fig|" prefix.
                    gene_id = re.sub(patric_gene_prefix_re, '', feature.split('/')[-1])
                    if not model.genes.has_id(gene_id):
                        gene = Gene(gene_id, subunit['role'])  # Use the role name as the Gene name
                        model.genes.append(gene)
                    feature_list.append(gene_id)

                #  Join multiple features using an OR relationship.
                if len(feature_list) > 1:
                    subunit_list.append('( {0} )'.format(' or '.join(feature_list)))
                else:
                    subunit_list.append(feature_list[0])

            # Join multiple protein subunits using an AND relationship.
            if len(subunit_list) > 0:
                if len(subunit_list) > 1:
                    protein_list.append('( {0} )'.format(' and '.join(subunit_list)))
                else:
                    protein_list.append(subunit_list[0])

        # If there is an association to a feature, add the rule to the reaction.
        if len(protein_list) > 0:
            # Join multiple proteins using an OR relationship.
            if len(protein_list) > 1:
                gpr_rule = '( {0} )'.format(' or '.join(protein_list))
            else:
                gpr_rule = protein_list[0]

            reaction.gene_reaction_rule = gpr_rule

    # Add a note with gapfill details. ModelSEED gapfill data is a dictionary where the key is the
    # gapfill solution ID and the value indicates if the reaction was added or reversed and the
    # direction of the reaction. For example: {u'gf.0': u'added:>'}
    if len(modelseed_reaction['gapfill_data']) > 0:
        reaction.notes['gapfill_data'] = modelseed_reaction['gapfill_data']

    # Add a note with likelihood value if available.
    if modelseed_reaction['id'] in likelihoods:
        reaction.notes['likelihood'] = likelihoods[modelseed_reaction['id']]

    # Finally, add the reaction to the model.
    model.add_reaction(reaction)

    return


def create_cobra_model_from_modelseed_model(reference, id_type='modelseed', validate=False):
    """ Create a COBRA model from a ModelSEED model.

        @param reference: Workspace reference to model
        @param id_type: Type of IDs ('modelseed' for _c or 'bigg' for '[c])
        @param validate: When True, check for common problems
        @return COBRApy Model object
    """

    # Validate the id_type parameter.
    if id_type == 'modelseed':
        cytosol_suffix = '_c'
    elif id_type == 'bigg':
        cytosol_suffix = '[c]'
    else:
        raise ValueError('id_type {0} is not supported'.format(id_type))

    # Get the ModelSEED model data.
    data = get_modelseed_model_data(reference)

    # Get the workspace object with the likelihoods and put the likelihood values in a dictionary
    # keyed by reaction ID. Calculating likelihoods is optional so the object may not exist.
    try:
        likelihood_data = get_workspace_object_data(join(reference, 'rxnprobs'))
        # for reaction_data in likelihood_data['reaction_probabilities']:
        #     likelihoods[reaction_data[0]] = reaction_data[1]
        likelihoods = {r[0]: r[1] for r in likelihood_data['reaction_probabilities']}
    except PatricClient.ObjectNotFoundError:
        likelihoods = dict()

    # Create a new COBRApy Model object.
    model = Model(data['id'], name=data['name'])

    # Add compartments to the COBRApy model.
    for index in range(len(data['modelcompartments'])):
        modelseed_compartment = data['modelcompartments'][index]
        cobra_id = _convert_compartment(modelseed_compartment['id'], id_type)
        model.compartments[cobra_id] = modelseed_compartment['label'][:-2]  # Strip _0 suffix from label

    # Create Metabolite objects for all of the compounds in the ModelSEED model.
    num_duplicates = 0
    for index in range(len(data['modelcompounds'])):
        duplicate = _add_metabolite(data['modelcompounds'][index], model, id_type=id_type)
        if duplicate:
            num_duplicates += 1

    # Report the number of duplicate metabolites.
    if num_duplicates > 0:
        warn('{0} duplicate metabolites were removed from model {1} of {2}'
             .format(num_duplicates, model.id, model.name))

    # Add all of the reactions to the COBRApy model.
    for index in range(len(data['modelreactions'])):
        _add_reaction(data['modelreactions'][index], model, id_type, likelihoods)

    # Add exchange reactions for metabolites in extracellular compartment.
    for index in range(len(model.metabolites)):
        metabolite = model.metabolites[index]
        if metabolite.compartment.startswith('e'):
            # Single reactant metabolite makes a system boundary reaction.
            reaction = Reaction(id='EX_' + metabolite.id,
                                name=metabolite.name + ' exchange',
                                lower_bound=-1000.0,
                                upper_bound=1000.0)  # @todo Should the upper bound be 0?
            reaction.add_metabolites({metabolite: -1.0})
            model.add_reaction(reaction)

    # A ModelSEED model must have an exchange reaction for the special biomass metabolite.
    metabolite = model.metabolites.get_by_id('cpd11416'+cytosol_suffix)
    reaction = Reaction(id='EX_' + metabolite.id,
                        name=metabolite.name,
                        lower_bound=-1000.0,
                        upper_bound=1000.0)
    reaction.add_metabolites({metabolite: -1.0})
    model.add_reaction(reaction)

    # Note that when ModelSEED exports to SBML, it includes exchange reactions for
    # cpd15302_c0 "Glycogen(n-1)" and cpd08636_c0 "4-5-dihydroxy-2-3-pentanedione".
    # No idea why exchange reactions for metabolites in the cytosol compartment are
    # added to the SBML file.
    # metabolite = model.metabolites.get_by_id('cpd08636'+cytosol_suffix)
    # reaction = Reaction(id='EX_' + metabolite.id,
    #                     name=metabolite.name,
    #                     lower_bound=-1000.0,
    #                     upper_bound=1000.0)
    # reaction.add_metabolites({metabolite: -1.0})
    # model.add_reaction(reaction)
    #
    # metabolite = model.metabolites.get_by_id('cpd15302'+cytosol_suffix)
    # reaction = Reaction(id='EX_' + metabolite.id,
    #                     name=metabolite.name,
    #                     lower_bound=-1000.0,
    #                     upper_bound=1000.0)
    # reaction.add_metabolites({metabolite: -1.0})
    # model.add_reaction(reaction)

    # Add the biomass reactions to the COBRApy model. ModelSEED models can have more than one biomass reaction
    # but the model does not identify which one to use as the objective so always use the first one.
    if len(data['biomasses']) > 1:
        warn('Found {0} biomass reactions and selected {1} as the objective'
             .format(len(data['biomasses']), data['biomasses'][0]['id']))
    for index in range(len(data['biomasses'])):
        biomass = data['biomasses'][index]
        biomass_metabolites = dict()
        for biomass_compound in biomass['biomasscompounds']:
            cobra_id = _convert_suffix(biomass_compound['modelcompound_ref'].split('/')[-1], id_type)
            metabolite = model.metabolites.get_by_id(cobra_id)
            biomass_metabolites[metabolite] = biomass_compound['coefficient']
        reaction = Reaction(id=biomass['id'],
                            name=biomass['name'],
                            lower_bound=0.0,
                            upper_bound=1000.0)
        reaction.add_metabolites(biomass_metabolites)
        if index == 0:
            reaction.objective_coefficient = 1.
        model.add_reaction(reaction)

    # If requested, validate the COBRApy model.
    if validate:
        # See if all of the reactions are mass balanced.
        num_unbalanced = 0
        for r in model.reactions:
            if not r.id.startswith('EX_'):  # Skip exchange reactions
                unbalanced = r.check_mass_balance()
                if len(unbalanced) > 0:
                    warn('Reaction {0} is unbalanced because {1}\n    {2}'
                         .format(r.id, unbalanced, r.build_reaction_string(use_metabolite_names=True)))
                    num_unbalanced += 1
        if num_unbalanced > 0:
            warn('Model {0} has {1} unbalanced reactions'.format(model.id, num_unbalanced))

    return model


def optimize_modelseed_model(reference, media_reference=None):
    """ Run flux balance analysis on a ModelSEED model.

        @param reference: Reference to model
        @param media_reference: Reference to media to optimize on (default is complete media)
        @return: Flux balance analysis data structure
    """

    # Set input parameters for method.
    params = dict()
    params['model'] = reference
    params['predict_essentiality'] = 1
    if media_reference is not None:
        params['media'] = media_reference

    # Run the server method.
    try:
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        job_id = ms_client.call('FluxBalanceAnalysis', params)
        _wait_for_job(job_id)
    except PatricClient.ServerError as e:
        references = [reference]
        if media_reference is not None:
            references.append(media_reference)
        handle_server_error(e, references)

    # The completed job does not have the reference to the fba object that
    # was just created so get the list of solutions. Last completed
    # solution is first in the list.
    solutions = get_modelseed_fba_solutions(reference)
    return solutions[0]


def reconstruct_modelseed_model(genome_id, source='PATRIC', template_reference=None, likelihood=False, name=None):
    """ Reconstruct a draft ModelSEED model for an organism.

        @param genome_id: Genome ID or workspace reference to genome
        @param source: Source of genome
        @param template_reference: Reference to template model
        @param likelihood: True to generate reaction likelihoods
        @param name: Name of output model
        @return: Model statistics dictionary
    """

    # Set input parameters for method.
    params = dict()
    if source == 'PATRIC':
        params['genome'] = 'PATRIC:' + genome_id
    elif source == 'RAST':
        params['genome'] = 'RAST:' + genome_id
    elif source == 'workspace':
        params['genome'] = genome_id
    else:
        raise ValueError('Source type {0} is not supported'.format(source))
    if name is None:
        name = genome_id
    params['output_file'] = name
    if template_reference is not None:
        params['templatemodel'] = template_reference
    if likelihood:
        params['probanno'] = 1
    else:
        params['probanno'] = 0
    params['gapfill'] = 0
    params['predict_essentiality'] = 0

    # Run the server method.
    try:
        ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
        job_id = ms_client.call('ModelReconstruction', params)
    except PatricClient.ServerError as e:
        references = None
        if template_reference is not None:
            references = [template_reference]
        handle_server_error(e, references)

    # The task structure has the workspace where the model is stored but not the name of the model.
    task = _wait_for_job(job_id)
    reference = task['workspace'] + name

    # Get the model statistics for the model.
    stats = get_modelseed_model_stats(reference)
    if stats['num_genes'] == 0:  # ModelSEED does not return an error if the genome ID is invalid
        warn('Model for genome ID {0} has no genes, verify genome ID is valid'.format(genome_id))
    return stats


def _wait_for_job(jobid):
    """ Wait for a job submitted to the ModelSEED app service to end.

        @param jobid: ID of submitted job
        @return Task structure with status of job
    """

    task = None
    done = False
    ms_client = PatricClient.PatricClient(modelseed_url, 'ProbModelSEED')
    while not done:
        jobs = ms_client.call('CheckJobs', {'jobs': [jobid]})
        if jobid in jobs:
            task = jobs[jobid]
            if task['status'] == 'failed':
                raise PatricClient.ServerError(task['errors'])
            elif task['status'] == 'completed':
                done = True
            else:
                sleep(3)
        else:
            raise PatricClient.JobError('Job {0} was not found'.format(jobid))
    return task


def get_workspace_object_meta(reference):
    """ Get the metadata for an object.

        @param reference: Reference to object in workspace
        @return ObjectMeta tuple with metadata
    """

    try:
        # The output from get() is a list of tuples.  When asking for metadata only,
        # the list entry is a tuple with only one element.
        ws_client = PatricClient.PatricClient(workspace_url, 'Workspace')
        metadata_list = ws_client.call('get', {'objects': [reference], 'metadata_only': 1})
        return metadata_list[0][0]
    except PatricClient.ServerError as e:
        handle_server_error(e, [reference])


def get_workspace_object_data(reference, json_data=True):
    """ Get the data for an object.

        @param reference: Reference to object in workspace
        @param json_data: When True, return data in JSON format
        @return Object data
    """

    data = None
    try:
        # The output from get() is a list of tuples. The first element in the
        # tuple is the metadata which has the url to the Shock node when the
        # data is stored in Shock. Otherwise the data is available in the second
        # element of the tuple.
        ws_client = PatricClient.PatricClient(workspace_url, 'Workspace')
        object_list = ws_client.call('get', {'objects': [reference]})
        if len(object_list[0][0][11]) > 0:
            data = PatricClient.shock_download(object_list[0][0][11], ws_client.token)
        else:
            data = object_list[0][1]
    except Exception as e:
        handle_server_error(e, [reference])

    if json_data:
        return json.loads(data)
    return data


def list_workspace_objects(folder, sort_key='folder', print_output=False):
    """ List the objects in the specified workspace folder.

        @param folder: Path to workspace folder
        @param sort_key: Name of field to use as sort key for output
        @param print_output: When True, print formatted output
        @return: List of object data for objects in folder or None if printed output
    """

    # Get the list of objects in the specified folder.
    try:
        ws_client = PatricClient.PatricClient(workspace_url, 'Workspace')
        output = ws_client.call('ls', {'paths': [folder], 'recursive': 1})
    except PatricClient.ServerError as e:
        handle_server_error(e, [folder])

    # See if the folder is in the returned data.
    if folder not in output:
        if print_output:
            print('No data for folder {0} was available'.format(folder))
        return None

    # Sort the objects by the specified key.
    reverse = False
    if sort_key == 'name':
        key = 0
    elif sort_key == 'folder':
        key = 2
    elif sort_key == 'date':
        key = 3
        reverse = True
    elif sort_key == 'type':
        key = 1
    else:
        raise ValueError('Sort key {0} is not supported'.format(sort_key))
    output[folder].sort(key=itemgetter(key), reverse=reverse)

    # Print details on the objects.
    if print_output:
        print('Contents of {0}:'.format(folder))
        for object_data in output[folder]:
            otype = '-'
            if object_data[1] == 'folder' or object_data[8]['is_folder']:
                otype = 'd'
            print('{0}{1}{2} {3:10}\t{4:>10}\t{5}\t{6:12}\t{7}{8}'
                  .format(otype, object_data[9], object_data[10], object_data[5], object_data[6],
                          object_data[3], object_data[1], object_data[2], object_data[0]))
        return None

    return output[folder]


def handle_server_error(e, references=None):
    """ Handle an error returned by a PATRIC service server.

        @param e: Exception returned by server
        @param references: List of workspace references for server method
        @return: Nothing
        @raise: ObjectNotFoundError
    """

    # Map a generic server error to a specific exception.
    if isinstance(e, PatricClient.ServerError):
        # Errors returned by ModelSEED server do not indicate which workspace reference had
        # a problem. One of the references in the list is the one that caused the problem.
        reference_string = ''
        if references is not None:
            reference_string = '"{0}"'.format('" or "'.join(references))
        if 'Object not found!' in e.message or 'Owner not specified in deletion!' in e.message:
            msg = 'An object was not found in workspace: {0}'.format(reference_string)
            raise PatricClient.ObjectNotFoundError(msg, e.data)

        if 'does not include at least a top level directory!' in e.message:
            msg = 'An object reference is missing a top level directory: {0}'.format(reference_string)
            raise PatricClient.ObjectNotFoundError(msg, e.data)

        if 'Path does not point to folder or object:' in e.message:
            msg = 'An object was not found in workspace: {0}'.format(reference_string)
            raise PatricClient.ObjectNotFoundError(msg, e.data)

        if 'User lacks permission to / for requested action!' in e.message:
            msg = 'User does not have permission to a directory in an object reference: {0}'.format(reference_string)
            raise PatricClient.ObjectNotFoundError(msg, e.data)

        if 'is not a valid object path!' in e.message:
            msg = 'An object reference is not a valid path: {0}'.format(reference_string)
            raise PatricClient.ObjectNotFoundError(msg, e.data)

        if 'No gapfilling needed on specified condition' in e.message:
            raise PatricClient.DuplicateGapfillError('Gap fill solution already available for specified media')

        if 'does not match specified type' in e.message:
            raise PatricClient.ObjectTypeError(e.message, e.data)

    # Raise the exception again.
    raise e
