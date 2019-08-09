from collections import defaultdict
from functools import reduce
import networkx as nx
from transform.pdc.pdc_harvester import load_all


class PDCTransformer():

    def __init__(self):
        # fetch data
        (self.all_cases, self.all_programs, self.all_studies) = load_all()
        print('all_cases', len(self.all_cases))
        print('all_programs', len(self.all_programs))
        print('all_studies', len(self.all_studies))
        print('all_files', reduce(lambda s, x: s + len(x.files), self.all_studies, 0))
        # assert that all cases exist within known projects
        program_project = set()
        for program in self.all_programs:
            for project in program.projects:
                program_project.add(project.project_submitter_id)
        self.case_project = defaultdict(int)
        for c in self.all_cases:
            self.case_project[c.project_submitter_id] += 1
        assert self.case_project.keys() == program_project

        # create summary by project
        def counter_dict():
            return defaultdict(int)

        self.sample_gdc_projects = defaultdict(counter_dict)
        for c in self.all_cases:
            for s in c.samples:
                if s.gdc_project_id in ['', None]:
                    s.gdc_project_id = 'N/A'
                self.sample_gdc_projects[c.project_submitter_id][s.gdc_project_id] += 1

    def graph(self):
        """Creates graph"""
        G = nx.MultiDiGraph()

        G.add_nodes_from(set([case.project_submitter_id for case in self.all_cases]), label='Project')
        G.add_nodes_from(set([case.case_submitter_id for case in self.all_cases]), label='Case')

        G.add_edges_from([(case.project_submitter_id, case.case_submitter_id) for case in self.all_cases], label='cases')

        def case_samples():
            for case in self.all_cases:
                for s in case.samples:
                    yield (case.case_submitter_id, s.sample_id)

        def case_diagnoses():
            for case in self.all_cases:
                for d in case.diagnoses:
                    yield (case.case_submitter_id, d.diagnosis_id)

        def case_demographics():
            for case in self.all_cases:
                for d in case.demographics:
                    yield (case.case_submitter_id, d.demographic_id)

        def sample_aliquots():
            for case in self.all_cases:
                for s in case.samples:
                    for a in s.aliquots:
                        yield (s.sample_id, a.aliquot_id)

        def aliquot_studies():
            for s in self.all_studies:
                for bs in s.biospecimens:
                    yield (bs.aliquot_id, s.study_submitter_id)

        def study_files():
            for s in self.all_studies:
                for f in s.files:
                    yield (s.study_submitter_id, f.file_id)

        G.add_nodes_from([case_sample[1] for case_sample in case_samples()], label='Sample')
        G.add_edges_from([case_sample for case_sample in case_samples()], label='samples')

        G.add_nodes_from([sample_aliquot[1] for sample_aliquot in sample_aliquots()], label='Aliquot')
        G.add_edges_from([sample_aliquot for sample_aliquot in sample_aliquots()], label='aliquots')

        G.add_nodes_from([case_diagnosis[1] for case_diagnosis in case_diagnoses()], label='Diagnosis')
        G.add_edges_from([case_diagnosis for case_diagnosis in case_diagnoses()], label='diagnoses')

        G.add_nodes_from([case_demographic[1] for case_demographic in case_demographics()], label='Demographic')
        G.add_edges_from([case_demographic for case_demographic in case_demographics()], label='demographics')

        G.add_nodes_from([aliquot_study[1] for aliquot_study in aliquot_studies()], label='Study')
        G.add_edges_from([aliquot_study for aliquot_study in aliquot_studies()], label='studies')

        G.add_nodes_from([study_file[1] for study_file in study_files()], label='File')
        G.add_edges_from([study_file for study_file in study_files()], label='files')

        return G
